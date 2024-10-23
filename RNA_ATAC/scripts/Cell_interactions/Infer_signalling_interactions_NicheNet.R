###############################################
##                                           ##
##  Infer signalling interactions (NicheNet) ##
##                                           ##
###############################################

source("/data/homes/louisc/Project_Babraham/RNA_ATAC/scripts/Settings.R")
library(nichenetr)
library(openxlsx)

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--sce',       type="character",                help='SingleCellExperiment')
p$add_argument('--sce_mc',       type="character",                help='SingleCellExperiment on metacell level')
p$add_argument('--metadata',       type="character",                help='Metadata metacells')
p$add_argument('--CellPhoneDB_data',       type="character",                help='CellPhoneDB data')
p$add_argument('--ncores',  type="integer",            default=4,      help='Amount of cores to use')
p$add_argument('--min_expr',  type="integer",            default=4,      help='Amount of cores to use')
p$add_argument('--min_coef',  type="double",            default=0.25,      help='Minimal regression coefficient')
p$add_argument('--max_pval',  type="double",            default=0.1,      help='Maximal p value')
p$add_argument('--DEG_overview',  type="character", help='Method for GRN building')
p$add_argument('--GRN_method',  type="character", help='Method for GRN building')
p$add_argument('--merged',  type="character", help='Merge cluster 4 and 5')
p$add_argument('--markers_TF',  type="character" ,     help='TF marker file')
p$add_argument('--plot_correlations',  type="logical" , default=FALSE,     help='Logical indicating if correlations have to be plotted')
p$add_argument('--outdir',       type="character",                help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

# # START TEST
# args <- list()
# args$sce <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/SingleCellExperiment_pseudobulk.rds"
# args$sce_mc <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/metacells/SingleCellExperiment_metacells_nodiff.rds"
# args$metadata <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/clustering/PeakMatrix/sample_metadata_nodiff_after_clustering.txt.gz"
# args$CellPhoneDB_data <- "/data/homes/louisc/Project_Babraham/pathways/CellPhoneDB_analysis/CellPhone_res_stat_LIGANDS_RECEPTORS.xlsx"
# args$ncores <- 10
# args$min_expr <- 4
# args$DEG_overview <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/RNA/DE_res/DEG_overview.txt"
# args$GRN_method <- "alt2"
# args$markers_TF <- "/data/homes/louisc/Project_Babraham/RNA_ATAC//pseudobulk/cluster/RNA/marker_TFs/TF.markers.clusters.txt"
# args$plot_correlations <- FALSE
# args$outdir <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/NicheNet/"
# # END TEST

# print(args)

###############################
## Load additional resources ##
###############################

ligand_target_matrix <- readRDS("/data/homes/louisc/Project_Babraham/RNA_ATAC/NicheNet/ligand_target_matrix.rds")

TF_markers <- read.table(args$markers_TF,header=T,sep="\t")$Gene

markers_TF <- NULL
for (i in 1:5){
  marker_tf_i <- read.table(sprintf("/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/RNA/marker_TFs/Marker_TF_cluster%s.txt",i),header=T,sep="\t")
  markers_TF <- rbind(markers_TF,cbind(marker_tf_i$Gene,rep(i,length(marker_tf_i$Gene))))
}
markers_TF <- data.frame(markers_TF)
colnames(markers_TF) <- c("TF","celltype")
markers_TF_undupl <- markers_TF[!duplicated(markers_TF$TF),]

##############
## Load sce ##
##############

sce <- readRDS(args$sce)

##########################
## CellPhoneDB genesets ##
##########################

dat_cpdb_ligand <- read.xlsx(args$CellPhoneDB_data,2)
colnames(dat_cpdb_ligand) <- c("cluster","goi")
dat_cpdb_receptor <- read.xlsx(args$CellPhoneDB_data,3)
colnames(dat_cpdb_receptor) <- c("cluster","goi")

dat_cpdb <- rbind(dat_cpdb_ligand,dat_cpdb_receptor)
dat_cpdb <- data.frame(dat_cpdb)
dat_cpdb$cluster <- gsub("\\|[1-5]$","",dat_cpdb$cluster)
print(dim(dat_cpdb))
print(head(dat_cpdb))

sum(duplicated(paste(dat_cpdb$goi,dat_cpdb$cluster,sep="_")))
dat_cpdb[duplicated(paste(dat_cpdb$goi,dat_cpdb$cluster,sep="_")),]$goi
dat_cpdb <- dat_cpdb[!duplicated(paste(dat_cpdb$goi,dat_cpdb$cluster,sep="_")),]
print(dim(dat_cpdb))
print(head(dat_cpdb))

list_cpdb <- list(list(dat_cpdb$goi[dat_cpdb$cluster==1],1),
                  list(dat_cpdb$goi[dat_cpdb$cluster==2],2),
                  list(dat_cpdb$goi[dat_cpdb$cluster==3],3),
                  list(dat_cpdb$goi[dat_cpdb$cluster==4],4),
                  list(dat_cpdb$goi[dat_cpdb$cluster==5],5),
                  list(unique(dat_cpdb$goi),"X"))
names(list_cpdb) <- c("CellPhoneDB_clusters1","CellPhoneDB_clusters2","CellPhoneDB_clusters3",
                      "CellPhoneDB_clusters4","CellPhoneDB_clusters5","CellPhoneDB_clustersX")

genesets <- list_cpdb

######################################
## Signalling molecules CellPhoneDB ##
######################################

dat_all <- read.table("/data/homes/louisc/Project_Babraham/RNA_ATAC/NicheNet/Heatmaps/CellPhoneDB_ligand_to_Expr_targets.txt",header=T,row.names=1,sep="\t")
print(dim(dat_all))
print(sum(dat_all>=0.1))

reshape_data <- function(i){
    if (i%%1000==0){print(i)}
    dat_reshaped <- cbind(rownames(dat_all),rep(colnames(dat_all)[i],nrow(dat_all)),dat_all[,i])
    return(dat_reshaped)
}
list_dat_reshaped <- mclapply(1:ncol(dat_all),reshape_data,mc.cores=args$ncores,mc.preschedule=T)
dat_reshaped <- do.call(rbind,list_dat_reshaped)
dat_reshaped <- data.frame(dat_reshaped)
colnames(dat_reshaped) <- c("from","to","weight")
dat_reshaped$weight <- as.numeric(dat_reshaped$weight)
print(head(dat_reshaped))
print(dim(dat_reshaped))

###################
## Functions GRN ##
###################

source("/data/homes/louisc/Project_Babraham/RNA_ATAC/scripts/GRN/build_GRN_functions.R")

##########################
## Load sce (metacells) ##
##########################

sce_mc <- readRDS(args$sce_mc)

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE]

sample_metadata <- data.frame(sample_metadata)
rownames(sample_metadata) <- sample_metadata$sample
sample_metadata_metacells <- sample_metadata[colnames(sce_mc),]

colData(sce_mc)$cluster <- sample_metadata_metacells$cluster

#######################################################
## Filter signalling connections for expressed genes ##
#######################################################

rna_tf.mtx <- logcounts(sce_mc)
rownames(rna_tf.mtx) <- toupper(rownames(rna_tf.mtx))
rna_targets.mtx <- logcounts(sce_mc)

dat_filtered <- dat_reshaped[dat_reshaped$weight>=0.1,]
head(dat_filtered)
dim(dat_filtered)

dat_filtered <- dat_filtered[dat_filtered$from%in%rownames(sce_mc),]
dat_filtered <- dat_filtered[dat_filtered$to%in%rownames(sce_mc),]
head(dat_filtered)
dim(dat_filtered)

genes <- intersect(unique(c(dat_filtered$to)),rownames(sce_mc))
dat_cpdb_both <- read.xlsx(args$CellPhoneDB_data,1)
dat_cpdb_both$gene_a[is.na(dat_cpdb_both$gene_a)] <- ""
dat_cpdb_both$gene_b[is.na(dat_cpdb_both$gene_b)] <- ""

dat_filtered$receptor <- NA
for (i in 1:nrow(dat_filtered)){
    l_i <- dat_filtered$from[i]
    if ((sum(dat_cpdb_both$gene_a==l_i)+sum(dat_cpdb_both$gene_b==l_i))==1){
        dat_filtered$receptor[i] <- c(dat_cpdb_both$gene_b[dat_cpdb_both$gene_a==l_i],dat_cpdb_both$gene_a[dat_cpdb_both$gene_b==l_i])
    } else if (sum(dat_cpdb_both$gene_a==l_i)>1){
        dat_filtered$receptor[i] <- paste( c(dat_cpdb_both$gene_b[dat_cpdb_both$gene_a==l_i],dat_cpdb_both$gene_a[dat_cpdb_both$gene_b==l_i]),collapse=",")
    }
}
colnames(dat_filtered) <- c("ligand","target","weight","receptor")
dat_filtered <- separate_rows(dat_filtered,receptor,sep=",")
dat_filtered$receptor_ori <- dat_filtered$receptor

########################
## Manual adaptations ##
########################

# for CXCL12 - CXCR4
dat_filtered$receptor[dat_filtered$ligand=="CXCL12"] <- "CXCR4"
# for HGF - MET
dat_filtered$receptor[dat_filtered$ligand=="HGF"] <- "MET"
# for NGF - NTRK1 and NGFR
lines_NGF <-dat_filtered[dat_filtered$ligand=="NGF",]
dat_filtered$receptor[dat_filtered$ligand=="NGF"] <- "NTRK1"
dat_filtered <- rbind(dat_filtered,lines_NGF) 
dat_filtered$receptor[dat_filtered$ligand=="NGF" & is.na(dat_filtered$receptor)] <- "NGFR"
# for PDGFB - PDGFRA and PDGFRB
lines_PDGFB <- dat_filtered[dat_filtered$ligand=="PDGFB",]
dat_filtered$receptor[dat_filtered$ligand=="PDGFB"] <- "PDGFRB"
dat_filtered <- rbind(dat_filtered,lines_PDGFB) 
dat_filtered$receptor[dat_filtered$ligand=="PDGFB" & is.na(dat_filtered$receptor)] <- "PDGFRA"
# TGFbeta_receptor1 to TGBR1
dat_filtered$receptor[dat_filtered$receptor=="TGFbeta_receptor1"] <- "TGFBR1"

# Split complexes
dat_filtered <- separate_rows(dat_filtered,receptor,sep="_")
dat_filtered <- dat_filtered[,c("ligand","receptor","receptor_ori","target","weight")]

dim(dat_filtered)
head(dat_filtered)
dat_filtered$receptor[(!dat_filtered$receptor%in%rownames(sce_mc)) | is.na(dat_filtered$receptor)] <- NA
dat_filtered <- dat_filtered[!duplicated(paste(dat_filtered$ligand,dat_filtered$receptor,dat_filtered$target,sep="_")),]
head(dat_filtered)
dim(dat_filtered)

########################
## Cluster annotation ##
########################

DEG_overview <- read.table(file=args$DEG_overview, header=T, sep="\t", quote="")

# max gene expression
cluster_max_expression <- cbind(DEG_overview[,1],apply(DEG_overview[,2:6],1,which.max))
cluster_max_expression <- data.frame(cluster_max_expression)
colnames(cluster_max_expression) <- c("gene","cluster")
rownames(cluster_max_expression) <- cluster_max_expression$gene

# annotate markers 
marker_gene_annot <- cbind(DEG_overview[,1],lapply(apply(DEG_overview[,grepl("clust[1-5]_spec",colnames(DEG_overview))],1,which),paste,collapse=","))
marker_gene_annot <- data.frame(marker_gene_annot)
colnames(marker_gene_annot) <- c("gene","cluster")
head(marker_gene_annot[marker_gene_annot$cluster!="",])
sum(marker_gene_annot$cluster!="")
rownames(marker_gene_annot) <- marker_gene_annot$gene
sce_mc_marker_gene_annot <- marker_gene_annot
marker_gene_annot$cluster <- gsub(",[0-9]$","",marker_gene_annot$cluster)

gene_cluster_annot <- data.frame(gene=unique(c(dat_filtered$ligand,dat_filtered$receptor,dat_filtered$target)),
                                cluster=NA)                                   
gene_cluster_annot$cluster[gene_cluster_annot$gene%in%marker_gene_annot$gene] <- marker_gene_annot[gene_cluster_annot$gene[gene_cluster_annot$gene%in%marker_gene_annot$gene],]$cluster
gene_cluster_annot$cluster[gene_cluster_annot$cluster==""] <- cluster_max_expression[gene_cluster_annot$gene[gene_cluster_annot$cluster==""],]$cluster
head(gene_cluster_annot)

############################
## Merge cluster 5 to 4  ##
############################

if(args$merged=="merged"){
  gene_cluster_annot$cluster[gene_cluster_annot$cluster==5] <- 4
}

#######################
## Select GRN method ##
#######################

print(paste0("Method for GRN coef calculation: ",args$GRN_method))

print(paste0("myfun_",args$GRN_method))
myfun_used <- get(paste0("myfun_",args$GRN_method))

print(head(myfun_used))

#######################################
## Build GRN for ligands & receptors ##
#######################################

sign_senders <- c(dat_filtered$ligand,dat_filtered$receptor)

comparisons <- c("ligand","receptor")
GRN_coef.dt_total <- dat_filtered
for (comp in comparisons){
    if (comp == "ligand"){
        tmp <- dat_filtered[,c(comp,"target")]
        tmp <- tmp[!duplicated(paste(unlist(tmp[,1]),unlist(tmp[,2]),sep="_")),]
    } else if (comp == "receptor"){
        tmp <- dat_filtered[,c(comp,"target")]
        tmp <- tmp[!is.na(tmp[,1]),]
        tmp <- tmp[!duplicated(paste(unlist(tmp[,1]),unlist(tmp[,2]),sep="_")),]
    }

    tmp <- tmp[,1:2]
    colnames(tmp) <- c("tf","gene")
    tmp <- data.table(tmp)
    print(tmp)
    print(dim(tmp))

    TFs <-  intersect(unique(tmp$tf),rownames(sce_mc))
    print(head(genes))
    genes <- intersect(unique(tmp$gene),rownames(sce_mc))
    print(head(genes))

    print(length(genes))
    list_results <- mclapply(genes,myfun_used,mc.cores=args$ncores,mc.preschedule=T)
    print(head(list_results))
    GRN_coef.dt <- do.call(rbind,list_results)
    GRN_coef.dt$beta <- as.numeric(GRN_coef.dt$beta)
    GRN_coef.dt$pvalue <- as.numeric(GRN_coef.dt$pvalue)

    GRN_coef.dt_total[paste0(comp,"_beta")] <- NA
    GRN_coef.dt_total[paste0(comp,"_pvalue")] <- NA

    mask_comp <- (paste(unlist(GRN_coef.dt_total[,comp]),GRN_coef.dt_total$target,sep="_")%in%paste(GRN_coef.dt$tf,GRN_coef.dt$gene,sep="_")) & 
                    (!duplicated(paste(unlist(GRN_coef.dt_total[,comp]),GRN_coef.dt_total$target,sep="_")))
    if (sum(mask_comp)!=nrow(GRN_coef.dt)){warning("wrong dimensions")}
    GRN_coef.dt_total[mask_comp,paste0(comp,"_beta")] <- GRN_coef.dt$beta
    GRN_coef.dt_total[mask_comp,paste0(comp,"_pvalue")] <- GRN_coef.dt$pvalue
}
GRN_coef.dt_total <- data.frame(GRN_coef.dt_total)

#######################################################
## Combine GRN for ligands-receptors-targets combo's ##
#######################################################

# duplicate values for ligand-tf pairs
list_mult_lt <- paste(GRN_coef.dt_total$ligand,GRN_coef.dt_total$target,sep="_")
list_mult_lt <- unique(list_mult_lt[duplicated(list_mult_lt)])
for (i in list_mult_lt){
    mask_mult_lt <- paste(GRN_coef.dt_total$ligand,GRN_coef.dt_total$target,sep="_")==i
    GRN_coef.dt_total[mask_mult_lt,c("ligand_beta","ligand_pvalue")] <- t(matrix(rep(c(GRN_coef.dt_total$ligand_beta[which(mask_mult_lt)[1]],GRN_coef.dt_total$ligand_pvalue[which(mask_mult_lt)[1]]),sum(mask_mult_lt)),nrow=2))
}

# duplicate values for receptor-tf pairs
list_mult_rc <- paste(GRN_coef.dt_total$receptor,GRN_coef.dt_total$target,sep="_")
list_mult_rc <- unique(list_mult_rc[duplicated(list_mult_rc) & !is.na(GRN_coef.dt_total$receptor)])
for (i in list_mult_rc){
    mask_mult_rc <- paste(GRN_coef.dt_total$receptor,GRN_coef.dt_total$target,sep="_")==i
    GRN_coef.dt_total[mask_mult_rc,c("receptor_beta","receptor_pvalue")] <- t(matrix(rep(c(GRN_coef.dt_total$receptor_beta[which(mask_mult_rc)[1]],GRN_coef.dt_total$receptor_pvalue[which(mask_mult_rc)[1]]),sum(mask_mult_rc)),nrow=2))
}

# add cluster annotation
types <- c("ligand","receptor","target")
c_names <- colnames(GRN_coef.dt_total)
for (type in types){
    l_marker <- sce_mc_marker_gene_annot[as.character(GRN_coef.dt_total[,type]),2]
    l_max_e <- cluster_max_expression[GRN_coef.dt_total[,type],2]
    l_marker[sapply(l_marker,is.null)] <- NA
    l_max_e[sapply(l_max_e,is.null)] <- NA
    GRN_coef.dt_total <- cbind(GRN_coef.dt_total,
                            unlist(l_marker),
                            unlist(l_max_e))
}
colnames(GRN_coef.dt_total) <- c(c_names,paste0(rep(types,each=2),rep(c("_marker_gene","_max_expr"),length(types))))

# add expression values
c_names <- colnames(GRN_coef.dt_total)
GRN_coef.dt_total <- cbind(GRN_coef.dt_total,logcounts(sce_mc)[GRN_coef.dt_total$ligand,1:5],logcounts(sce_mc)[GRN_coef.dt_total$receptor,1:5],logcounts(sce_mc)[GRN_coef.dt_total$target,1:5])
colnames(GRN_coef.dt_total) <- c(c_names,paste(rep(c("ligand","receptor","target"),each=5),rep("cluster",15),rep(1:5,3),sep="_"))
print(sprintf("%s/GRN_coef_signalling_%s.txt",args$outdir,args$GRN_method))
write.table(GRN_coef.dt_total,file=sprintf("%s/GRN_coef_signalling_%s_%s.txt",args$outdir,args$GRN_method,args$merged),
            col.names=T,row.names=F,sep="\t",quote=F)

##################################
## Select relevant interactions ##
##################################

## select relevant interactions
## rules:
#1) p-value for ligand < args$max_pval
#2) absolute beta value for ligand > args$min_coef
#3) at least one receptor or one subunit of the receptor complex should not contradict the beta value for the ligand (ie show a beta value of the same sign or be not significant)

mask_contr_receptors <- ((GRN_coef.dt_total$ligand_beta>args$min_coef) & (GRN_coef.dt_total$ligand_pvalue<args$max_pval)) & ((GRN_coef.dt_total$receptor_beta<args$min_coef) & (GRN_coef.dt_total$receptor_pvalue<args$max_pval)) |
    ((GRN_coef.dt_total$ligand_beta< -args$min_coef) & (GRN_coef.dt_total$ligand_pvalue<args$max_pval)) & ((GRN_coef.dt_total$receptor_beta> -args$min_coef) & (GRN_coef.dt_total$receptor_pvalue<args$max_pval))

# prev
# mask_contr_receptors <- ((GRN_coef.dt_total$ligand_beta>args$min_coef) & (GRN_coef.dt_total$ligand_pvalue<args$max_pval)) & ((GRN_coef.dt_total$receptor_beta<args$min_coef) & (GRN_coef.dt_total$receptor_pvalue<args$max_pval)) |
#     ((GRN_coef.dt_total$ligand_beta<args$min_coef) & (GRN_coef.dt_total$ligand_pvalue<args$max_pval)) & ((GRN_coef.dt_total$receptor_beta>args$min_coef) & (GRN_coef.dt_total$receptor_pvalue<args$max_pval))
blacklist_ligands <- NULL
for (i in unique(paste(GRN_coef.dt_total$ligand,GRN_coef.dt_total$target,sep="_"))){
    if (sum(paste(GRN_coef.dt_total$ligand,GRN_coef.dt_total$target,sep="_")==i)==sum(mask_contr_receptors[paste(GRN_coef.dt_total$ligand,GRN_coef.dt_total$target,sep="_")==i])){
        blacklist_ligands <- c(blacklist_ligands,i)
    }
}
GRN_coef.dt_total[paste(GRN_coef.dt_total$ligand,GRN_coef.dt_total$target,sep="_")%in%blacklist_ligands,1:9]

mask_rel_interactions <- (abs(GRN_coef.dt_total$ligand_beta)>args$min_coef) & (GRN_coef.dt_total$ligand_pvalue<args$max_pval) & !(paste(GRN_coef.dt_total$ligand,GRN_coef.dt_total$target,sep="_")%in%blacklist_ligands) &
    !duplicated(paste(GRN_coef.dt_total$ligand,GRN_coef.dt_total$target,sep="_"))
print(paste0("Number of relevant interactions: ",sum(mask_rel_interactions)))

GRN_coef.dt_to_add <- GRN_coef.dt_total[mask_rel_interactions,c("ligand","target","ligand_marker_gene","ligand_max_expr")]
colnames(GRN_coef.dt_to_add) <- c("from","to","fr_marker","fr_max_expr")

##############
## Save GRN ##
##############

print(sprintf("%s/GRN_coef_signalling_%s_to_add.txt",args$outdir,args$GRN_method))
write.table(GRN_coef.dt_to_add,
            file=sprintf("%s/GRN_coef_signalling_%s_to_add_%s.txt",args$outdir,args$GRN_method,args$merged),
            col.names=T,row.names=F,sep="\t",quote=F)
