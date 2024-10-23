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
p$add_argument('--sce',       type="character",                help='SingleCellExperiment - pseudobulk')
p$add_argument('--sce_mc',       type="character",                help='SingleCellExperiment - metacells')
p$add_argument('--CellPhoneDB_data',       type="character",                help='CellPhoneDB data')
p$add_argument('--min_expr',  type="integer",            default=4,      help='Amount of cores to use')
p$add_argument('--GRN_method',  type="character", help='Method for GRN building')
p$add_argument('--GRN_type',  type="character", help='Method for GRN building')
p$add_argument('--merged',  type="character", help='Merge cluster 4 and 5')
p$add_argument('--DEG_overview',  type="character", help='Method for GRN building')
p$add_argument('--markers_TF',  type="character" ,     help='TF marker file')
p$add_argument('--ligands_to_consider',  type="character" ,     help='Ligands to consider')
p$add_argument('--targets_to_consider',  type="character" ,     help='Targets to consider')
p$add_argument('--plot_correlations',  type="logical", default=FALSE, help='Do you want to plot correlations plots')
p$add_argument('--outdir',       type="character",                help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

# # START TEST
# args <- list()
# args$sce <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/SingleCellExperiment_pseudobulk.rds"
# args$sce_mc <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/metacells/SingleCellExperiment_metacells_nodiff.rds"
# args$CellPhoneDB_data <- "/data/homes/louisc/Project_Babraham/pathways/CellPhoneDB_analysis/CellPhone_res_stat_LIGANDS_RECEPTORS.xlsx"
# args$min_expr <- 4
# args$min_coef <- 0.25
# args$max_pval <- 0.10
# args$GRN_method <- "alt5"
# args$GRN_type <- "TFandSign"
# args$merged <- "merged"
# args$DEG_overview <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/RNA/DE_res/DEG_overview.txt"
# args$markers_TF <- "/data/homes/louisc/Project_Babraham/RNA_ATAC//pseudobulk/cluster/RNA/marker_TFs/TF.markers.clusters.txt"
# args$outdir <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/NicheNet"
# # ligands_to_consider <- c("Expr","CellPhoneDB")
# args$ligands_to_consider <- c("CellPhoneDB")
# args$targets_to_consider <- c("markerTF","DEGenes","Expr")
# # END TEST

print(args)


###################
## GRN_functions ##
###################

source("/data/homes/louisc/Project_Babraham/RNA_ATAC/scripts/GRN/build_GRN_functions.R")

###############################
## Load additional resources ##
###############################

ligand_target_matrix <- readRDS("/data/homes/louisc/Project_Babraham/RNA_ATAC/NicheNet/ligand_target_matrix.rds")

#####################
## Define settings ##
#####################

dir.create(paste0(args$outdir,"/Heatmaps"))
args$outdir <- paste0(args$outdir,"/Heatmaps")

print(paste0("myfun_",args$GRN_method))
myfun_used <- get(paste0("myfun_",args$GRN_method))

#####################
## Load marker TFs ##
#####################

TF_markers <- read.table(args$markers_TF,header=T,sep="\t")$Gene

markers_TF <- NULL
for (i in 1:5){
  marker_tf_i <- read.table(sprintf("/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/RNA/marker_TFs/Marker_TF_cluster%s.txt",i),header=T,sep="\t")
  markers_TF <- rbind(markers_TF,cbind(marker_tf_i$Gene,rep(i,length(marker_tf_i$Gene))))
}
markers_TF <- data.frame(markers_TF)
colnames(markers_TF) <- c("TF","celltype")
markers_TF_undupl <- markers_TF[!duplicated(markers_TF$TF),]

###########################################
## Load RNA expression data - pseudobulk ##
###########################################

sce <- readRDS(args$sce)

##############
## DE genes ##
##############

DEG_overview <- read.table(file=args$DEG_overview, header=T, sep="\t", quote="")

c_names <- colnames(DEG_overview)
DEG_overview_bis <- DEG_overview

for (i in 1:5){
    cols_logFC <- which(grepl(i,c_names) & grepl("vs",c_names) & grepl("logFC",c_names))
    cols_FDR <-  which(grepl(i,c_names) & grepl("vs",c_names) & grepl("FDR",c_names))
    sign <- c(1,-1)[c(substr(gsub("[A-Z,a-z]+_","",c_names[cols_logFC]),1,1)==i)+1]

    DEG_overview_bis <- cbind(DEG_overview_bis,
                             (rowSums((t(t(DEG_overview_bis[,cols_logFC])*sign)>1) & (DEG_overview_bis[,cols_FDR]<0.01))>=1))
    c_names <- c(c_names,paste0("DE_clust",i))
}
colnames(DEG_overview_bis) <- c_names
for (i in 43:47){
    DEG_overview_bis[is.na(DEG_overview_bis[,i]),i] <- FALSE
}
head(DEG_overview_bis[rowSums(DEG_overview_bis[,43:46])>0,]) 

# Check total of DE genes
sum(DEG_overview_bis$DE_clust1 | DEG_overview_bis$DE_clust2 | DEG_overview_bis$DE_clust3 | DEG_overview_bis$DE_clust4 | DEG_overview_bis$DE_clust5)

# Make list 
DE_genes <- unique(c(DEG_overview$Gene[DEG_overview_bis$DE_clust1],DEG_overview$Gene[DEG_overview_bis$DE_clust2],
                     DEG_overview$Gene[DEG_overview_bis$DE_clust3],DEG_overview$Gene[DEG_overview_bis$DE_clust4],
                     DEG_overview$Gene[DEG_overview_bis$DE_clust5]))

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
dim(dat_cpdb)
head(dat_cpdb)

sum(duplicated(paste(dat_cpdb$goi,dat_cpdb$cluster,sep="_")))
dat_cpdb[duplicated(paste(dat_cpdb$goi,dat_cpdb$cluster,sep="_")),]$goi
dat_cpdb <- dat_cpdb[!duplicated(paste(dat_cpdb$goi,dat_cpdb$cluster,sep="_")),]
dim(dat_cpdb)
head(dat_cpdb)

list_cpdb <- list(list(dat_cpdb$goi[dat_cpdb$cluster==1],1),
                  list(dat_cpdb$goi[dat_cpdb$cluster==2],2),
                  list(dat_cpdb$goi[dat_cpdb$cluster==3],3),
                  list(dat_cpdb$goi[dat_cpdb$cluster==4],4),
                  list(dat_cpdb$goi[dat_cpdb$cluster==5],5),
                  list(unique(dat_cpdb$goi),"X"))
names(list_cpdb) <- c("CellPhoneDB_clusters1","CellPhoneDB_clusters2","CellPhoneDB_clusters3",
                      "CellPhoneDB_clusters4","CellPhoneDB_clusters5","CellPhoneDB_clustersX")

genesets <- list_cpdb

#####################
## Expressed genes ##
#####################

expressed_genes <- rownames(sce)[rowMeans(logcounts(sce)[sce$celltype%in%c(1:5),])>=args$min_expr]
for (cluster in 1:5){
    sce_cluster <- sce[,sce$celltype %in% cluster]
    expressed_genes <- unique(c(expressed_genes,rownames(sce_cluster)[logcounts(sce_cluster)>=args$min_expr]))
}

##############
## Heatmaps ##
##############

for (ltc in args$ligands_to_consider){
    for (ttc in args$targets_to_consider){
        print(ltc)
        print(ttc)

        if (ltc=="Expr"){
            ligands <- expressed_genes
        } else if (ltc=="CellPhoneDB") {
            ligands <- genesets[["CellPhoneDB_clustersX"]][[1]]
        } else {
            stop("Not yet implemented...")
        } 

        if (ttc=="markerTF"){
            targets <- unique(markers_TF$TF)
        } else if (ttc=="Expr"){
            targets <- expressed_genes
        } else if (ttc=="DEGenes"){
           targets <- DE_genes
        } else {
            stop("Not yet implemented...")
        }

        ltm_sub <- ligand_target_matrix[rownames(ligand_target_matrix)%in%targets,
                                        colnames(ligand_target_matrix)%in%ligands] 

        print(dim(ltm_sub))
        print(summary(as.numeric(ltm_sub)))

        jpeg(paste0(args$outdir,"/hist_",ltc,"_ligand_to_",ttc,"_targets.jpg"))
        hist(ltm_sub,breaks=50)
        dev.off()                                  

        dist_ligands = dist(ltm_sub %>% t(), method = "euclidian")
        hclust_ligands = hclust(dist_ligands, method = "ward.D2")
        if (ttc == "markerTF"){
            ltm_sub <- ltm_sub[TF_markers,hclust_ligands$order]
        } else {
            dist_targets <- dist(ltm_sub, method = "euclidian")
            h_clust_targets <- hclust(dist_targets, method = "ward.D2")
            ltm_sub <- ltm_sub[h_clust_targets$order,hclust_ligands$order]
        }

        write.table(ltm_sub %>% t(), file=paste0(args$outdir,"/",ltc,"_ligand_to_",ttc,"_targets_",args$GRN_method,"_",args$merged,".txt"),
                    row.names=TRUE,col.names=TRUE,sep="\t")

        p_ltm_sub <- ltm_sub %>% t() %>% make_heatmap_ggplot("Ligands","Targets", color = "mediumvioletred", x_axis_position = "top",
                                                                            legend_title = "Prior interaction potential") 

        rel_size_row <- min(max(c(1,floor(ncol(ltm_sub)/nrow(ltm_sub)*0.8))),5)
        rel_size_col <- min(max(c(1,floor(nrow(ltm_sub)/ncol(ltm_sub)*0.8))),5)

        jpeg(paste0(args$outdir,"/heatmap_",ltc,"_ligand_to_",ttc,"_targets.jpg"),
             width=1024*rel_size_col,height=1024*rel_size_row)
        print(p_ltm_sub)
        dev.off()       

        cutoffs <- c(0.05,0.075,0.10)
        for (cutoff in cutoffs){
            print(paste0("Cutoff: ",cutoff))
            ltm_sub_cutoff <- ltm_sub[,colSums(ltm_sub>cutoff)>0]
            ltm_sub_cutoff <- ltm_sub_cutoff[rowSums(ltm_sub_cutoff>cutoff)>0,]
            print(dim(ltm_sub_cutoff))

            if(nrow(ltm_sub_cutoff)<length(TF_markers)){
                dat_NA <- matrix(rep(0,(length(TF_markers)-nrow(ltm_sub_cutoff))*ncol(ltm_sub_cutoff)),ncol=ncol(ltm_sub_cutoff))
                rownames(dat_NA) <- TF_markers[!(TF_markers%in%rownames(ltm_sub_cutoff))]
                colnames(dat_NA) <- colnames(ltm_sub_cutoff)
                ltm_sub_cutoff <- rbind(ltm_sub_cutoff,dat_NA)
            }

            dist_ligands = dist(ltm_sub_cutoff %>% t(), method = "euclidian")
            hclust_ligands = hclust(dist_ligands, method = "ward.D2")
            if (ttc == "markerTF"){
                ltm_sub_cutoff <- ltm_sub_cutoff[TF_markers,hclust_ligands$order]
            } else {
                dist_targets = dist(ltm_sub_cutoff, method = "euclidian")
                hclust_targets = hclust(dist_targets, method = "ward.D2")
                ltm_sub_cutoff <- ltm_sub_cutoff[hclust_targets$order,hclust_ligands$order]
            }

            rel_size_row <- min(max(c(1,floor(ncol(ltm_sub_cutoff)/nrow(ltm_sub_cutoff)*0.8))),5)
            rel_size_col <- min(max(c(1,floor(nrow(ltm_sub_cutoff)/ncol(ltm_sub_cutoff)*0.8))),5)

            p_ltm_sub_cutoff <- ltm_sub_cutoff %>% t() %>% make_heatmap_ggplot("Ligands","Targets", color = "mediumvioletred", x_axis_position = "top",
                                                                                legend_title = "Prior interaction potential")

            jpeg(paste0(args$outdir,"/heatmap_",ltc,"_ligand_to_",ttc,"_targets_cutoff_",cutoff,".jpg"),
                 width=1024*rel_size_col,height=1024*rel_size_row)
            print(p_ltm_sub_cutoff)
            dev.off()         
        }        
    }
} 