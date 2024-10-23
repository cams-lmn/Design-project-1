##########################################
##                                      ##
##  Build_GRN_metacells_trajectories.R  ##
##                                      ##
##########################################

source("/data/homes/louisc/Project_Babraham/RNA_ATAC/scripts/Settings.R")

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--trajectory_name',  type="character",              help='Trajectory name')
p$add_argument('--sce',       type="character",                help='SingleCellExperiment')
p$add_argument('--metadata',       type="character",                help='Metadata metacells')
p$add_argument('--sign_connections',       type="character",                help='Metadata metacells')
p$add_argument('--tf2gene_virtual_chip',       type="character",                help='Links between tfs and genes based on virtual chipseq')
p$add_argument('--max_distance',  type="integer",            default=1e5,      help='Maximum distance for a linkage between a peak and a gene')
p$add_argument('--min_chip_score',  type="double",            default=0.25,      help='Minimum in silico ChIP-seq score')
p$add_argument('--min_coef',  type="double",            default=0.24,      help='Minimal coefficient for linear regression')
p$add_argument('--max_pval',  type="double",            default=0.10,      help='Maximal pvalue for linear regression')
p$add_argument('--ncores',  type="integer",            default=4,      help='Amount of cores to use')
p$add_argument('--plot_correlations',  type="logical", default=TRUE, help='Do you want to plot correlations plots')
p$add_argument('--GRN_method',  type="character", help='Method for GRN building')
p$add_argument('--GRN_type',  type="character", help='Method for GRN building')
p$add_argument('--merged',  type="character", help='Merge cluster 4 and 5')
p$add_argument('--markers_TF',  type="character", help='TF marker file')
p$add_argument('--DEG_overview',  type="character", help='Method for GRN building')
p$add_argument('--outdir', type="character", help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

# I/O
# args <- list()
# args$trajectory_name <- "N2P"
# args$sce <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/metacells/SingleCellExperiment_metacells_nodiff.rds"
# args$metadata <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/clustering/PeakMatrix/sample_metadata_nodiff_after_clustering.txt.gz"
# args$sign_connections <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/NicheNet/GRN_coef_signalling_alt7_to_add_merged.txt"
# args$tf2gene_virtual_chip <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/gene_regulatory_network/metacells/CISBP/TF2gene_after_virtual_chip.txt.gz"
# args$max_distance <- 5e4
# args$min_chip_score <- 0.06
# args$ncores <- 10
# args$plot_correlations <- TRUE
# args$GRN_method <- "conv"
# args$GRN_method <- "TFandSign"
# args$DEG_overview <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/RNA/DE_res/DEG_overview.txt"
# args$merged <- "merged"
# args$markers_TF <-  "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/RNA/marker_TFs/TF.markers.clusters.txt"
# args$outdir <-  sprintf("/data/homes/louisc/Project_Babraham/RNA_ATAC/gene_regulatory_network/metacells/trajectories/N2P/%s/%s/%s",args$GRN_type,args$GRN_method,args$merged)

print(args)

###################
## GRN_functions ##
###################

source("/data/homes/louisc/Project_Babraham/RNA_ATAC/scripts/GRN/build_GRN_functions.R")

#####################
## Define settings ##
#####################

print(paste0("Method for GRN coef calculation: ",args$GRN_method))

if(grepl("full",args$GRN_method)){
  args$plot_correlations <- FALSE
}

dir.create(args$outdir, showWarnings = F, recursive = T)

print(paste0("myfun_",args$GRN_method))
myfun_used <- get(paste0("myfun_",args$GRN_method))

if (args$plot_correlations){
  dir.create(sprintf("%s/plots_cor_score%s",args$outdir,args$min_chip_score))
}

opts$celltype_trajectory_dic <- list(
  "N2P" = c(1,2,3,4,5)
)

stopifnot(args$trajectory_name%in%names(opts$celltype_trajectory_dic))
celltypes.to.plot <- opts$celltype_trajectory_dic[[args$trajectory_name]]

##############################
## Load RNA expression data ##
##############################

sce_mc <- readRDS(args$sce)

########################
## Cluster annotation ##
########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE]

sample_metadata <- data.frame(sample_metadata)
rownames(sample_metadata) <- sample_metadata$sample
sample_metadata_metacells <- sample_metadata[colnames(sce_mc),]

colData(sce_mc)$cluster <- sample_metadata_metacells$cluster

#####################
## Load marker TFs ##
#####################

# TF markers with clusters
markers_TF <- NULL
for (i in celltypes.to.plot){
  marker_tf_i <- read.table(sprintf("/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/RNA/marker_TFs/Marker_TF_cluster%s.txt",i),header=T,sep="\t")
  markers_TF <- rbind(markers_TF,cbind(marker_tf_i$Gene,rep(i,length(marker_tf_i$Gene))))
}
markers_TF <- data.frame(markers_TF)
colnames(markers_TF) <- c("TF","celltype")
markers_TF_undupl <- markers_TF[!duplicated(markers_TF$TF),]
print("Marker TFs")
print(nrow(markers_TF_undupl))
print(head(markers_TF_undupl))

####################################################
## Load TF2gene links based on in silico ChIP-seq ##
####################################################

tf2gene_chip.dt <- fread(args$tf2gene_virtual_chip) %>%
  .[chip_score>=args$min_chip_score & dist<=args$max_distance] %>% 
  .[,c("tf","gene")] %>% unique # Only keep TF-gene links

print(sprintf("Number of TFs: %s",length(unique(tf2gene_chip.dt$tf))))
print(sprintf("Number of genes: %s",length(unique(tf2gene_chip.dt$gene))))
print(head(tf2gene_chip.dt))
tf2gene_chip.dt_ori <- tf2gene_chip.dt

###########################
## Load signalling links ##
###########################

print("Loading signalling interactions")

sign2TF <- fread(args$sign_connections)
print(dim(sign2TF))
print(head(sign2TF))

sign_senders <- unique(sign2TF[sign2TF$to%in%markers_TF_undupl$TF,]$from)

print(sign_senders)

#############
## Combine ##
#############

tf2gene_chip.dt <- rbind(tf2gene_chip.dt,sign2TF[,1:2],use.names=FALSE)

##########################
## Filter TFs and genes ##
##########################

TFs <- intersect(unique(tf2gene_chip.dt$tf),toupper(rownames(sce_mc)))
print(head(TFs))
genes <- intersect(unique(tf2gene_chip.dt$gene),rownames(sce_mc))

if (grepl("full",args$outdir)){
  tf2gene_chip.dt <- tf2gene_chip.dt[tf%in%TFs & gene%in%genes,]
} else {
  tf2gene_chip.dt <- tf2gene_chip.dt[tf%in%TFs & gene%in%TFs,]
}

# Fetch RNA expression matrices
rna_tf.mtx <- logcounts(sce_mc)[unique(c(tf2gene_chip.dt$tf,sign2TF$from)),]
rownames(rna_tf.mtx) <- toupper(rownames(rna_tf.mtx))
rna_targets.mtx <- logcounts(sce_mc)[unique(c(tf2gene_chip.dt$gene,sign2TF$to)),]
print(dim(rna_tf.mtx))
print(head(rna_tf.mtx[,1:6]))

# Final set
print("After selection of TF and genes")
print(dim(tf2gene_chip.dt))
print(sprintf("Number of TFs: %s",length(TFs)))
print(sprintf("Number of genes: %s",length(genes)))

############################
## Get cluster annotation ##
############################

DEG_overview <- read.table(file=args$DEG_overview, header=T, sep="\t", quote="")

# annotate markers 
marker_gene_annot <- cbind(DEG_overview[,1],lapply(apply(DEG_overview[,grepl("clust[1-5]_spec",colnames(DEG_overview))],1,which),paste,collapse=","))
marker_gene_annot <- data.frame(marker_gene_annot)
colnames(marker_gene_annot) <- c("gene","cluster")
rownames(marker_gene_annot) <- marker_gene_annot$gene
marker_gene_annot <- marker_gene_annot[marker_gene_annot$cluster!="",]
marker_gene_annot$cluster <- gsub(",[0-9]$","",marker_gene_annot$cluster)

# max gene expressions
cluster_max_expression <- cbind(DEG_overview[,1],apply(DEG_overview[,2:6],1,which.max))
cluster_max_expression <- data.frame(cluster_max_expression)
colnames(cluster_max_expression) <- c("gene","cluster")
rownames(cluster_max_expression) <- cluster_max_expression$gene

# marker genes
Gene_markers <- DEG_overview$Gene[rowSums((DEG_overview[,grepl("clust[0-5]_spec",colnames(DEG_overview))]) & (!is.na(DEG_overview[,grepl("clust[0-5]_spec",colnames(DEG_overview))])))>0]

gene_cluster_annot <- data.frame(gene=unique(c(TFs,genes)),
                                  cluster=NA)                                   
gene_cluster_annot$cluster[gene_cluster_annot$gene%in%marker_gene_annot$gene] <- marker_gene_annot[gene_cluster_annot$gene[gene_cluster_annot$gene%in%marker_gene_annot$gene],]$cluster
gene_cluster_annot$cluster[is.na(gene_cluster_annot$cluster)] <- cluster_max_expression[gene_cluster_annot$gene[is.na(gene_cluster_annot$cluster)],]$cluster
rownames(gene_cluster_annot) <- gene_cluster_annot$gene

############################
## Merge cluster 5 to 4  ##
############################

if(args$merged=="merged"){
  gene_cluster_annot$cluster[gene_cluster_annot$cluster==5] <- 4
}

####################
## run regression ##
####################

tmp <- tf2gene_chip.dt
print(dim(tmp))
list_results <- mclapply(genes,myfun_used,mc.cores=args$ncores,mc.preschedule=T)
print(head(list_results))
GRN_coef.dt <- do.call(rbind,list_results)
colnames(GRN_coef.dt) <- c("sender","target","beta","pvalue","clusters_for_calc")

###############################
## Add origin of interaction ##
###############################

GRN_coef.dt$origin <- "in-silico_ChipSeq"
GRN_coef.dt$origin[GRN_coef.dt$sender%in%sign2TF$from] <- "signalling"

print(table(GRN_coef.dt$origin))

#####################
## Add sender info ##
#####################

GRN_coef.dt$sender_class <- "other"
GRN_coef.dt$sender_class[GRN_coef.dt$sender%in%sign2TF$from] <- "signalling"
GRN_coef.dt$sender_class[(GRN_coef.dt$sender%in%TFs) & !(GRN_coef.dt$sender%in%sign2TF$from)] <- "TF"

GRN_coef.dt$sender_marker <- "not-marker"
GRN_coef.dt$sender_marker[(GRN_coef.dt$sender%in%Gene_markers)] <- "marker" 

GRN_coef.dt$sender_cluster_marker <- NA
GRN_coef.dt$sender_cluster_marker[GRN_coef.dt$sender%in%marker_gene_annot$gene] <- marker_gene_annot[GRN_coef.dt$sender[GRN_coef.dt$sender%in%marker_gene_annot$gene],]$cluster

GRN_coef.dt$sender_cluster_max_expr <- NA
GRN_coef.dt$sender_cluster_max_expr[GRN_coef.dt$sender%in%cluster_max_expression$gene] <- cluster_max_expression[GRN_coef.dt$sender[GRN_coef.dt$sender%in%cluster_max_expression$gene],]$cluster

#####################
## Add target info ##
#####################

GRN_coef.dt$target_class <- "other"
GRN_coef.dt$target_class[GRN_coef.dt$target%in%sign2TF$from] <- "signalling"
GRN_coef.dt$target_class[(GRN_coef.dt$target%in%TFs) & !(GRN_coef.dt$target%in%sign2TF$from)] <- "TF"

GRN_coef.dt$target_marker <- "not-marker"
GRN_coef.dt$target_marker[(GRN_coef.dt$target%in%Gene_markers)] <- "marker" 

GRN_coef.dt$target_cluster_marker <- NA
GRN_coef.dt$target_cluster_marker[GRN_coef.dt$target%in%marker_gene_annot$gene] <- marker_gene_annot[GRN_coef.dt$target[GRN_coef.dt$target%in%marker_gene_annot$gene],]$cluster

GRN_coef.dt$target_cluster_max_expr <- NA
GRN_coef.dt$target_cluster_max_expr[GRN_coef.dt$target%in%cluster_max_expression$gene] <- cluster_max_expression[GRN_coef.dt$target[GRN_coef.dt$target%in%cluster_max_expression$gene],]$cluster

#######################
## Consensus cluster ##
#######################

GRN_coef.dt$sender_cluster <- NA
GRN_coef.dt$sender_cluster[(is.na(GRN_coef.dt$sender_cluster)) & !(is.na(GRN_coef.dt$sender_cluster_marker))] <- GRN_coef.dt$sender_cluster_marker[(is.na(GRN_coef.dt$sender_cluster)) & !(is.na(GRN_coef.dt$sender_cluster_marker))]
GRN_coef.dt$sender_cluster[(is.na(GRN_coef.dt$sender_cluster)) & !(is.na(GRN_coef.dt$sender_cluster_max_expr))] <- GRN_coef.dt$sender_cluster_max_expr[(is.na(GRN_coef.dt$sender_cluster)) & !(is.na(GRN_coef.dt$sender_cluster_max_expr))]

GRN_coef.dt$target_cluster <- NA
GRN_coef.dt$target_cluster[(is.na(GRN_coef.dt$target_cluster)) & !(is.na(GRN_coef.dt$target_cluster_marker))] <- GRN_coef.dt$target_cluster_marker[(is.na(GRN_coef.dt$target_cluster)) & !(is.na(GRN_coef.dt$target_cluster_marker))]
GRN_coef.dt$target_cluster[(is.na(GRN_coef.dt$target_cluster)) & !(is.na(GRN_coef.dt$target_cluster_max_expr))] <- GRN_coef.dt$target_cluster_max_expr[(is.na(GRN_coef.dt$target_cluster)) & !(is.na(GRN_coef.dt$target_cluster_max_expr))]

##############
## Save GRN ##
##############

print("Final GRN coef")
print(dim(GRN_coef.dt))
print(head(GRN_coef.dt))

for (i in 5:ncol(GRN_coef.dt)){
  print(colnames(GRN_coef.dt)[i])
  print(table(GRN_coef.dt[,i]))
}

GRN_coef.dt <- data.table(GRN_coef.dt)

fwrite(GRN_coef.dt, file.path(args$outdir,sprintf('global_chip_GRN_sign_%s_coef_score%s_%s.txt.gz',args$GRN_method,args$min_chip_score,args$merged)), sep="\t")

########################
## Create cleaned GRN ##
########################

print("GRN cleaned version")

GRN_coef.dt$beta <- as.numeric(GRN_coef.dt$beta)
GRN_coef.dt$pvalue <- as.numeric(GRN_coef.dt$pvalue)

## filter on coef & pvalue
################################

print("Filter coef & pvalue")
GRN_coef.dt_clean <- GRN_coef.dt %>%
  .[pvalue<args$max_pval & abs(beta)>=args$min_coef]
print(dim(GRN_coef.dt_clean))
print(head(GRN_coef.dt_clean))

## Remove self activatory int 
################################

print("Remove self activatory int")
GRN_coef.dt_clean <- GRN_coef.dt_clean[sender!=target]
print(dim(GRN_coef.dt_clean))
print(head(GRN_coef.dt_clean))

#################################
## Determine key nodes for GRN ##
#################################

# TF connections that have no incoming activatory signals (except cluster 1) 
print("Nodes that have no incoming activatory signals")
GRN_coef.dt_act <- GRN_coef.dt_clean[sender%in%c(markers_TF_undupl$TF,sign_senders) & beta>=args$min_coef]
# print(head(GRN_coef.dt_act))
# print(sort(unique(GRN_coef.dt_act$target)))

genes_no_act <- markers_TF_undupl$TF[!((markers_TF_undupl$TF%in%unique(GRN_coef.dt_act$target)) |
                  (markers_TF_undupl$celltype==1))]

print(sort(genes_no_act))

# Signalling connections that have no marker TF targets
print("Signalling nodes that have no marker TF as downstream target")
ligand_no_marker_target <- sign2TF$from
ligand_no_marker_target <- ligand_no_marker_target[!(ligand_no_marker_target%in%sign2TF[sign2TF$to%in%markers_TF_undupl$TF,]$from)]
print(sort(ligand_no_marker_target))
print("Signalling nodes that have  no incoming activatory signals")
ligand_no_act <- sign2TF[!duplicated(sign2TF$from),]
ligand_no_act <- ligand_no_act$from[!((ligand_no_act$from%in%unique(GRN_coef.dt_act$target)) |
                  (ligand_no_act$fr_max_expr==1))]
print(sort(ligand_no_act))

print("Determine key components GRN")

key_components <- unique(c(markers_TF_undupl$TF,sign2TF$from))
print(length(key_components))
print(key_components)
print("filtering...")
key_components <- key_components[!(key_components%in%c(genes_no_act,ligand_no_act,ligand_no_marker_target))]    
print(length(key_components))
print(key_components)

## Select key components
################################

print("Select key components")
GRN_coef.dt_clean <- GRN_coef.dt_clean[(sender%in%key_components) & (target%in%key_components)]
print(dim(GRN_coef.dt_clean))
print(head(GRN_coef.dt_clean))

####################
## Save clean GRN ##
####################

fwrite(GRN_coef.dt_clean, file.path(args$outdir,sprintf('global_chip_GRN_sign_%s_coef_score%s_%s_cleaned.txt.gz',args$GRN_method,args$min_chip_score,args$merged)), sep="\t")

##################
# Print warnings #
##################

warnings()