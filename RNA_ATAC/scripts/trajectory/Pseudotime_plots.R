########################
##                    ##
##  Pseudotime plots  ##
##                    ##
########################

source("/data/homes/louisc/Project_Babraham/RNA_ATAC/scripts/Settings.R")
library("viridis")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--trajectory_name',  type="character",              help='Trajectory name')
p$add_argument('--pseudotime',    type="character",    help='Pseudotime file')
p$add_argument('--sce_metacell',    type="character",    help='Expression data on metacell level')
p$add_argument('--atac_metacell_TSS',    type="character",    help='Accessibility data on metacell level (TSS)')
p$add_argument('--atac_metacell_Distal',    type="character",    help='Accessibility data on metacell level (Distal)')
p$add_argument('--goi',          type="character",  nargs="+",    help='Genes of interest')
p$add_argument('--umap_coord',    type="character",    help='UMAP coordinates for visualisation')
p$add_argument('--grn_coef',       type="character",                help='GRN coefficients')
p$add_argument('--DE_gene_clusters',    type="character",    help='Expression data on metacell level')
p$add_argument('--DA_gene_clusters',    type="character",    help='Expression data on metacell level')
p$add_argument('--min_coef',  type="double",            default=0.25,      help='Minimal coefficient for linear regression')
p$add_argument('--max_pval',  type="double",            default=0.10,      help='Maximal pvalue for linear regression')
p$add_argument('--DEG_overview',    type="character",    help='Differential expression analysis overview')
p$add_argument('--DAG_overview',    type="character",    help='Differential accessibility analysis overview')
p$add_argument('--incl_samples',        type="character",                               help='Suffix')
# p$add_argument('--TF',    type="character",    help='Transcription factor file')
# p$add_argument('--seed',            type="integer",     default=42,             help='Random seed')
p$add_argument('--gap_rows',            type="integer",   nargs="+",    default=NULL,     help='Random seed')
p$add_argument('--goi_ordered',  type="logical", default=FALSE, help='Is the goi set ordered?')
# p$add_argument('--matrix',          type="character",  nargs="+",    help='Matrix to use')
# p$add_argument('--mofa_metadata',    type="character",    help='MOFA metadata')
# p$add_argument('--factors',    type="character",    help='MOFA factors')
p$add_argument('--outdir',   type="character",    help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

# # START TEST
# args <- list()
# args$trajectory_name <- "N2P"
# args$pseudotime <- sprintf("/data/homes/louisc/Project_Babraham/RNA_ATAC/trajectories/%s/%s_trajectory_nodiff.txt.gz",args$trajectory_name,args$trajectory_name)
# args$sce_metacell <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/metacells/SingleCellExperiment_metacells_nodiff.rds"
# args$atac_metacell_TSS <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/metacells/GeneScoreMatrix_TSS_summarized_experiment_metacells_nodiff.rds"
# args$atac_metacell_Distal <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/metacells/GeneScoreMatrix_distal_summarized_experiment_metacells_nodiff.rds"
# # args$goi <- c("KLF4", "KLF5", "SPIC", "TFCP2L1", "DNMT3L", "FGF4", "KLF17", "SOX15", "DPPA3", "DPPA5", "SALL4", "TFAP2C", "FBP1", "ARGFX",
# #               "DLL3", "FOXH1", "NODAL", "GDF3", "SPP1", "PRDM14", "DPPA2", "DPPA4", "ETV4", "ETV5", "SALL2", "SOX11", "ZIC2", "SFRP2",
# #               "MYC", "FZD7", "FGF2", "HES1", "FZD2", "OTX2", "FST", "TCF15", "CDH2", "TCF7L1")
# args$goi <- c("SPIC" , "DPPA5" , "DPPA3" , "KLF4" , "DNMT3L" , "ARGFX" , "TFCP2L1" , "PRDM14" , "FBP1" , "KLF17" , "KLF5" , 
#               "GDF3" , "FGF4" , "NODAL" , "DLL3" , "GBX2" , "ETV4" , "PODXL" , "SFRP2" , "ETV5" , "FGF2" , "HES1" , "SALL4" , 
#               "ZIC2" , "SOX11", "SALL2", "FST" , "MYC" , "SPP1" , "FZD7" , "TCF15" , "CDH2" , "TCF7L1" , "SALL1")
# args$gap_rows <- which(args$goi%in%c("NODAL","GBX2"))
# args$gap_rows <- NULL
# args$goi_ordered <- TRUE
# args$outdir <- sprintf("/data/homes/louisc/Project_Babraham/RNA_ATAC/trajectories/%s",args$trajectory_name)
# args$umap_coord <- "/data/homes/louisc/Project_Babraham/RNA_ATAC//mofa/pdf/umap_nodiff.txt.gz"
# args$DE_gene_clusters <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/RNA/DE_res/genes_gene_clustering_ordered.txt"
# args$DA_gene_clusters <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/GeneScoreMatrix_distal/DE_res/genes_gene_clustering_ordered.txt"
# args$grn_coef <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/gene_regulatory_network/metacells/trajectories/N2P/full/global_chip_GRN_coef_score0.06.txt.gz"
# args$min_coef <- 0.25
# args$max_pval <- 0.10
# args$DEG_overview <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/RNA/DE_res/DEG_overview.txt"
# args$DAG_overview <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/GeneScoreMatrix_distal/DE_res/DEG_overview.txt"
# # END TEST

print(args)

####################
## Initialisation ##
####################

# Set default zip commando
Sys.setenv(R_ZIPCMD="zip")

# cols_scale <- (brewer.pal(9, "YlOrRd"))
cols_scale <- inferno(8)

ann_col_colors <- list(
  dos=opts$color_scheme,
  cluster=opts$celltype.colors
)

###############
## Functions ##
###############

count.gaps <- function(x,clusters){
  clusters %>%
    filter(cluster<=x) %>%
    nrow -> y
  return(y)
}

#####################
## Load Pseudotime ##
#####################

print("Pseudotime")

pseudotime.dt <- fread(args$pseudotime)

pseudotime.dt$PT <- pseudotime.dt$PC1+pseudotime.dt$DC1*100
pseudotime.dt <- pseudotime.dt %>% setorder(PT)

print(dim(pseudotime.dt))
print(str(pseudotime.dt))
print(head(pseudotime.dt))

ann_col_colors$dos <- ann_col_colors$dos[1:length(unique(pseudotime.dt$dos))]
names(ann_col_colors$dos) <- levels(as.factor(pseudotime.dt$dos))[sort(as.numeric(gsub("^d","",levels(as.factor(pseudotime.dt$dos)))),index.return=T)$ix]
ann_col_colors$cluster <- ann_col_colors$cluster[1:length(unique(pseudotime.dt$cluster))]
names(ann_col_colors$cluster) <- levels(as.factor(pseudotime.dt$cluster))

print(ann_col_colors)

##########################
## Load RNA & ATAC data ##
##########################

print("RNA & ATAC data")

sce <- readRDS(args$sce_metacell)
atac.sce.TSS <- readRDS(args$atac_metacell_TSS)
atac.sce.Distal <- readRDS(args$atac_metacell_Distal)

print(sce)
print(atac.sce.TSS)
print(atac.sce.Distal)

assay(atac.sce.TSS,"logcounts") <- log2(1e6*(sweep(assay(atac.sce.TSS),2,colSums(assay(atac.sce.TSS)),"/"))+1)
assay(atac.sce.Distal,"logcounts") <- log2(1e6*(sweep(assay(atac.sce.Distal),2,colSums(assay(atac.sce.Distal)),"/"))+1)

if (sum(colnames(assay(sce,"logcounts"))==colnames(assay(atac.sce.TSS,"logcounts")))!=ncol(sce)){
  stop("RNAseq & ATACseq not metacell info not in same order")
}
if (sum(colnames(assay(sce,"logcounts"))==colnames(assay(atac.sce.Distal,"logcounts")))!=ncol(sce)){
  stop("RNAseq & ATACseq not metacell info not in same order")
}

#############################
## Check availability goi  ##
#############################

if (sum(args$goi%in%rownames(sce))!=length(args$goi)){
  print("Some genes of interest are not present in expression data")
  print(args$goi[!(args$goi%in%rownames(sce))])

  stop("Check goi!")
}

#################################################
## Construct data object for genes of interest ##
#################################################

print("Data object goi")

val_per_TF <- 3

dat_all <- NULL
dat_all_z <- NULL
dat_all_minmax <- NULL
for (i in 1:length(args$goi)){
  data_goi <- rbind(assay(sce,"logcounts")[args$goi[i],],
                    assay(atac.sce.TSS,"logcounts")[args$goi[i],],
                    assay(atac.sce.Distal,"logcounts")[args$goi[i],]) %>% data.frame
  colnames(data_goi) <- colnames(sce)

  data_goi_sorted <- data_goi[,pseudotime.dt$cell]
  dat_all <- rbind(dat_all,data_goi_sorted)

  # z-scores
  dat_all_z <- rbind(dat_all_z,t(apply(data_goi_sorted, 1, z_score)))

  # minmax normalisation
  dat_all_minmax <- rbind(dat_all_minmax,t(apply(data_goi_sorted, 1, minmax.10pct.normalisation)))
}
rownames(dat_all) <- paste0(rep(args$goi,each=val_per_TF),rep("_",length(args$goi)*val_per_TF),rep(c("exprs","atac_TSS","atac_distal"),length(args$goi)))
rownames(dat_all_z) <- paste0(rep(args$goi,each=val_per_TF),rep("_",length(args$goi)*val_per_TF),rep(c("exprs","atac_TSS","atac_distal"),length(args$goi)))
rownames(dat_all_minmax) <- paste0(rep(args$goi,each=val_per_TF),rep("_",length(args$goi)*val_per_TF),rep(c("exprs","atac_TSS","atac_distal"),length(args$goi)))

print(dim(dat_all))
print(dat_all[1:6,1:6])

#############
## Heatmap ##
#############

dat_all_z %>%
  pheatmap(cluster_rows = F,
           cluster_cols = F,
           show_rownames = T,
           show_colnames = F,
           scale = "none",
           fontsize_row = 8,
           gaps_row = seq(from=3,to=nrow(dat_all_z),by=3),
           color = colorRampPalette(cols_scale)(100),
           filename = sprintf("%s/metacell_plot_per_gene_z.pdf",args$outdir), width = 10, height = 10)

dat_all_minmax %>%
  pheatmap(cluster_rows = F,
           cluster_cols = F,
           show_rownames = T,
           show_colnames = F,
           scale = "none",
           fontsize_row = 8,
           gaps_row = seq(from=3,to=nrow(dat_all_minmax),by=3),
           color = colorRampPalette(cols_scale)(100),
           filename = sprintf("%s/metacell_plot_per_gene_minmax.pdf",args$outdir), width = 10, height = 10)

# ############################
# ## Heatmap - per modality ##
# ############################

print("Heatmap - per modality")

# expression

dat_all_minmax_exprs <- dat_all_minmax[grepl("exprs",rownames(dat_all_minmax)),]
rownames(dat_all_minmax_exprs) <- gsub("_exprs","",rownames(dat_all_minmax_exprs))

if (args$goi_ordered){
  order_goi <- args$goi
} else {
  dist_dat_all_minmax_exprs = dist(dat_all_minmax_exprs, method = "euclidian")
  hclust_dat_all_minmax_exprs  = hclust(dist_dat_all_minmax_exprs, method = "ward.D2")
  order_goi <- hclust_dat_all_minmax_exprs$order
}

dat_all_minmax_exprs <- dat_all_minmax_exprs[order_goi,]

ann_col_data <- data.frame(pseudotime.dt[,c("dos","cluster")])
for (i in 1:ncol(ann_col_data)){
  ann_col_data[,i] <- as.factor(ann_col_data[,i])
}
rownames(ann_col_data) <- pseudotime.dt$cell

print(dim(dat_all_minmax_exprs))
print(dim(ann_col_data))
print(str(ann_col_data))
print(head(ann_col_data))

print(sum(rownames(ann_col_data)==colnames(dat_all_minmax_exprs)))

dat_all_minmax_exprs  %>%
  pheatmap(cluster_rows = F,
           cluster_cols = F,
           show_rownames = T,
           show_colnames = F,
           annotation_col=ann_col_data,
           annotation_colors=ann_col_colors,
           scale = "none",
           gaps_row = args$gap_rows,
           fontsize_row = 8,
           color = colorRampPalette(cols_scale)(100),
           filename = sprintf("%s/metacell_exprs_plot_per_gene_minmax.pdf",args$outdir), width = 10, height = 10)

dat_all_minmax_exprs  %>%
  pheatmap(cluster_rows = F,
           cluster_cols = F,
           show_rownames = T,
           show_colnames = F,
           annotation_col=ann_col_data,
           annotation_colors=ann_col_colors,
           scale = "none",
           gaps_row = args$gap_rows,
           fontsize_row = 8,
           color = colorRampPalette(cols_scale)(100),
           filename = sprintf("%s/metacell_exprs_plot_per_gene_minmax.png",args$outdir), width = 10, height = 10)

dat_all_minmax_exprs  %>%
  pheatmap(cluster_rows = F,
           cluster_cols = F,
           show_rownames = F,
           show_colnames = F,
           annotation_col=ann_col_data,
           annotation_colors=ann_col_colors,
           legend= F,
           scale = "none",
           gaps_row = args$gap_rows,
           fontsize_row = 8,
           color = colorRampPalette(cols_scale)(100),
           filename = sprintf("%s/metacell_exprs_plot_per_gene_minmax_nolabs.pdf",args$outdir), width = 10, height = 10)

dat_all_minmax_exprs  %>%
  pheatmap(cluster_rows = F,
           cluster_cols = F,
           show_rownames = F,
           show_colnames = F,
           annotation_col=ann_col_data,
           annotation_colors=ann_col_colors,
           legend= F,
           scale = "none",
           gaps_row = args$gap_rows,
           fontsize_row = 8,
           color = colorRampPalette(cols_scale)(100),
           filename = sprintf("%s/metacell_exprs_plot_per_gene_minmax_nolabs.png",args$outdir), width = 10, height = 10)

# accessibility - distal

dat_all_minmax_access_distal <- dat_all_minmax[grepl("atac_distal",rownames(dat_all_minmax)),]
rownames(dat_all_minmax_access_distal) <- gsub("_atac_distal","",rownames(dat_all_minmax_access_distal))

dat_all_minmax_access_distal <- dat_all_minmax_access_distal[order_goi,]

dat_all_minmax_access_distal  %>%
  pheatmap(cluster_rows = F,
           cluster_cols = F,
           show_rownames = T,
           show_colnames = F,
           annotation_col=ann_col_data,
           annotation_colors=ann_col_colors,
           scale = "none",
           gaps_row = args$gap_rows,
           fontsize_row = 8,
           color = colorRampPalette(cols_scale)(100),
           filename = sprintf("%s/metacell_access_distal_plot_per_gene_minmax.pdf",args$outdir), width = 10, height = 10)

dat_all_minmax_access_distal  %>%
  pheatmap(cluster_rows = F,
           cluster_cols = F,
           show_rownames = T,
           show_colnames = F,
           annotation_col=ann_col_data,
           annotation_colors=ann_col_colors,
           scale = "none",
           gaps_row = args$gap_rows,
           fontsize_row = 8,
           color = colorRampPalette(cols_scale)(100),
           filename = sprintf("%s/metacell_access_distal_plot_per_gene_minmax.png",args$outdir), width = 10, height = 10)

dat_all_minmax_access_distal  %>%
  pheatmap(cluster_rows = F,
           cluster_cols = F,
           show_rownames = F,
           show_colnames = F,
           annotation_col=ann_col_data,
           annotation_colors=ann_col_colors,
           legend= F,
           scale = "none",
           gaps_row = args$gap_rows,
           fontsize_row = 8,
           color = colorRampPalette(cols_scale)(100),
           filename = sprintf("%s/metacell_access_distal_plot_per_gene_minmax_nolabs.pdf",args$outdir), width = 10, height = 10)

dat_all_minmax_access_distal  %>%
  pheatmap(cluster_rows = F,
           cluster_cols = F,
           show_rownames = F,
           show_colnames = F,
           annotation_col=ann_col_data,
           annotation_colors=ann_col_colors,
           legend= F,
           scale = "none",
           gaps_row = args$gap_rows,
           fontsize_row = 8,
           color = colorRampPalette(cols_scale)(100),
           filename = sprintf("%s/metacell_access_distal_plot_per_gene_minmax_nolabs.png",args$outdir), width = 10, height = 10)


dat_all_minmax_access_TSS <- dat_all_minmax[grepl("atac_TSS",rownames(dat_all_minmax)),]
rownames(dat_all_minmax_access_TSS) <- gsub("_atac_TSS","",rownames(dat_all_minmax_access_TSS))

dat_all_minmax_access_TSS <- dat_all_minmax_access_TSS[order_goi,]

dat_all_minmax_access_TSS  %>%
  pheatmap(cluster_rows = F,
           cluster_cols = F,
           show_rownames = T,
           show_colnames = F,
           annotation_col=ann_col_data,
           annotation_colors=ann_col_colors,
           scale = "none",
           gaps_row = args$gap_rows,
           fontsize_row = 8,
           color = colorRampPalette(cols_scale)(100),
           filename = sprintf("%s/metacell_access_TSS_plot_per_gene_minmax.pdf",args$outdir), width = 10, height = 10)

dat_all_minmax_access_TSS  %>%
  pheatmap(cluster_rows = F,
           cluster_cols = F,
           show_rownames = T,
           show_colnames = F,
           annotation_col=ann_col_data,
           annotation_colors=ann_col_colors,
           scale = "none",
           gaps_row = args$gap_rows,
           fontsize_row = 8,
           color = colorRampPalette(cols_scale)(100),
           filename = sprintf("%s/metacell_access_TSS_plot_per_gene_minmax.png",args$outdir), width = 10, height = 10)

dat_all_minmax_access_TSS  %>%
  pheatmap(cluster_rows = F,
           cluster_cols = F,
           show_rownames = F,
           show_colnames = F,
           annotation_col=ann_col_data,
           annotation_colors=ann_col_colors,
           legend= F,
           scale = "none",
           gaps_row = args$gap_rows,
           fontsize_row = 8,
           color = colorRampPalette(cols_scale)(100),
           filename = sprintf("%s/metacell_access_TSS_plot_per_gene_minmax_nolabs.pdf",args$outdir), width = 10, height = 10)

dat_all_minmax_access_TSS  %>%
  pheatmap(cluster_rows = F,
           cluster_cols = F,
           show_rownames = F,
           show_colnames = F,
           annotation_col=ann_col_data,
           annotation_colors=ann_col_colors,
           legend= F,
           scale = "none",
           gaps_row = args$gap_rows,
           fontsize_row = 8,
           color = colorRampPalette(cols_scale)(100),
           filename = sprintf("%s/metacell_access_TSS_plot_per_gene_minmax_nolabs.png",args$outdir), width = 10, height = 10)

print("Heatmaps done!")

###################
## Data object DE/DA genes ##
###################

print("Data object DE/DA genes")

# Gene list
############

# genes_clusters <- read.table(args$DE_gene_clusters)
# colnames(genes_clusters) <- "Gene"
# ix_genes_clusters <- sort(genes_clusters$Gene,index.return=T)$ix

DEG_overview <- read.table(args$DEG_overview,header=T)
DAG_overview <- read.table(args$DAG_overview,header=T)

DE_genes_clusters <- read.table(args$DE_gene_clusters)
colnames(DE_genes_clusters) <- "Gene"
DE_ix_genes_clusters <- sort(DE_genes_clusters$Gene,index.return=T)$ix

DE_gene_clusters <- merge(DE_genes_clusters,DEG_overview)[,c("Gene","gene_cluster")][sort(DE_ix_genes_clusters,index.return=T)$ix,]
dim(DE_gene_clusters)
head(DE_gene_clusters)

DA_genes_clusters <- read.table(args$DA_gene_clusters)
colnames(DA_genes_clusters) <- "Gene"
DA_ix_genes_clusters <- sort(DA_genes_clusters$Gene,index.return=T)$ix

DA_gene_clusters <- merge(DA_genes_clusters,DAG_overview)[,c("Gene","gene_cluster")][sort(DA_ix_genes_clusters,index.return=T)$ix,]
dim(DE_gene_clusters)
head(DE_gene_clusters)

# Number of genes DA (distal), but not DE
sum(DAG_overview$Gene[!is.na(DAG_overview$gene_cluster)]%in%DEG_overview$Gene[!is.na(DEG_overview$gene_cluster)])
sum(!(DAG_overview$Gene[!is.na(DAG_overview$gene_cluster)]%in%DEG_overview$Gene[!is.na(DEG_overview$gene_cluster)]))
sum(!(DEG_overview$Gene[!is.na(DEG_overview$gene_cluster)]%in%DAG_overview$Gene[!is.na(DAG_overview$gene_cluster)]))

genelist <- unique(c(DAG_overview$Gene[!is.na(DAG_overview$gene_cluster)],DEG_overview$Gene[!is.na(DEG_overview$gene_cluster)]))
genelist <- genelist[(genelist%in%rownames(atac.sce.Distal)) & (genelist%in%rownames(sce))]

# DE_gene_clusters <- merge(genes_clusters,DEG_overview)[,c("Gene","gene_cluster")][sort(ix_genes_clusters,index.return=T)$ix,]
# dim(DE_gene_clusters)
# head(DE_gene_clusters)

# Data obj
############

val_per_TF <- 2

dat_all_full <- NULL
dat_all_full_z <- NULL
dat_all_full_minmax <- NULL
genes_to_remove <- NULL
print(length(genelist))
for (i in 1:length(genelist)){
  if (i%%1000==0){print(i)}
  if ((genelist[i]%in%rownames(atac.sce.Distal)) & (genelist[i]%in%rownames(sce))){
      data_goi <- rbind(assay(sce,"logcounts")[genelist[i],],
                    assay(atac.sce.Distal,"logcounts")[genelist[i],]) %>% data.frame
  }
  colnames(data_goi) <- colnames(sce)

  if (sum(rowSums(data_goi)==0)>0){
    genes_to_remove <- c(genes_to_remove,genelist[i])
    next()
  }

  data_goi_sorted <- data_goi[,pseudotime.dt$cell]
  dat_all_full <- rbind(dat_all_full,data_goi_sorted)

  # z-scores
  dat_all_full_z <- rbind(dat_all_full_z,t(apply(data_goi_sorted, 1, z_score)))

  # minmax normalisation
  dat_all_full_minmax <- rbind(dat_all_full_minmax,t(apply(data_goi_sorted, 1, minmax.10pct.normalisation)))
}
genelist <- genelist[!(genelist%in%genes_to_remove)]
rownames(dat_all_full) <- paste0(rep(genelist,each=val_per_TF),rep("_",length(genelist)*val_per_TF),rep(c("exprs","atac_distal"),length(genelist)))
rownames(dat_all_full_z) <- paste0(rep(genelist,each=val_per_TF),rep("_",length(genelist)*val_per_TF),rep(c("exprs","atac_distal"),length(genelist)))
rownames(dat_all_full_minmax) <- paste0(rep(genelist,each=val_per_TF),rep("_",length(genelist)*val_per_TF),rep(c("exprs","atac_distal"),length(genelist)))


to.plot_all <- data.frame(cbind(dat_all_full[grepl("exprs",rownames(dat_all_full)),],
                        dat_all_full[grepl("atac_distal",rownames(dat_all_full)),]))

to.plot_z <- data.frame(cbind(dat_all_full_z[grepl("exprs",rownames(dat_all_full_z)),],
                        dat_all_full_z[grepl("atac_distal",rownames(dat_all_full_z)),]))

to.plot_minmax <- data.frame(cbind(dat_all_full_minmax[grepl("exprs",rownames(dat_all_full_minmax)),],
                             dat_all_full_minmax[grepl("atac_distal",rownames(dat_all_full_minmax)),]))


#############################
## Lineplots ifo metacells ##
#############################

print("Lineplots")
                             
# Lineplot
############

to.lineplot_minmax <- NULL
for (i in unique(DE_gene_clusters$gene_cluster)){
  to.lineplot_minmax <- rbind(to.lineplot_minmax,
                              cbind(rep(1:ncol(dat_all_full_minmax),3),
                                    c(colMeans(dat_all_full_minmax[grepl("exprs",rownames(dat_all_full_minmax)) &
                                          (gsub("_.*$","",rownames(dat_all_full_minmax))%in%DE_gene_clusters$Gene[DE_gene_clusters$gene_cluster==i]),],na.rm=T),
                                      colMeans(dat_all_full_minmax[grepl("atac_TSS",rownames(dat_all_full_minmax)) &
                                          (gsub("_.*$","",rownames(dat_all_full_minmax))%in%DE_gene_clusters$Gene[DE_gene_clusters$gene_cluster==i]),],na.rm=T),
                                      colMeans(dat_all_full_minmax[grepl("atac_distal",rownames(dat_all_full_minmax)) &
                                          (gsub("_.*$","",rownames(dat_all_full_minmax))%in%DE_gene_clusters$Gene[DE_gene_clusters$gene_cluster==i]),],na.rm=T)),
                                    c(rep("exprs",ncol(dat_all_full_minmax)),rep("atac_TSS",ncol(dat_all_full_minmax)),rep("atac_distal",ncol(dat_all_full_minmax))),
                                    rep(i,ncol(dat_all_full_minmax)*3)))
}
colnames(to.lineplot_minmax) <- c("x_val","av_val","modality","gene_cluster")
to.lineplot_minmax <- data.frame(to.lineplot_minmax)
to.lineplot_minmax$x_val <- as.integer(to.lineplot_minmax$x_val)
to.lineplot_minmax$av_val <- as.numeric(to.lineplot_minmax$av_val)
to.lineplot_minmax$modality <- as.factor(to.lineplot_minmax$modality)
to.lineplot_minmax$gene_cluster <- as.factor(to.lineplot_minmax$gene_cluster)

p <- ggplot(data=to.lineplot_minmax, aes(x=x_val, y=av_val, color=modality)) +
    geom_smooth(span=0.05) +
    ylim(0,1) +
    facet_wrap(~gene_cluster,  ncol=2)

pdf(sprintf("%s/metacell_lineplot_per_gene_cluster_minmax.pdf",args$outdir), width=7, height=5)
print(p)
dev.off()

# Lineplot RNA
############

to.lineplot_RNA_minmax <- NULL
for (i in unique(DE_gene_clusters$gene_cluster)){
  to.lineplot_RNA_minmax <- rbind(to.lineplot_RNA_minmax,
                              cbind(rep(1:ncol(dat_all_full_minmax)),
                                    c(colMeans(dat_all_full_minmax[grepl("exprs",rownames(dat_all_full_minmax)) &
                                          (gsub("_.*$","",rownames(dat_all_full_minmax))%in%DE_gene_clusters$Gene[DE_gene_clusters$gene_cluster==i]),],na.rm=T)),
                                    c(rep("exprs",ncol(dat_all_full_minmax))),
                                    rep(i,ncol(dat_all_full_minmax))))
}
colnames(to.lineplot_RNA_minmax) <- c("x_val","av_val","modality","gene_cluster")
to.lineplot_RNA_minmax <- data.frame(to.lineplot_RNA_minmax)
to.lineplot_RNA_minmax$x_val <- as.integer(to.lineplot_RNA_minmax$x_val)
to.lineplot_RNA_minmax$av_val <- as.numeric(to.lineplot_RNA_minmax$av_val)
to.lineplot_RNA_minmax$modality <- as.factor(to.lineplot_RNA_minmax$modality)
to.lineplot_RNA_minmax$gene_cluster <- as.factor(to.lineplot_RNA_minmax$gene_cluster)

p <- ggplot(data=to.lineplot_RNA_minmax, aes(x=x_val, y=av_val,color=1)) +
    geom_smooth(span=0.05) +
    ylim(0,1) +
    facet_wrap(~gene_cluster,  ncol=2)

output_plot(p,sprintf("%s/metacell_exprs_lineplot_per_gene_cluster_minmax",args$outdir), width=7, height=5)

# plot per cluster
dir.create(sprintf("%s/plots_per_gene_cluster",args$outdir))

for (cluster in unique(DE_gene_clusters$gene_cluster)){
p <- ggplot(data=to.lineplot_RNA_minmax[to.lineplot_RNA_minmax$gene_cluster==cluster,], aes(x=x_val, y=av_val,color=1)) +
    geom_smooth(span=0.05) +
    ylim(0,1) +  
    theme_classic() +
    theme(
      axis.text =  element_text(size=rel(0.8)),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.title = element_blank(),
      legend.position = "top",
      legend.text = element_text(size=rel(0.75))
    )

    output_plot(p,sprintf("%s/plots_per_gene_cluster/metacell_exprs_lineplot_per_gene_cluster%s_minmax",args$outdir,cluster), width=7, height=5)
}

# Lineplot ATAC
############

to.lineplot_ATAC_minmax <- NULL
for (i in unique(DA_gene_clusters$gene_cluster)){
  to.lineplot_ATAC_minmax <- rbind(to.lineplot_ATAC_minmax,
                              cbind(rep(1:ncol(dat_all_full_minmax)),
                                    c(colMeans(dat_all_full_minmax[grepl("atac_distal",rownames(dat_all_full_minmax)) &
                                          (gsub("_.*$","",rownames(dat_all_full_minmax))%in%DA_gene_clusters$Gene[DA_gene_clusters$gene_cluster==i]),],na.rm=T)),
                                    c(rep("exprs",ncol(dat_all_full_minmax))),
                                    rep(i,ncol(dat_all_full_minmax))))
}
colnames(to.lineplot_ATAC_minmax) <- c("x_val","av_val","modality","gene_cluster")
to.lineplot_ATAC_minmax <- data.frame(to.lineplot_ATAC_minmax)
to.lineplot_ATAC_minmax$x_val <- as.integer(to.lineplot_ATAC_minmax$x_val)
to.lineplot_ATAC_minmax$av_val <- as.numeric(to.lineplot_ATAC_minmax$av_val)
to.lineplot_ATAC_minmax$modality <- as.factor(to.lineplot_ATAC_minmax$modality)
to.lineplot_ATAC_minmax$gene_cluster <- as.factor(to.lineplot_ATAC_minmax$gene_cluster)

p <- ggplot(data=to.lineplot_ATAC_minmax, aes(x=x_val, y=av_val,color=1)) +
    geom_smooth(span=0.05) +
    ylim(0,1) +
    facet_wrap(~gene_cluster,  ncol=2)

output_plot(p,sprintf("%s/metacell_access_lineplot_per_gene_cluster_minmax",args$outdir), width=7, height=5)

# plot per cluster
for (cluster in unique(DA_gene_clusters$gene_cluster)){
p <- ggplot(data=to.lineplot_ATAC_minmax[to.lineplot_ATAC_minmax$gene_cluster==cluster,], aes(x=x_val, y=av_val,color=1)) +
    geom_smooth(span=0.05) +
    ylim(0,1) +  
    theme_classic() +
    theme(
      axis.text =  element_text(size=rel(0.8)),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.title = element_blank(),
      legend.position = "top",
      legend.text = element_text(size=rel(0.75))
    )

    output_plot(p,sprintf("%s/plots_per_gene_cluster/metacell_access_lineplot_per_gene_cluster%s_minmax",args$outdir,cluster), width=7, height=5)
}


# ##########
# ## UMAP ##
# ##########

# # umap coord
# umap.dt <- fread(args$umap_coord)
# dim(umap.dt)
# umap.dt <- umap.dt[umap.dt$sample%in%pseudotime.dt$cell,]
# umap.dt <- umap.dt[sort(sort(colnames(dat_all),index.return=T)$ix,index.return=T)$ix,]
# dim(umap.dt)
# head(umap.dt)
# sum(umap.dt$sample==colnames(dat_all))

# # add expression & atac data
# dat_all_merge <- data.frame(t(dat_all)) %>% rownames_to_column(var="sample")
# dim(umap.dt)
# to.plot <- merge(umap.dt,dat_all_merge)
# dim(to.plot)

# dir.create(sprintf("%s/metacell_plot_per_gene",args$outdir))

# for (mod in c("exprs","atac_TSS","atac_distal")){
#   for (i in 1:length(args$goi)){

#     p <- ggplot(to.plot, aes(x=V1, y=V2, fill=get(paste(args$goi[i],mod,sep="_")))) +
#       geom_point(size=3, shape=21,color='transparent') +
#       scale_fill_viridis_c(option="inferno") +
#       ggtitle(paste0(args$goi[i]," - ",toupper(mod))) +
#       labs(fill = mod) +
#       theme_classic() +
#       ggplot_theme_NoAxes()

#     pdf(sprintf("%s/metacell_plot_per_gene/%s_%s_umap.pdf",args$outdir,args$goi[i],mod), width=7, height=5)
#     print(p)
#     dev.off()
#   }
# }

# # add expression & atac data - minmax
# dat_all_minmax_merge <- data.frame(t(dat_all_minmax)) %>% rownames_to_column(var="sample")
# dim(umap.dt)
# to.plot_minmax <- merge(umap.dt,dat_all_minmax_merge)
# dim(to.plot_minmax)

# dir.create(sprintf("%s/metacell_plot_per_gene_minmax",args$outdir))

# for (mod in c("exprs","atac_TSS","atac_distal")){
#   for (i in 1:length(args$goi)){

#     p <- ggplot(to.plot_minmax, aes(x=V1, y=V2, fill=get(paste(args$goi[i],mod,sep="_")))) +
#       geom_point(size=3, shape=21,color='transparent') +
#       scale_fill_viridis_c(option="inferno") +
#       ggtitle(paste0(args$goi[i]," - ",toupper(mod))) +
#       labs(fill = mod) +
#       theme_classic() +
#       ggplot_theme_NoAxes()

#     pdf(sprintf("%s/metacell_plot_per_gene_minmax/%s_%s_umap.pdf",args$outdir,args$goi[i],mod), width=7, height=5)
#     print(p)
#     dev.off()
#   }
# }


# # All genes
# ################

# # make sure clustering object has no NA values
# Rclusterpp.hclust(to.plot_z) -> clustering.genes

# n_clusters <- 10

# cutree(clustering.genes, k=n_clusters) %>%
#   as_tibble(rownames = "Gene") %>%
#   dplyr::rename(cluster = value) -> clusters.n

# lapply(c(1:n_clusters), count.gaps, clusters.n) %>%
#   do.call("rbind",.) %>%
#   as.vector -> clusters.gaps.n

# extensions <- c("pdf","png")
# for (ext in extensions){
#   to.plot_minmax %>%
#     rownames_to_column(var="Gene") %>%
#     full_join(clusters.n) %>%
#     arrange(cluster) %>%
#     column_to_rownames(var = "Gene") %>%
#     select(-cluster) %>%
#     pheatmap(cluster_rows = F,
#             cluster_cols = F,
#             show_rownames = F,
#             show_colnames = F,
#             scale = "none",
#             fontsize_row = 8,
#             gaps_row = clusters.gaps.n,
#             gaps_col = c(1:(val_per_TF-1))*ncol(sce),
#             color = colorRampPalette(cols_scale)(100),
#             filename = sprintf("%s/DEorDA_heatmap_%scluster_minmax.%s",args$outdir,n_clusters,ext), width = 10, height = 10)
# }

# # Split data
# ############

# DE_genes_clusters <- read.table(args$DE_gene_clusters)
# colnames(DE_genes_clusters) <- "Gene"
# DE_ix_genes_clusters <- sort(DE_genes_clusters$Gene,index.return=T)$ix

# DE_gene_clusters <- merge(DE_genes_clusters,DEG_overview)[,c("Gene","gene_cluster")][sort(DE_ix_genes_clusters,index.return=T)$ix,]
# dim(DE_gene_clusters)
# head(DE_gene_clusters)

# DA_genes_clusters <- read.table(args$DA_gene_clusters)
# colnames(DA_genes_clusters) <- "Gene"
# DA_ix_genes_clusters <- sort(DA_genes_clusters$Gene,index.return=T)$ix

# DA_gene_clusters <- merge(DA_genes_clusters,DAG_overview)[,c("Gene","gene_cluster")][sort(DA_ix_genes_clusters,index.return=T)$ix,]
# dim(DE_gene_clusters)
# head(DE_gene_clusters)

# genelist_up <- unique(c(DE_gene_clusters$Gene[DE_gene_clusters$gene_cluster%in%c("E","F")],
#                         DA_gene_clusters$Gene[DA_gene_clusters$gene_cluster%in%c("E","F")]))
# length(genelist_up)
# genelist_down <- unique(c(DE_gene_clusters$Gene[DE_gene_clusters$gene_cluster%in%c("A","B","C")],
#                           DA_gene_clusters$Gene[DA_gene_clusters$gene_cluster%in%c("A","B","C")]))
# length(genelist_down)

# to.plot_z_up <- to.plot_z[gsub("_exprs","",rownames(to.plot_z))%in%genelist_up,]
# to.plot_z_down <- to.plot_z[gsub("_exprs","",rownames(to.plot_z))%in%genelist_down,]

# to.plot_minmax_up <- to.plot_minmax[gsub("_exprs","",rownames(to.plot_minmax))%in%genelist_up,]
# to.plot_minmax_down <- to.plot_minmax[gsub("_exprs","",rownames(to.plot_minmax))%in%genelist_down,]


# # Downreg & Upreg
# ############

# # make sure clustering object has no NA values
# Rclusterpp.hclust(to.plot_z_up) -> clustering.genes

# cutree(clustering.genes, k=n_clusters) %>%
#   as_tibble(rownames = "Gene") %>%
#   dplyr::rename(cluster = value) -> clusters.n

# lapply(c(1:n_clusters), count.gaps, clusters.n) %>%
#   do.call("rbind",.) %>%
#   as.vector -> clusters.gaps.n

# extensions <- c("pdf","png")
# for (ext in extensions){
#   to.plot_minmax_up %>%
#     rownames_to_column(var="Gene") %>%
#     full_join(clusters.n) %>%
#     arrange(cluster) %>%
#     column_to_rownames(var = "Gene") %>%
#     select(-cluster) %>%
#     pheatmap(cluster_rows = F,
#             cluster_cols = F,
#             show_rownames = F,
#             show_colnames = F,
#             scale = "none",
#             fontsize_row = 8,
#             gaps_col = c(1:(val_per_TF-1))*ncol(sce),
#             gaps_row = clusters.gaps.n,
#             color = colorRampPalette(cols_scale)(100),
#             filename = sprintf("%s/DEorDA_upreg_heatmap_%scluster_minmax.%s",args$outdir,n_clusters,ext), width = 10, height = 10)
# }


# # make sure clustering object has no NA values
# Rclusterpp.hclust(to.plot_z_down) -> clustering.genes

# cutree(clustering.genes, k=n_clusters) %>%
#   as_tibble(rownames = "Gene") %>%
#   dplyr::rename(cluster = value) -> clusters.n

# lapply(c(1:n_clusters), count.gaps, clusters.n) %>%
#   do.call("rbind",.) %>%
#   as.vector -> clusters.gaps.n

# extensions <- c("pdf","png")
# for (ext in extensions){
#   to.plot_minmax_down %>%
#   rownames_to_column(var="Gene") %>%
#   full_join(clusters.n) %>%
#   arrange(cluster) %>%
#   column_to_rownames(var = "Gene") %>%
#   select(-cluster) %>%
#   pheatmap(cluster_rows = F,
#            cluster_cols = F,
#            show_rownames = F,
#            show_colnames = F,
#            scale = "none",
#            fontsize_row = 8,
#            gaps_row = clusters.gaps.n,
#            gaps_col = c(1:(val_per_TF-1))*ncol(sce),
#            color = colorRampPalette(cols_scale)(100),
#            filename = sprintf("%s/DEorDA_downreg_heatmap_%scluster_minmax.%s",args$outdir,n_clusters,ext), width = 10, height = 10)

# }

# DE and DA
############

geneset_oi <- "DEandDA"
dir.create(sprintf("%s/%s",args$outdir,geneset_oi))

genelist_overlap <- DAG_overview$Gene[!is.na(DAG_overview$gene_cluster)][DAG_overview$Gene[!is.na(DAG_overview$gene_cluster)]%in%DEG_overview$Gene[!is.na(DEG_overview$gene_cluster)]]
length(genelist_overlap)

to.plot_z_overlap <- to.plot_z[gsub("_exprs","",rownames(to.plot_z))%in%genelist_overlap,]
dim(to.plot_z_overlap)

to.plot_minmax_overlap <- to.plot_minmax[gsub("_exprs","",rownames(to.plot_minmax))%in%genelist_overlap,]

Rclusterpp.hclust(to.plot_z_overlap) -> clustering.genes

extensions <- c("pdf","png")
for (n_clusters in 21:23){
  cutree(clustering.genes, k=n_clusters) %>%
    as_tibble(rownames = "Gene") %>%
    dplyr::rename(cluster = value) -> clusters.n

  lapply(c(1:n_clusters), count.gaps, clusters.n) %>%
    do.call("rbind",.) %>%
    as.vector -> clusters.gaps.n

  for (ext in extensions){
    to.plot_minmax_overlap %>%
      rownames_to_column(var="Gene") %>%
      full_join(clusters.n) %>%
      arrange(cluster) %>%
      column_to_rownames(var = "Gene") %>%
      select(-cluster) %>%
      pheatmap(cluster_rows = F,
              cluster_cols = F,
              show_rownames = F,
              show_colnames = F,
              scale = "none",
              fontsize_row = 8,
              gaps_col = c(1:(val_per_TF-1))*ncol(sce),
              gaps_row = clusters.gaps.n,
              color = colorRampPalette(cols_scale)(100),
              filename = sprintf("%s/%s/heatmap_%sclusters_minmax.%s",args$outdir,geneset_oi,n_clusters,ext), width = 12, height = 20)
  }
}

# 21 clusters
n_clusters <- 21
dir.create(sprintf("%s/%s/%sclusters",args$outdir,geneset_oi,n_clusters))

cutree(clustering.genes, k=n_clusters) %>%
  as_tibble(rownames = "Gene") %>%
  dplyr::rename(cluster = value) -> clusters.n

# order 
req_order <- c(7,3,6,9,14,16,10,2,12,17,5,15,18,4,8,11,13,21,1,20,19)
#req_order <- c(18,8,2,11,3,12,19,1,4,7,16,9,15,5,13,6,14,10,21,20,17)
print(cbind(clusters.n$cluster,req_order[clusters.n$cluster])[!duplicated(clusters.n$cluster),])
clusters.n$cluster <- req_order[clusters.n$cluster]
print(table(clusters.n$cluster))

# reorder to remove small clusters
req_order <- 1:length(req_order) 
clust_to_remove <- names(table(clusters.n$cluster))[table(clusters.n$cluster)<10]
print(clust_to_remove)
req_order[as.numeric(clust_to_remove)] <- 100 + req_order[as.numeric(clust_to_remove)] 
req_order <- sort(req_order,index.return=T)$ix
print(req_order)
req_order <- order(req_order)
print(req_order)
print(cbind(clusters.n$cluster,req_order[clusters.n$cluster])[!duplicated(clusters.n$cluster),])
clusters.n$cluster <- req_order[clusters.n$cluster]
print(table(clusters.n$cluster))

lapply(c(1:n_clusters), count.gaps, clusters.n) %>%
  do.call("rbind",.) %>%
  as.vector -> clusters.gaps.n

dat_21clust <- to.plot_minmax_overlap %>%
  rownames_to_column(var="Gene") %>%
    full_join(clusters.n) %>%
      arrange(cluster) %>%
      column_to_rownames(var = "Gene")

for (ext in extensions){
  # heatmap
  to.plot_minmax_overlap %>%
    rownames_to_column(var="Gene") %>%
    full_join(clusters.n) %>%
    arrange(cluster) %>%
    column_to_rownames(var = "Gene") %>%
    select(-cluster) %>%
    pheatmap(cluster_rows = F,
        cluster_cols = F,
        show_rownames = F,
        show_colnames = F,
        scale = "none",
        fontsize_row = 8,
        gaps_col = c(1:(val_per_TF-1))*ncol(sce),
        gaps_row = clusters.gaps.n,
        color = colorRampPalette(cols_scale)(100),
        filename = sprintf("%s/%s/%sclusters/heatmap_%sclusters_minmax_sorted.%s",args$outdir,geneset_oi,n_clusters,n_clusters,ext), width = 12, height = 20)

  # heatmap - no labs
  to.plot_minmax_overlap %>%
      rownames_to_column(var="Gene") %>%
      full_join(clusters.n) %>%
      arrange(cluster) %>%
      column_to_rownames(var = "Gene") %>%
      select(-cluster) %>%
      pheatmap(cluster_rows = F,
          cluster_cols = F,
          show_rownames = F,
          show_colnames = F,
          legend= F,
          scale = "none",
          fontsize_row = 8,
          gaps_col = c(1:(val_per_TF-1))*ncol(sce),
          gaps_row = clusters.gaps.n,
          color = colorRampPalette(cols_scale)(100),
          filename = sprintf("%s/%s/%sclusters/heatmap_%sclusters_minmax_sorted_nolabs.%s",args$outdir,geneset_oi,n_clusters,n_clusters,ext), width = 12, height = 20)

  to.plot_minmax_overlap_nosmall <- to.plot_minmax_overlap %>%
      rownames_to_column(var="Gene") %>%
      full_join(clusters.n) %>%
      arrange(cluster)

  clust_to_remove <- names(table(clusters.n$cluster))[table(clusters.n$cluster)<10]
  to.plot_minmax_overlap_nosmall <- to.plot_minmax_overlap_nosmall[!(to.plot_minmax_overlap_nosmall$cluster%in%clust_to_remove),]

  print(clusters.gaps.n)
  clusters.gaps.n_nosmall <- clusters.gaps.n
  len <- length(clusters.gaps.n_nosmall)
  for (i in 1:length(clust_to_remove)){
    clust <-as.numeric(clust_to_remove[i])
    if (clust == len){
      clusters.gaps.n_nosmall[clust] <- NA
    } else {
      clusters.gaps.n_nosmall[(clust+1):len] <- clusters.gaps.n_nosmall[(clust+1):len] - table(clusters.n$cluster)[clust_to_remove[i]]
      clusters.gaps.n_nosmall[clust] <- NA
    }
  }
  print(clusters.gaps.n_nosmall)
  clusters.gaps.n_nosmall <- clusters.gaps.n_nosmall[!is.na(clusters.gaps.n_nosmall)]
  print(clusters.gaps.n_nosmall)

  # heatmap - no small clusters
  to.plot_minmax_overlap_nosmall %>%
    arrange(cluster) %>%
    column_to_rownames(var = "Gene") %>%
    .[!(cluster%in%names(table(clusters.n$cluster))[table(clusters.n$cluster)<10])] %>%
    select(-cluster) %>%
    pheatmap(cluster_rows = F,
        cluster_cols = F,
        show_rownames = F,
        show_colnames = F,
        scale = "none",
        fontsize_row = 8,
        gaps_col = c(1:(val_per_TF-1))*ncol(sce),
        gaps_row = clusters.gaps.n_nosmall,
        color = colorRampPalette(cols_scale)(100),
        filename = sprintf("%s/%s/%sclusters/heatmap_%sclusters_minmax_sorted_nosmall.%s",args$outdir,geneset_oi,n_clusters,n_clusters,ext), width = 12, height = 20)

  # heatmap - no small clusters - no labs
  to.plot_minmax_overlap_nosmall %>%
      arrange(cluster) %>%
      column_to_rownames(var = "Gene") %>%
      .[!(cluster%in%names(table(clusters.n$cluster))[table(clusters.n$cluster)<10]) ] %>%
      select(-cluster) %>%
      pheatmap(cluster_rows = F,
          cluster_cols = F,
          show_rownames = F,
          show_colnames = F,
          legend= F,
          scale = "none",
          fontsize_row = 8,
          gaps_col = c(1:(val_per_TF-1))*ncol(sce),
          gaps_row = clusters.gaps.n_nosmall,
          color = colorRampPalette(cols_scale)(100),
          filename = sprintf("%s/%s/%sclusters/heatmap_%sclusters_minmax_sorted_nosmall_nolabs.%s",args$outdir,geneset_oi,n_clusters,n_clusters,ext), width = 12, height = 20)
}

av_per_cluster <- NULL
for (i in 1:n_clusters){
  av_per_cluster <- rbind(av_per_cluster,
                          dat_21clust[dat_21clust$cluster==i,] %>% select(-cluster) %>% colMeans)
}
rownames(av_per_cluster) <- 1:n_clusters
colnames(av_per_cluster) <- paste0(rep(1:(ncol(av_per_cluster)/2),2),rep(c("_RNA","_ATAC"),each=ncol(av_per_cluster)/2))
av_per_cluster <- av_per_cluster %>% melt()
colnames(av_per_cluster) <- c("cluster","metacell","value")
av_per_cluster$modality <- gsub("^[0-9]+_","",av_per_cluster$metacell)
av_per_cluster$metacell <- gsub("_[A-Z]+$","",av_per_cluster$metacell)
av_per_cluster$cluster <- as.factor(av_per_cluster$cluster)
av_per_cluster$metacell <- as.numeric(av_per_cluster$metacell)
av_per_cluster$modality <- as.factor(av_per_cluster$modality)
head(av_per_cluster)

p <- ggplot(data=av_per_cluster, aes(x=metacell, y=value, color=modality)) +
    geom_smooth(span=0.05, size=0.5) +
    ylim(0,1) +
    facet_wrap(~cluster,  nrow=21,ncol=1) +  
    theme_classic() +
    theme(
      axis.text =  element_text(size=rel(0.35)),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.title = element_blank(),
      legend.position = "none",
      strip.text = element_text(size=rel(0.35))
    )

output_plot(p,sprintf("%s/%s/%sclusters/lineplot_%sclusters_minmax_sorted",args$outdir,geneset_oi,n_clusters,n_clusters), width=1.5, height=20)

av_per_cluster_nosmall <- av_per_cluster[!(av_per_cluster$cluster%in%names(table(clusters.n$cluster))[table(clusters.n$cluster)<10]),]
p <- ggplot(data=av_per_cluster_nosmall, aes(x=metacell, y=value, color=modality)) +
    geom_smooth(span=0.05, size=0.5) +
    ylim(0,1) +
    facet_wrap(~cluster,  nrow=21,ncol=1) +  
    theme_classic() +
    theme(
      axis.text =  element_text(size=rel(0.35)),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.title = element_blank(),
      legend.position = "none",
      strip.text = element_text(size=rel(0.35))
    )

output_plot(p,sprintf("%s/%s/%sclusters/lineplot_%sclusters_minmax_sorted_nosmall",args$outdir,geneset_oi,n_clusters,n_clusters), width=1.5, height=20)


# pdf(sprintf("%s/metacell_plot_per_gene_%scluster_minmax_lineplot_av_DEandDA.pdf",args$outdir,n_clusters), width=7, height=5)
# plot(av_per_cluster$metacell[av_per_cluster$cluster==1 & av_per_cluster$modality=="RNA"],av_per_cluster$value[av_per_cluster$cluster==1 & av_per_cluster$modality=="RNA"],type="l")
# dev.off()

# # Table
table_21clusters <- cbind(gsub("_exprs","",rownames(dat_21clust)),dat_21clust$cluster)
table_21clusters <- data.frame(table_21clusters)
colnames(table_21clusters) <- c("Gene","Cluster")

aggregate_bins <- function(values,bins){
  agg_val <- NULL
  for (i in 1:(length(bins)-1)){
    if (i!=(length(bins))-1){
      agg_val <- c(agg_val,mean(values[bins[i]:(bins[i+1]-1)]))
    }
    else {
      agg_val <- c(agg_val,mean(values[bins[i]:bins[i+1]]))
    }
  }
  return(agg_val)
}


table_dat_RNA <- NULL
table_dat_ATAC <- NULL
for (i in 1:nrow(table_21clusters)){
  table_dat_RNA <- rbind(as.numeric(dat_all_full[grepl(paste0("^",table_21clusters$Gene[i],"_"),rownames(dat_all_full)) & grepl("exprs",rownames(dat_all_full)),]))
  table_dat_ATAC <- rbind(as.numeric(dat_all_full[grepl(paste0("^",table_21clusters$Gene[i],"_"),rownames(dat_all_full)) & grepl("atac",rownames(dat_all_full)),]))
}
table_21clusters_RNA <- cbind(table_21clusters,table_dat_RNA)
table_21clusters_ATAC <- cbind(table_21clusters,table_dat_ATAC)
colnames(table_21clusters_RNA) <- c("Gene","Cluster",paste0(rep("exprs_mc_",ncol(table_dat_RNA)),1:ncol(table_dat_RNA)))
colnames(table_21clusters_ATAC) <- c("Gene","Cluster",paste0(rep("atac_mc_",ncol(table_dat_ATAC)),1:ncol(table_dat_ATAC)))
write.table(table_21clusters_RNA,file=sprintf("%s/%s/%sclusters/table_%sclusters_DEandDA_RNA.txt",args$outdir,geneset_oi,n_clusters,n_clusters),col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
write.table(table_21clusters_ATAC,file=sprintf("%s/%s/%sclusters/table_%sclusters_DEandDA_ATAC.txt",args$outdir,geneset_oi,n_clusters,n_clusters),col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")


# bins <- c(seq(1,261,18),261)
# bins[1:(length(bins)-1)]-bins[2:length(bins)]
# table_dat_RNA <- NULL
# table_dat_ATAC <- NULL
# for (i in 1:nrow(table_21clusters)){
#   table_dat_RNA <- rbind(table_dat_RNA,aggregate_bins(as.numeric(dat_all_full[grepl(paste0("^",table_21clusters$Gene[i],"_"),rownames(dat_all_full)) & grepl("exprs",rownames(dat_all_full)),]),bins))
#   table_dat_ATAC <- rbind(table_dat_ATAC,aggregate_bins(as.numeric(dat_all_full[grepl(paste0("^",table_21clusters$Gene[i],"_"),rownames(dat_all_full)) & grepl("atac",rownames(dat_all_full)),]),bins))
# }
# table_21clusters <- cbind(table_21clusters,table_dat_RNA,table_dat_ATAC)
# colnames(table_21clusters) <- c("Gene","Cluster",paste0(rep("exprs_mc_bin_",ncol(table_dat_RNA)),1:ncol(table_dat_RNA)),
#                                                  paste0(rep("atac_mc_bin_",ncol(table_dat_ATAC)),1:ncol(table_dat_ATAC)))
# write.table(table_21clusters,file=sprintf("%s/%s/%sclusters/table_%sclusters_DEandDA.txt",args$outdir,geneset_oi,n_clusters,n_clusters),col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")


# # DE and not DA
# ############

# geneset_oi <- "DEandnotDA"
# dir.create(sprintf("%s/%s",args$outdir,geneset_oi))

# genelist_overlap <- DEG_overview$Gene[!is.na(DEG_overview$gene_cluster)][!(DEG_overview$Gene[!is.na(DEG_overview$gene_cluster)]%in%DAG_overview$Gene[!is.na(DAG_overview$gene_cluster)])]
# length(genelist_overlap)

# to.plot_z_overlap <- to.plot_z[gsub("_exprs","",rownames(to.plot_z))%in%genelist_overlap,]
# dim(to.plot_z_overlap)

# to.plot_minmax_overlap <- to.plot_minmax[gsub("_exprs","",rownames(to.plot_minmax))%in%genelist_overlap,]

# Rclusterpp.hclust(to.plot_z_overlap) -> clustering.genes

# extensions <- c("pdf","png")
# for (n_clusters in c(5,10,15:25,30)){
#   cutree(clustering.genes, k=n_clusters) %>%
#     as_tibble(rownames = "Gene") %>%
#     dplyr::rename(cluster = value) -> clusters.n

#   lapply(c(1:n_clusters), count.gaps, clusters.n) %>%
#     do.call("rbind",.) %>%
#     as.vector -> clusters.gaps.n

#   for (ext in extensions){
#     to.plot_minmax_overlap %>%
#       rownames_to_column(var="Gene") %>%
#       full_join(clusters.n) %>%
#       arrange(cluster) %>%
#       column_to_rownames(var = "Gene") %>%
#       select(-cluster) %>%
#       pheatmap(cluster_rows = F,
#               cluster_cols = F,
#               show_rownames = F,
#               show_colnames = F,
#               scale = "none",
#               fontsize_row = 8,
#               gaps_col = c(1:(val_per_TF-1))*ncol(sce),
#               gaps_row = clusters.gaps.n,
#               color = colorRampPalette(cols_scale)(100),
#               filename = sprintf("%s/%s/heatmap_%sclusters_minmax.%s",args$outdir,geneset_oi,n_clusters,ext), width = 10, height = 10)
#   }
# }


# # DA and not DE
# ############

# geneset_oi <- "DAandnotDE"
# dir.create(sprintf("%s/%s",args$outdir,geneset_oi))

# genelist_overlap <- DAG_overview$Gene[!is.na(DAG_overview$gene_cluster)][!(DAG_overview$Gene[!is.na(DAG_overview$gene_cluster)]%in%DEG_overview$Gene[!is.na(DEG_overview$gene_cluster)])]
# length(genelist_overlap)

# to.plot_z_overlap <- to.plot_z[gsub("_exprs","",rownames(to.plot_z))%in%genelist_overlap,]
# dim(to.plot_z_overlap)

# to.plot_minmax_overlap <- to.plot_minmax[gsub("_exprs","",rownames(to.plot_minmax))%in%genelist_overlap,]

# Rclusterpp.hclust(to.plot_z_overlap) -> clustering.genes

# extensions <- c("pdf","png")
# for (n_clusters in c(5,10,15:25,30)){
#   cutree(clustering.genes, k=n_clusters) %>%
#     as_tibble(rownames = "Gene") %>%
#     dplyr::rename(cluster = value) -> clusters.n

#   lapply(c(1:n_clusters), count.gaps, clusters.n) %>%
#     do.call("rbind",.) %>%
#     as.vector -> clusters.gaps.n

#   for (ext in extensions){
#     to.plot_minmax_overlap %>%
#       rownames_to_column(var="Gene") %>%
#       full_join(clusters.n) %>%
#       arrange(cluster) %>%
#       column_to_rownames(var = "Gene") %>%
#       select(-cluster) %>%
#       pheatmap(cluster_rows = F,
#               cluster_cols = F,
#               show_rownames = F,
#               show_colnames = F,
#               scale = "none",
#               fontsize_row = 8,
#               gaps_col = c(1:(val_per_TF-1))*ncol(sce),
#               gaps_row = clusters.gaps.n,
#               color = colorRampPalette(cols_scale)(100),
#               filename = sprintf("%s/%s/heatmap_%sclusters_minmax.%s",args$outdir,geneset_oi,n_clusters,ext), width = 10, height = 10)
#   }
# }


print("Reimplement part with GRN coef!")

# # GSA clusters
# ############
# geneset_oi <- "DEandDA"
# n_clusters <- 21

# dim(dat_21clust)
# dat_21clust[1:6,1:6]
# head(dat_21clust[,ncol(dat_21clust)])

# GRN_coef.dt <- fread(args$grn_coef) %>%
#   .[,gene:=toupper(gene)] %>% # .[gene%in%unique(marker_TFs_filt.dt$gene) & tf%in%unique(marker_TFs_filt.dt$gene)] %>%
#   .[pvalue<args$max_pval & abs(beta)>=args$min_coef]

# head(GRN_coef.dt)

# cluster_oi <- c(13,16)

# for (i in 1:length(cluster_oi)){
#   cluster <- cluster_oi[i]
#   dir.create(sprintf("%s/%s/%sclusters/cluster%s",args$outdir,geneset_oi,n_clusters,cluster))

#   dir.create(sprintf("%s/%s/%sclusters/cluster%s/GSA_ind_cluster",args$outdir,geneset_oi,n_clusters,cluster))

#   # Reference genes
#   GenesRefFileName <- sprintf("%s/%s/%sclusters/cluster%s/GSA_ind_cluster/ORA_ref.txt",args$outdir,geneset_oi,n_clusters,cluster)
#   GenesRef <- gsub("_exprs$","",rownames(dat_21clust))
#   write.table(GenesRef,GenesRefFileName,
#               sep="\n",col.names = F,row.names = F,quote = F)

#   # Reference genes
#   GenesORAFileName <- sprintf("%s/%s/%sclusters/cluster%s/GSA_ind_cluster/ORA_goi.txt",args$outdir,geneset_oi,n_clusters,cluster)
#   GenesORA <-  gsub("_exprs$","",rownames(dat_21clust))[dat_21clust[,ncol(dat_21clust)]==cluster]
#   write.table(GenesORA,GenesORAFileName,
#               sep="\n",col.names = F,row.names = F,quote = F)

#   dir_name <- sprintf("%s/%s/%sclusters/cluster%s/GSA_ind_cluster",args$outdir,geneset_oi,n_clusters,cluster)

#   enrichResult_MF <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
#                                  enrichDatabase="geneontology_Molecular_Function_noRedundant",
#                                  interestGeneFile=GenesORAFileName, interestGeneType="genesymbol",
#                                  referenceGeneFile=GenesRefFileName, referenceGeneType="genesymbol",
#                                  dagColor="continuous",
#                                  sigMethod="top", minNum=5, reportNum=500,
#                                  projectName="ORA_MF", outputDirectory=dir_name)

#   enrichResult_BP <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
#                                  enrichDatabase="geneontology_Biological_Process_noRedundant",
#                                  interestGeneFile=GenesORAFileName, interestGeneType="genesymbol",
#                                  referenceGeneFile=GenesRefFileName, referenceGeneType="genesymbol",
#                                  dagColor = "continuous",
#                                  sigMethod="top", minNum=5, reportNum = 500,
#                                  outputDirectory =dir_name, projectName = "ORA_BP")

#   enrichResult_CC <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
#                                  enrichDatabase="geneontology_Cellular_Component_noRedundant",
#                                  interestGeneFile=GenesORAFileName, interestGeneType="genesymbol",
#                                  referenceGeneFile=GenesRefFileName, referenceGeneType="genesymbol",
#                                  dagColor = "continuous",
#                                  sigMethod="top", minNum=5, reportNum = 500,
#                                  outputDirectory =dir_name, projectName = "ORA_CC")

#   enrichResult_KEGG <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
#                                    enrichDatabase="pathway_KEGG",
#                                    interestGeneFile=GenesORAFileName, interestGeneType="genesymbol",
#                                    referenceGeneFile=GenesRefFileName, referenceGeneType="genesymbol",
#                                    dagColor = "continuous",
#                                    sigMethod="top", minNum=5, reportNum = 500,
#                                    outputDirectory =dir_name, projectName = "ORA_KEGG")


#   l_MF <- nrow(enrichResult_MF)
#   l_BP <- nrow(enrichResult_BP)
#   l_CC <- nrow(enrichResult_CC)
#   l_KEGG <- nrow(enrichResult_KEGG)

#   if (is.null(l_MF)){
#     l_MF <- 0
#   }
#   if (is.null(l_BP)){
#     l_BP <- 0
#   }
#   if (is.null(l_CC)){
#     l_CC <- 0
#   }
#   if (is.null(l_KEGG)){
#     l_KEGG <- 0
#   }

#   l_tot <- l_MF+l_BP+l_CC+l_KEGG
#   print(l_MF)
#   print(l_BP)
#   print(l_CC)
#   print(l_KEGG)
#   print(l_tot)

#   if (l_tot!=0){
#     col2 <- c(rep("MF",l_MF),rep("BP",l_BP),rep("CC",l_CC),rep("KEGG",l_KEGG))

#     enrichResult_all <- cbind(col2,rbind(enrichResult_MF,enrichResult_BP,enrichResult_CC,enrichResult_KEGG))
#     write.table(enrichResult_all,file=sprintf("%s/%s/%sclusters/cluster%s/GSA_cluster%s.txt",args$outdir,geneset_oi,n_clusters,cluster,cluster),
#                   col.names=TRUE,row.names=FALSE, quote=FALSE,sep="\t")
#   }

#   # TF enrichment
#   ############

#   tf_ORA <- GRN_coef.dt$tf[GRN_coef.dt$gene%in%GenesORA]
#   tf_Ref <- GRN_coef.dt$tf[GRN_coef.dt$gene%in%GenesRef]

#   dat_TF_enrichment <- NULL
#   counter <- 1
#   for (i in 1:length(unique(tf_ORA))){
#     tf_i <- unique(tf_ORA)[i]
#     X_chisq <- cbind(c(sum(tf_ORA==tf_i),sum(tf_Ref==tf_i)),
#                     c(length(GenesORA)-sum(tf_ORA==tf_i),length(GenesRef)-sum(tf_Ref==tf_i)))
#     pval_i <- chisq.test(X_chisq)$p.value
#     rate_i <- c(X_chisq[,1]/(X_chisq[,1]+X_chisq[,2]))[1]/c(X_chisq[,1]/(X_chisq[,1]+X_chisq[,2]))[2]

#     dat_TF_enrichment <- rbind(dat_TF_enrichment,c(tf_i,X_chisq[,1],X_chisq[,2],rate_i,pval_i))

#     # if ((pval_i < 0.05) & (rate_i > 1)){
#     #   print(paste0(counter," ",tf_i))
#     #   print(rate_i)
#     #   print(pval_i)
#     #   counter <- counter + 1
#     # }
#   }
#   dat_TF_enrichment <- data.frame(dat_TF_enrichment,stringsAsFactors=FALSE)
#   for (i in 2:ncol(dat_TF_enrichment)){dat_TF_enrichment[,i]<-as.numeric(dat_TF_enrichment[,i])}
#   colnames(dat_TF_enrichment) <- c("TF","n_goi","n_bg","tot_goi","tot_bg","rate","pval")
#   dat_TF_enrichment$p_adj <- dat_TF_enrichment$pval*nrow(dat_TF_enrichment)
#   dat_TF_enrichment$p_adj[dat_TF_enrichment$p_adj>1] <- 1
#   print(dat_TF_enrichment[(dat_TF_enrichment$p_adj<0.05) & (dat_TF_enrichment$rate>1),])
#   write.table(dat_TF_enrichment,file=sprintf("%s/%s/%sclusters/cluster%s/TF_enrichment_cluster%s.txt",args$outdir,geneset_oi,n_clusters,cluster,cluster),
#               col.names=TRUE,row.names=FALSE, quote=FALSE,sep="\t")
# }

# # CpG content
# ############

# CpG_content <- cbind(gsub("_exprs$","",rownames(dat_21clust)),dat_21clust$cluster)
# CpG_content <- data.frame(CpG_content,stringsAsFactors=F)
# colnames(CpG_content) <- c("gene","cluster")
# CpG_content$cluster <- as.factor(as.numeric(CpG_content$cluster))
# saveRDS(CpG_content,file="CpG_content.rds")

# library(biomaRt)
# # Import gene structure
# if (!file.exists(paste0(args$outdir,"/Genes.rds"))){
#   ensembl <- useEnsembl(biomart="ensembl", version=107, dataset="hsapiens_gene_ensembl")
#   genes_info <- getBM(attributes = c('external_gene_name','chromosome_name', 'ensembl_gene_id','start_position','end_position','transcription_start_site'),
#                       mart = ensembl,filter="external_gene_name",values=CpG_content$gene)
#   save(genes_info,file="Genes_struct.Rda")
#   print(head(genes_info))
# } else {
#   genes_info <- readRDS(paste0(args$outdir,"/Genes.rds"))
#   print(head(genes_info))
# }

# genes_info <- genes_info[!duplicated(genes_info$external_gene_name),]

# genes_info$sense <- c("AS","S")[(abs(genes_info$start_position-genes_info$transcription_start_site) < abs(genes_info$end_position-genes_info$transcription_start_site))+1]

# CpG_content$chr <- NA
# CpG_content$prom_start <- NA
# CpG_content$prom_end <- NA
# for (i in 1:nrow(CpG_content)){
#   mask_genes_info <- genes_info$external_gene_name==CpG_content$gene[i]
#   if (sum(mask_genes_info)!=0){
#     CpG_content$chr[i] <- genes_info$chromosome_name[mask_genes_info]
#     if (genes_info$sense[mask_genes_info]=="S"){
#       CpG_content$prom_start[i] <- genes_info$transcription_start_site[mask_genes_info]-2000
#       CpG_content$prom_end[i] <- genes_info$transcription_start_site[mask_genes_info]+500
#     } else {
#       CpG_content$prom_start[i] <- genes_info$transcription_start_site[mask_genes_info]-500
#       CpG_content$prom_end[i] <- genes_info$transcription_start_site[mask_genes_info]+2000
#     }
#   } else {
#     print(paste0(CpG_content$gene[i]," not in ensembl annotation!"))
#   }
# }

# HS.genome.ensembl <- getGenome( db       = "ensembl",
#                                 organism = "Homo sapiens",
#                                 release = 107,
#                                 path     = file.path("_ncbi_downloads","genomes") ,
#                                 assembly_type = "primary_assembly")


# #########
# #######

# markers_TF <- NULL
# for (i in 1:5){
#   marker_tf_i <- read.table(sprintf("/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/RNA/marker_TFs/Marker_TF_cluster%s.txt",i),header=T,sep="\t")
#   markers_TF <- rbind(markers_TF,cbind(marker_tf_i$Gene,rep(i,length(marker_tf_i$Gene))))
# }
# markers_TF <- data.frame(markers_TF)
# colnames(markers_TF) <- c("TF","celltype")
# markers_TF_undupl <- markers_TF[!duplicated(markers_TF$TF),]

# GRN_coef.dt <- fread(args$grn_coef) %>%
#   .[,gene:=toupper(gene)]

# dir.create("/data/homes/louisc/Project_Babraham/GRN_goi")
# genes_of_interest <- c("ESRRB","ARID5B","ETV1")
# for (goi in genes_of_interest){
#   GRN_coef.dt_goi <- GRN_coef.dt[GRN_coef.dt$tf==goi,]
#   GRN_coef.dt_goi <- GRN_coef.dt_goi[sort(abs(GRN_coef.dt_goi$beta),decreasing=T,index.return=T)$ix,]
#   head(GRN_coef.dt_goi)
#   write.table(GRN_coef.dt_goi,file=paste0("/data/homes/louisc/Project_Babraham/GRN_goi/GRN_coef_",goi,".txt"),col.names=T,row.names=T,sep="\t",quote=F)
#   GRN_coef.dt_goi_mTF <- GRN_coef.dt_goi[GRN_coef.dt_goi$gene%in%markers_TF_undupl$TF]
#   head(GRN_coef.dt_goi_mTF)
#   write.table(GRN_coef.dt_goi_mTF,file=paste0("/data/homes/louisc/Project_Babraham/GRN_goi/GRN_coef_",goi,"_mTF.txt"),col.names=T,row.names=T,sep="\t",quote=F)
# }

######################
## Completion token ##
######################

args$incl_samples <- "nodiff"

print("Completion token")
print(sprintf("%s/pseudotime_plots_completed_%s_%s.txt",args$outdir,args$trajectory_name,args$incl_samples))

file.create(sprintf("%s/pseudotime_plots_completed_%s_%s.txt",args$outdir,args$trajectory_name,args$incl_samples))