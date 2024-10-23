#####################
##                 ##
##  Markers_TFs.R  ##
##                 ##
#####################

source("/data/homes/louisc/Project_Babraham/RNA_ATAC/scripts/Settings.R")
library("viridis")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--DEG_overview',    type="character",    help='Differential expression analysis overview')
p$add_argument('--sce',       type="character",                help='SingleCellExperiment')
p$add_argument('--sce_mc',       type="character",                help='SingleCellExperiment on metacell level')
p$add_argument('--TF',    type="character",    help='Transcription factor file')
p$add_argument('--seed',            type="integer",     default=42,             help='Random seed')
p$add_argument('--ncores',            type="integer",     default=1,             help='Number of cores')
p$add_argument('--matrix',          type="character",  nargs="+",    help='Matrix to use')
p$add_argument('--mofa_metadata',    type="character",    help='MOFA metadata')
p$add_argument('--factors',    type="character",    help='MOFA factors')
p$add_argument('--outdir',   type="character",    help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

# START TEST
args <- list()
args$sce <- "/data/homes/louisc/Project_Babraham/RNA//SingleCellExperiment.rds"
args$sce_mc <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/metacells/SingleCellExperiment_metacells_nodiff.rds"
args$DEG_overview <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/RNA/DE_res/DEG_overview.txt"
args$TF <- "/data/homes/louisc/Project_Babraham/TFs.txt"
args$incl_samples <- "nodiff"
args$mofa_metadata <-  sprintf("/data/homes/louisc/Project_Babraham/RNA_ATAC/mofa/sample_metadata_%s.txt.gz",args$incl_samples)
args$factors <- sprintf("/data/homes/louisc/Project_Babraham/RNA_ATAC/mofa/pdf/factors_%s.txt.gz",args$incl_samples)
args$ncores <- 10
args$seed <-  42
args$outdir <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/RNA/marker_TFs"
# END TEST

print(args)

#####################
## Define settings ##
#####################

print(args$outdir)

unlink(args$outdir,recursive = T,force = T)
dir.create(args$outdir)



###############
## Load data ##
###############

DEG_overview_RNA <- fread(args$DEG_overview)

TFs <- fread(args$TF,header=F)$V1

sce <- readRDS(args$sce)
sce_mc <- readRDS(args$sce_mc)

# DEG_overview_GeneScore_distal <- fread("/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/GeneScoreMatrix_distal/DE_res/DEG_overview.txt")

########################
## Filter for real TF ##
########################

TF_annotation <- fread("/data/homes/louisc/Project_Babraham/Homo_sapiens_TF.txt")[-1,1:3]
print(head(TF_annotation))

print("Filtering for 'fake TFs'")
print(length(TFs))
TFs <- TFs[TFs%in%TF_annotation$Symbol]
print(length(TFs))

###########################
## Marker TF per cluster ##
###########################

DEG_overview_TF <- data.frame(DEG_overview_RNA[DEG_overview_RNA$Gene%in%TFs])
dim(DEG_overview_TF)

marker_TF <- NULL
for (i in 1:5){
  print(paste0("cluster ",i))
  DE_TF_cluster_i <- DEG_overview_TF[(DEG_overview_TF[,colnames(DEG_overview_TF)==paste0("clust",i,"_spec")]) &
                                       !is.na(DEG_overview_TF[,colnames(DEG_overview_TF)==paste0("clust",i,"_spec")]),]
  print(dim(DE_TF_cluster_i))
  print(DE_TF_cluster_i$Gene)
  marker_TF <- rbind(marker_TF,cbind(DE_TF_cluster_i$Gene,rep(i,nrow(DE_TF_cluster_i))))
  # print(head(DE_TF_cluster_i))
  
  write.table(DE_TF_cluster_i,file=sprintf("%s/Marker_TF_cluster%s.txt",args$outdir,i),col.names=T,row.names = F,quote = F,sep="\t")
}

##########
## Umap ##
##########

factors.dt <- fread(args$factors)
factors.dt <- factors.dt %>% column_to_rownames("cell")

# Run
set.seed(args$seed)
umap_embedding <- uwot::umap(factors.dt, n_neighbors=25, min_dist=0.50, metric="cosine")

mofa_metadata <- fread(args$mofa_metadata)

#################
## Plots Maria ##
#################

theme_set(theme_bw(base_size=14))

### Colors options

cols1 <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
cols2 <- colorRampPalette(c("#DBD8D7","#F9F9F9","#FCBBA2","#FB8B6A","#F86546","#E93529","#9E050D"))(100)
cols3 <- colorRampPalette(rev(rocket(8)))(100)

cols_scale <- colorRampPalette(inferno(8))(100)

# Read tables 
########################

read_table(paste0(args$outdir,"/Marker_TF_cluster1.txt")) %>%
  as_tibble() %>%
  dplyr::select(1:6) -> TF.cluster.1
print(head(TF.cluster.1))

read_table(paste0(args$outdir,"/Marker_TF_cluster2.txt")) %>%
  as_tibble() %>%
  dplyr::select(1:6) -> TF.cluster.2

read_table(paste0(args$outdir,"/Marker_TF_cluster3.txt")) %>%
  as_tibble() %>%
  dplyr::select(1:6) -> TF.cluster.3

read_table(paste0(args$outdir,"/Marker_TF_cluster4.txt")) %>%
  as_tibble() %>%
  dplyr::select(1:6) -> TF.cluster.4

read_table(paste0(args$outdir,"/Marker_TF_cluster5.txt")) %>%
  as_tibble() %>%
  dplyr::select(1:6) -> TF.cluster.5

TF.cluster.1 %>%
  full_join(TF.cluster.2) %>%
  full_join(TF.cluster.3) %>%
  full_join(TF.cluster.4) %>%
  full_join(TF.cluster.5) %>%
  dplyr::rename(cluster1 = log_av_1,
                cluster2 = log_av_2,
                cluster3 = log_av_3,
                cluster4 = log_av_4,
                cluster5 = log_av_5) -> TF.markers

print("TF.markers")
print(TF.markers)


# Make groups of TFs 
###############################

TF.cluster.1 %>%
  anti_join(TF.cluster.2) %>%
  anti_join(TF.cluster.3) %>%
  anti_join(TF.cluster.4) %>%
  anti_join(TF.cluster.5) %>%
  mutate(Cluster = "1") %>%
  mutate(Cluster.all = "1") %>%
  dplyr::rename(cluster1 = log_av_1,
                cluster2 = log_av_2,
                cluster3 = log_av_3,
                cluster4 = log_av_4,
                cluster5 = log_av_5) -> TF.cluster.1.only


TF.cluster.1 %>%
  inner_join(TF.cluster.2) %>%
  anti_join(TF.cluster.3) %>%
  anti_join(TF.cluster.4) %>%
  anti_join(TF.cluster.5) %>%
  mutate(Cluster = "1") %>%
  mutate(Cluster.all = "1,2") %>%
  dplyr::rename(cluster1 = log_av_1,
                cluster2 = log_av_2,
                cluster3 = log_av_3,
                cluster4 = log_av_4,
                cluster5 = log_av_5) -> TF.cluster.1.2


TF.cluster.2 %>%
  anti_join(TF.cluster.1) %>%
  anti_join(TF.cluster.3) %>%
  anti_join(TF.cluster.4) %>%
  anti_join(TF.cluster.5) %>%
  mutate(Cluster = "2") %>%
  mutate(Cluster.all = "2") %>%
  dplyr::rename(cluster1 = log_av_1,
                cluster2 = log_av_2,
                cluster3 = log_av_3,
                cluster4 = log_av_4,
                cluster5 = log_av_5) -> TF.cluster.2.only


TF.cluster.2 %>%
  inner_join(TF.cluster.3) %>%
  anti_join(TF.cluster.1) %>%
  anti_join(TF.cluster.4) %>%
  anti_join(TF.cluster.5) %>%
  mutate(Cluster = "2") %>%
  mutate(Cluster.all = "2,3") %>%
  dplyr::rename(cluster1 = log_av_1,
                cluster2 = log_av_2,
                cluster3 = log_av_3,
                cluster4 = log_av_4,
                cluster5 = log_av_5) -> TF.cluster.2.3


TF.cluster.3 %>%
  anti_join(TF.cluster.2) %>%
  anti_join(TF.cluster.1) %>%
  anti_join(TF.cluster.4) %>%
  anti_join(TF.cluster.5) %>%
  mutate(Cluster = "3") %>%
  mutate(Cluster.all = "3") %>%
  dplyr::rename(cluster1 = log_av_1,
                cluster2 = log_av_2,
                cluster3 = log_av_3,
                cluster4 = log_av_4,
                cluster5 = log_av_5) -> TF.cluster.3.only


TF.cluster.4 %>%
  anti_join(TF.cluster.2) %>%
  anti_join(TF.cluster.3) %>%
  anti_join(TF.cluster.1) %>%
  anti_join(TF.cluster.5) %>%
  mutate(Cluster = "4") %>%
  mutate(Cluster.all = "4") %>%
  dplyr::rename(cluster1 = log_av_1,
                cluster2 = log_av_2,
                cluster3 = log_av_3,
                cluster4 = log_av_4,
                cluster5 = log_av_5) -> TF.cluster.4.only


TF.cluster.4 %>%
  inner_join(TF.cluster.5) %>%
  anti_join(TF.cluster.1) %>%
  anti_join(TF.cluster.2) %>%
  anti_join(TF.cluster.3) %>%
  mutate(Cluster = "4") %>%
  mutate(Cluster.all = "4,5") %>%
  dplyr::rename(cluster1 = log_av_1,
                cluster2 = log_av_2,
                cluster3 = log_av_3,
                cluster4 = log_av_4,
                cluster5 = log_av_5) -> TF.cluster.4.5


TF.cluster.5 %>%
  anti_join(TF.cluster.1) %>%
  anti_join(TF.cluster.2) %>%
  anti_join(TF.cluster.3) %>%
  anti_join(TF.cluster.4) %>%
  mutate(Cluster = "5") %>%
  mutate(Cluster.all = "5") %>%
  dplyr::rename(cluster1 = log_av_1,
                cluster2 = log_av_2,
                cluster3 = log_av_3,
                cluster4 = log_av_4,
                cluster5 = log_av_5) -> TF.cluster.5.only

# Annotation 
##########################

ann_cols = list(cluster = c(cluster1="#daf7a6", cluster2="#ffc300", cluster3="#a68004", cluster4="#c70039", cluster5="#900c3f"))

annotation = data.frame(cluster = c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5"),
                        cluster1 = c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5"))

# Heatmaps
########################

TF.cluster.1.only %>%
  full_join(TF.cluster.1.2) %>%
  full_join(TF.cluster.2.only) %>%
  full_join(TF.cluster.2.3) %>%
  full_join(TF.cluster.3.only) %>%
  full_join(TF.cluster.4.only) %>%
  full_join(TF.cluster.4.5) %>%
  full_join(TF.cluster.5.only) -> TF.markers

TF.markers %>%
  group_by(Cluster.all) %>%
  dplyr::count()

print("TF.markers")
print(dim(TF.markers))
print(TF.markers)

TF.markers %>%
  dplyr::select(-Cluster, -Cluster.all) %>%
  column_to_rownames(var = "Gene") %>%
  pheatmap(cluster_rows = F,
           cluster_cols = F,
           show_rownames = T,
           show_colnames = F,
           scale = "row",
           fontsize_row = 8,
           annotation_col = annotation %>% column_to_rownames(var="cluster1"),
           annotation_colors = ann_cols,
           color = colorRampPalette(cols3)(100)
           ,breaks=seq(from=-1,to=1.5,length.out=101),
           gaps_row = c(11,18,28,32,36,41,63),
           filename = paste0(args$outdir,"/10X.multiome.transition.TF.markers.clusters.pdf"), width = 10, height = 10)

TF.markers %>%
  write_tsv(paste0(args$outdir,"/TF.markers.clusters.txt"))

#####################################
## Expression plots for marker TFs ##
#####################################

dir.create(sprintf("%s/expr_plots",args$outdir))
dir.create(sprintf("%s/expr_mc_plots",args$outdir))

# to.plot <- umap_embedding %>% as.data.table %>%
#   .[,sample:=rownames(factors.dt)] %>%
#   merge(mofa_metadata[,c("sample","batch")] %>% as.data.table)

to.plot <- umap_embedding[,1:4]

to.plot_mc <- to.plot[to.plot$sample%in%colnames(sce_mc),]
rownames(to.plot_mc) <- to.plot_mc$sample
to.plot_mc <- to.plot_mc[colnames(sce_mc),]

color_scheme <- list(cols1,cols2,cols3)
names(color_scheme) <- c("RdBl","GrRd","rocket")
for (j in 1:length(color_scheme)){
  print(names(color_scheme)[j])
  dir.create(sprintf("%s/expr_plots/%s",args$outdir,names(color_scheme)[j]))
  dir.create(sprintf("%s/expr_mc_plots/%s",args$outdir,names(color_scheme)[j]))

  myfun <- function(i){
    print(paste0(TF.markers$Gene[i],": marker for cluster ", TF.markers$Cluster.all[i]))

    #print("Cell level")
    logcounts_TF <- assay(sce,"logcounts")[TF.markers$Gene[i],] %>% data.frame %>%
      rownames_to_column("sample")
    colnames(logcounts_TF)[2] <- "log_expr"
    to.plot_i <- to.plot %>% merge(logcounts_TF)
    
    to.plot_i$log_expr[to.plot_i$log_expr==0] <- NA

    to.plot_i$minmax_expr <- minmax.10pct.normalisation(to.plot_i$log_expr)

    #print(head(to.plot_i))
    
    p <- ggplot(to.plot_i, aes(x=V1, y=V2)) +
      geom_point(aes(x=V1, y=V2), size=1, color="grey", alpha=0.4, data=to.plot_i[is.na(to.plot_i$minmax_expr),]) +
      geom_point(aes(x=V1, y=V2, fill=minmax_expr), size=2, stroke=0.1, shape=21, alpha=1.0, data=to.plot_i[!is.na(to.plot_i$minmax_expr),]) +
      # geom_point(size=1, shape=21, stroke=0.01) +
      # ggrastr::geom_point_rast(size=1, shape=21, stroke=0.1) +
      theme_classic() +
      scale_fill_gradientn(colours = color_scheme[[j]]) +
      # theme(legend.position="none") +
      ggplot_theme_NoAxes() 
    
    output_plot(p,sprintf("%s/expr_plots/%s/%s_marker_cluster%s",args$outdir,names(color_scheme)[j],TF.markers$Gene[i],TF.markers$Cluster.all[i]), width=7, height=5, UMAP=TRUE)

    #print("Metacell level")
    logcounts_mc_TF <- assay(sce_mc,"logcounts")[TF.markers$Gene[i],] %>% data.frame %>%
      rownames_to_column("sample")
    colnames(logcounts_mc_TF)[2] <- "log_expr"

    to.plot_mc_i <- to.plot_mc %>% merge(logcounts_mc_TF)

    to.plot_mc_i$log_expr[to.plot_mc_i$log_expr==0] <- NA

    to.plot_mc_i$minmax_expr <- minmax.10pct.normalisation(to.plot_mc_i$log_expr)

    #print(head(to.plot_mc_i))

    p <- ggplot(to.plot_i, aes(x=V1, y=V2)) +
      geom_point(aes(x=V1, y=V2), size=1, color="grey", alpha=0.4, data=to.plot_i) +
      geom_point(aes(x=V1, y=V2, fill=minmax_expr), size=4, stroke=0.1, color="black", shape=21, alpha=1.0,  data=to.plot_mc_i[is.na(to.plot_mc_i$minmax_expr),]) +
      geom_point(aes(x=V1, y=V2, fill=minmax_expr), size=4, stroke=0.1, shape=21, alpha=1.0, data=to.plot_mc_i[!is.na(to.plot_mc_i$minmax_expr),]) +
      # geom_point(size=1, shape=21, stroke=0.01) +
      # ggrastr::geom_point_rast(size=1, shape=21, stroke=0.1) +
      theme_classic() +
      scale_fill_gradientn(colours = color_scheme[[j]]) +
      # theme(legend.position="none") +
      ggplot_theme_NoAxes() 
    
    output_plot(p,sprintf("%s/expr_mc_plots/%s/%s_marker_cluster%s",args$outdir,names(color_scheme)[j],TF.markers$Gene[i],TF.markers$Cluster.all[i]), width=7, height=5, UMAP=TRUE)
  }
  res <- mclapply(1:nrow(TF.markers),myfun,mc.cores =args$ncores,mc.preschedule=T)
}

# Cluster 4vs5
mask_c4vs5_RNA <- (abs(DEG_overview_RNA$logFC_4vs5)>1) & (DEG_overview_RNA$FDR4vs5<0.01) & 
  (!is.na(DEG_overview_RNA$logFC_4vs5)) & (!is.na(DEG_overview_RNA$FDR4vs5))
mask_c4vs5_TF <- (abs(DEG_overview_TF$logFC_4vs5)>1) & (DEG_overview_TF$FDR4vs5<0.01) & 
  (!is.na(DEG_overview_TF$logFC_4vs5)) & (!is.na(DEG_overview_TF$FDR4vs5))
sort(DEG_overview_RNA$Gene[mask_c4vs5_RNA])
sort(DEG_overview_TF$Gene[mask_c4vs5_TF])
DEG_overview_TF[mask_c4vs5_TF,]

