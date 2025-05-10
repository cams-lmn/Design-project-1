####################
##                ##
##  Clustering.R  ##
##                ##
####################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',             type="character",                               help='SingleCellExperiment file')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--samples',         type="character",  default="all",              nargs='+',     help='Samples')
p$add_argument('--features',        type="integer",    default=1000,                help='Number of features')
p$add_argument('--npcs',            type="integer",    default=30,                  help='Number of PCs')
p$add_argument('--n_neighbors',     type="integer",    default=30,     help='(UMAP) Number of neighbours')
p$add_argument('--min_dist',        type="double",     default=0.3,     help='(UMAP) Minimum distance')
p$add_argument('--n_neighbors_clusters',     type="integer",    default=15,     help='Louvain clustering number of neighbours')
p$add_argument('--resolution',        type="double",     default=0.25,     help='Louvain clustering resolution')
p$add_argument('--colour_by',       type="character",  default="celltype",  nargs='+',  help='Metadata columns to colour the UMAP by')
p$add_argument('--seed1',            type="integer",    default=42,                  help='Random seed')
p$add_argument('--seed2',            type="integer",    default=4242,                  help='Random seed')
p$add_argument('--seed3',            type="integer",    default=424242,                  help='Random seed')
p$add_argument('--seed4',            type="integer",    default=42424242,                  help='Random seed')
p$add_argument('--outdir',          type="character",                               help='Output directory')
p$add_argument('--clusterdir',          type="character",                               help='Cluster directory')
p$add_argument('--vars_to_regress', type="character",                nargs='+',     help='Metadata columns to regress out')
p$add_argument('--batch_variable',type="character",   default="None",                            help='Metadata column to apply batch correction on')
p$add_argument("--incl_samples",  type="character",    help='Which samples should be included')
p$add_argument("--sort_samples",  type="logical",default=T, help="Should samples be sorted?")
p$add_argument("--find_markers",  type="logical",default=T, help="Find markers?")
p$add_argument('--n_plots_per_cluster',            type="integer",    default=5,                  help='Amount of violin plots per cluster')
p$add_argument("--reorder_clusters",  type="integer",   default= NULL,  nargs='+', help="New numerical order for clusters")
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

dir.create(args$outdir, showWarnings = F, recursive = T)
dir.create(args$clusterdir, showWarnings = F, recursive = T)

args$filter_differentiated <- FALSE

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE & sample%in%args$samples]

if (args$sort_samples){
  table(sample_metadata$sample)[args$samples]
} else {
  table(sample_metadata$sample)
}

#####################
## Parse arguments ##
#####################

if (args$batch_variable=="None") {
  args$batch_variable <- NULL
}

###################
## Sanity checks ##
###################

stopifnot(args$colour_by %in% colnames(sample_metadata))

if (length(args$vars_to_regress)>0) {
  stopifnot(args$vars_to_regress%in%colnames(sample_metadata))
}

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(args$sce, cells=sample_metadata$cell, normalise = TRUE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

#####################
## Cell annotation ##
#####################

load("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA/scripts/clustering/mouse.rnaseq.rda")

ref <- mouse.rnaseq
pred <- SingleR(test = sce, ref = ref$data, labels = ref$main_types)

sample_metadata$annotation <- pred$labels

############################
## Regress out covariates ##
############################

if (length(args$vars_to_regress)>0) {
  print(sprintf("Regressing out variables: %s", paste(args$vars_to_regress,collapse=" ")))
  logcounts(sce_filt) <- RegressOutMatrix(
    mtx = logcounts(sce_filt),
    covariates = colData(sce_filt)[,args$vars_to_regress,drop=F]
  )
}

#######################
## Feature selection ##
#######################

decomp_full <- modelGeneVar(sce)
decomp <- decomp_full[decomp_full$mean > 0.01,]
hvgs <- decomp[order(decomp$FDR),] %>% head(n=args$features) %>% rownames

# Subset SingleCellExperiment
sce_filt <- sce[hvgs,]

###################
## PCA: sce full ##
###################

set.seed(args$seed1)

sce <- runPCA(sce, ncomponents = args$npcs, ntop=args$features)

################
## Clustering ##
################

set.seed(args$seed4)
snn.gr <- buildSNNGraph(sce, use.dimred="PCA", k=args$n_neighbors_clusters)
clusters <- factor(igraph::cluster_louvain(snn.gr,resolution=args$resolution)$membership)

sample_metadata$cluster <- clusters

if (!is.null(args$reorder_cluster)){
  print("Reordering clusters...")
  print(length(args$reorder_clusters))
  print(args$reorder_clusters)
  print(length(unique(sample_metadata$cluster)))
  print(unique(sample_metadata$cluster))
  if (length(args$reorder_clusters)!=length(unique(sample_metadata$cluster))){
    warning("Vector for cluster reordering of wrong length")
  } else {
    print(table(sample_metadata$cluster,
                args$reorder_clusters[as.numeric(sample_metadata$cluster)]))
    sample_metadata$cluster <- args$reorder_clusters[as.numeric(sample_metadata$cluster)]
  }
}

sample_metadata$cluster <- as.factor(sample_metadata$cluster)

colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame


###################
## PCA: sce filt ##
###################

set.seed(args$seed1)

sce_filt <- runPCA(sce_filt, ncomponents = args$npcs, ntop=args$features)

##########
## UMAP ##
##########

# Run
set.seed(args$seed2)
sce_filt <- runUMAP(sce_filt, dimred="PCA", n_neighbors = args$n_neighbors, min_dist = args$min_dist)

# Fetch UMAP coordinates
umap.dt <- reducedDim(sce_filt,"UMAP") %>% as.data.table %>% 
  .[,cell:=colnames(sce_filt)] %>%
  setnames(c("UMAP1","UMAP2","cell"))

##########
## plot ##
##########

pt.size <- ifelse(ncol(sce)>=1e4,0.8,1.2)

to.plot <- reducedDim(sce_filt,"UMAP") %>% as.data.table %>% 
  .[,cell:=colnames(sce_filt)] %>%
  merge(sample_metadata, by="cell")

p <- ggplot(to.plot, aes_string(x="V1", y="V2", fill="cluster")) +
  geom_point(size=pt.size, shape=21, stroke=0.05) +
  theme_classic() +
  ggplot_theme_NoAxes()

# Save UMAP plot
outfile <- file.path(args$outdir,sprintf("umap_features%d_pcs%d_neigh%d_neighclust%d_dist%s_res%s_cluster_%s.pdf",
                                           args$features, args$npcs, args$n_neighbors, args$n_neighbors_cluster, 
                                         args$min_dist, args$resolution, args$incl_samples))

pdf(outfile, width=7, height=5)
print(p)
dev.off()

p <- ggplot(to.plot, aes_string(x="V1", y="V2", fill="annotation")) +
  geom_point(size=pt.size, shape=21, stroke=0.05) +
  theme_classic() +
  ggplot_theme_NoAxes()

# Save annotation plot
outfile <- file.path(args$clusterdir,sprintf("annotated_clusters.pdf"))

pdf(outfile, width=7, height=5)
print(p)
dev.off()

##################
## Find markers ##
##################

if (args$find_markers){
  # create output dir
  dir.create(sprintf("%s/cluster_neighclust%s_res%s_%s", args$clusterdir, args$n_neighbors_clusters, 
                     args$resolution, args$incl_samples), showWarnings = F)
  fwrite(sample_metadata, sprintf("%s/cluster_neighclust%s_res%s_%s/sample_metadata_after_clustering.txt.gz",
                                  args$clusterdir,args$n_neighbors_clusters,args$resolution,args$incl_samples),
           sep="\t", na="NA", quote=F)

  
  # recreate clustering plot with labels
  label.df <- data.frame(cluster=levels(to.plot$cluster),
                         label=paste0("cluster: ",levels(to.plot$cluster)))
  label.df_2 <- to.plot %>% 
    group_by(cluster) %>% 
    summarize(V1 = quantile(V1,0.95), V2 = quantile(V2,0.95)) %>% 
    left_join(label.df)
  
  p <- ggplot(to.plot, aes_string(x="V1", y="V2", fill="cluster")) +
    geom_point(size=pt.size, shape=21, stroke=0.05) +
    theme_classic() +
    ggplot_theme_NoAxes() +
    ggrepel::geom_label_repel(data = label.df_2, aes(label = label))
  
  outfile <- file.path(sprintf("%s/cluster_neighclust%s_res%s_%s/umap_features%d_pcs%d_neigh%d_neighclust%d_dist%s_res%s_cluster.pdf",
                                 args$clusterdir, args$n_neighbors_clusters, args$resolution, args$incl_samples ,args$features, args$npcs, 
                                 args$n_neighbors, args$n_neighbors_cluster, args$min_dist, args$resolution))
  
  pdf(outfile, width=7, height=5)
  print(p)
  dev.off()
  
  # find markers 
  markers <- findMarkers(sce,groups=clusters,assay.type="logcounts",add.summary=T)
  
  # write output per cluster
  for (i in 1:length(markers)){
    print(paste0("cluster ",i," started"))
    
    outdir_marker <- sprintf("%s/cluster_neighclust%s_res%s_%s/cluster%s",
                             args$clusterdir,args$n_neighbors_clusters,args$resolution,args$incl_samples,i)
    
    
    dir.create(outdir_marker, showWarnings = F)
    
    write.table(markers[[i]],sprintf("%s/markers_cluster%s.txt",outdir_marker,i),col.names = T, row.names = T, sep="\t",quote = F)
    
    gene_list <- rownames(markers[[i]])
    gene_list <- gene_list[markers[[i]]$summary.logFC>0]
    gene_list <- gene_list[1:args$n_plots_per_cluster]
    p <- plotExpression(sce,gene_list,
                        x = "cluster",exprs_values = "logcounts")
    
    outfile_jpg <- file.path(sprintf("%s/violin_cluster%s_neighclust%s_res%s.jpg",
                                       outdir_marker, i, args$n_neighbors_clusters, args$resolution))
    outfile_pdf <- file.path(sprintf("%s/violin_cluster%s_neighclust%s_res%s.pdf",
                                       outdir_marker, i, args$n_neighbors_clusters, args$resolution))

    jpeg(outfile_jpg,width = 700, height = 500, units = "px")
    print(p)
    dev.off()
    
    pdf(outfile_pdf, width=7, height=5)
    print(p)
    dev.off()
  }
}