####################
##                ##
##  clustering.R  ##
##                ##
####################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/ATAC/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',    type="character",    help='ArchR directory')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--matrix',          type="character",  default="PeakMatrix",   help='Matrix to use')
p$add_argument('--samples',         type="character",                nargs='+',     help='Samples')
p$add_argument('--genome',          type="character", default="hg38",      help='Genome')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
p$add_argument('--nfeatures',       type="integer",    default=1000,               help='Number of features')
p$add_argument('--ndims',           type="integer",    default=30,                  help='Number of LSI dimensions')
p$add_argument('--batch_variable',  type="character",                               help='Metadata column to apply batch correction on')
p$add_argument('--batch_method',    type="character",  default="MNN",               help='Batch correctin method ("Harmony" or "MNN")')
p$add_argument('--n_neighbors',     type="integer",    default=30,   nargs='+',     help='(UMAP) Number of neighbours')
p$add_argument('--min_dist',        type="double",     default=0.3,  nargs='+',     help='(UMAP) Minimum distance')
p$add_argument('--res',        type="double",     default=0.3,  nargs='+',     help='Resolution for louvain clustering')
p$add_argument('--seed',            type="integer",    default=42,                  help='Random seed')
p$add_argument('--seed3',            type="integer",    default=42,                  help='Random seed 3')
p$add_argument('--outdir',          type="character",                               help='Output directory')
p$add_argument('--clusterdir',          type="character",                               help='Cluster directory')
p$add_argument("--sort_samples",  type="logical",default=TRUE, help="Should samples be sorted?")
p$add_argument("--incl_samples",  type="character",    help='Which samples should be included')
p$add_argument("--reorder_clusters",  type="integer", nargs="+", default=NULL   , help='Should clusters be reordered')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# args$metadata <- "/data/homes/louisc/Project_Babraham/ATAC/archR/qc/sample_metadata_after_qc.txt.gz"
# args$nfeatures <- 25000
# args$matrix <- "PeakMatrix"
# args$incl_samples <- "nodiff"
# args$sort_samples <- TRUE
# args$ndims <- 15
# args$seed <- 42
# args$n_neighbors <- 25
# args$min_dist <- 0.25
# args$genome <- "hg38"
# args$batch_variable <- "SeqRun"
# args$batch_method <- "Harmony"
# args$threads <- 1
# args$seed3 <- 42
# args$res <- 0.15
# args$samples <- c("d0","d1","d3","d5","d7","d10","d14","d18","DE","NE","AME_E","AME_L")
# args$outdir <- sprintf("/data/homes/louisc/Project_Babraham/ATAC/archR/dimensionality_reduction/%s/batch_correction_%s_%s",args$matrix,args$batch_variable,args$batch_method)
# args$archr_directory <- "/data/homes/louisc/Project_Babraham/ATAC/archR"
# args$clusterdir <- sprintf("/data/homes/louisc/Project_Babraham/ATAC/archR/clustering/%s/batch_correction_%s_%s",args$matrix,args$batch_variable,args$batch_method)
# args$reorder_clusters <- c(1,2,3,4,5,6)
# args$test <- "Test mode still on!"
## END TEST ##

print(args)

#####################
## Define settings ##
#####################

# Options
opts$lsi.iterations = 2
opts$lsi.cluster.resolution = 2

if(grepl("nodiff",args$incl_samples)){
  args$filter_differentiated <- TRUE
} else {
  args$filter_differentiated <- FALSE
}

dir.create(args$clusterdir, showWarnings = F)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE] %>%
  .[,log_nFrags_atac:=log10(nFrags_atac)]

#############################
## Filter naieve to primed ##
#############################

if (args$filter_differentiated){
  print("Removing differentiated cells...")
  sample_metadata <- sample_metadata[grepl("d[0-9]+",sample_metadata$sample),]
  args$samples <- args$samples[grepl("d[0-9]+",args$samples)]
  if (args$sort_samples){
    table(sample_metadata$sample)[args$samples]
  } else {
    table(sample_metadata$sample)
  }
}     

#################
## Filter plot ##
#################

if(!is.null(args$max_point_per_sample)){
  sample_metadata_new <- NULL
  for (j in 1:length(unique(sample_metadata$sample))){
    set.seed(args$seed3*j)
    if (sum(sample_metadata$sample==unique(sample_metadata$sample)[j])>args$max_point_per_sample){
      sample_metadata_new <- rbind(sample_metadata_new,sample_metadata[sample_metadata$sample==unique(sample_metadata$sample)[j],][sample(1:sum(sample_metadata$sample==unique(sample_metadata$sample)[j]),args$max_point_per_sample,replace=F),])
    } else {
      sample_metadata_new <- rbind(sample_metadata_new,sample_metadata[sample_metadata$sample==unique(sample_metadata$sample)[j],])
    }
  }
  sample_metadata <- sample_metadata_new
  if (args$sort_samples){
    table(sample_metadata$sample)[args$samples]
  } else {
    table(sample_metadata$sample)
  }
}


########################
## Load ArchR project ##
########################

# source(here::here("atac/archR/load_archR_project.R"))

setwd(args$archr_directory)

addArchRGenome(args$genome)
addArchRThreads(threads = args$threads)

ArchRProject <- loadArchRProject(args$archr_directory)
ArchRProject.filt <- ArchRProject[sample_metadata$cell]
ArchRProject.filt

###################
## Sanity checks ##
###################

# stopifnot(args$matrix %in% getAvailableMatrices(ArchRProject))
# 
# if (length(args$batch_variable)>0) {
#   stopifnot(args$batch_variable%in%colnames(sample_metadata))
#   if (length(unique(sample_metadata[[args$batch_variable]]))==1) {
#     message(sprintf("There is a single level for %s, no batch correction applied",args$batch_variable))
#     args$batch_variable <- NULL
#   } else {
#     library(batchelor)
#   }
# }
# 
# if (length(args$vars_to_regress)>0) {
#   stopifnot(args$vars_to_regress%in%colnames(sample_metadata))
# }

#############################
## Add sequencing Run info ##
#############################

if (args$batch_variable=="SeqRun"){
  sample_metadata$SeqRun <- "Run2"
  sample_metadata$SeqRun[sample_metadata$sample=="d0"] <- "Jasmin"
  sample_metadata$SeqRun[sample_metadata$sample%in%c("d10","d14","d18","DE")] <- "Run1"
} else if (args$batch_variable=="None"){
  args$batch_variable <- NULL
}


###########################
## Update ArchR metadata ##
###########################

print(head(sample_metadata))
print(head(rownames(getCellColData(ArchRProject.filt))))

sample_metadata.to.archr <- sample_metadata %>%
  .[cell%in%rownames(getCellColData(ArchRProject.filt))] %>% 
  as.data.frame() %>% tibble::column_to_rownames("cell")

# stopifnot(all(rownames(sample_metadata.to.archr) == rownames(getCellColData(ArchRProject.filt))))
for (i in args$batch_variable) {
  ArchRProject.filt <- addCellColData(
    ArchRProject.filt,
    data = sample_metadata.to.archr[[i]],
    name = i,
    cells = rownames(sample_metadata.to.archr),
    force = TRUE
  )
  print(table(getCellColData(ArchRProject.filt,i)[[1]]))
}

#######################
## Feature selection ##
#######################

if (args$matrix=="PeakMatrix_filt") {
  
  # Load peak variability estimates
  peak.variability.dt <- fread(io$archR.peak.variability)
  
  # Define highly variable peaks
  peaks <- peak.variability.dt %>%
    .[,peak:=sub("_",":",peak)] %>% .[,peak:=sub("_","-",peak)] %>%
    setorder(-variance_pseudobulk) %>%
    head(n=args$nfeatures) %>%
    .$peak
  
  # Subset peaks in the ArchR object
  names(ArchRProject.filt@peakSet) <- sprintf("%s:%s-%s",seqnames(ArchRProject.filt@peakSet), start(ArchRProject.filt@peakSet), end(ArchRProject.filt@peakSet))
  ArchRProject.filt@peakSet[peaks,] <- ArchRProject.filt@peakSet[peaks,]
  
  ArchRProject.filt <- addFeatureMatrix(
    input = ArchRProject.filt,
    features = ArchRProject.filt@peakSet,
    matrixName = "PeakMatrix_filt",
    binarize = TRUE
  )
  
}

###########################
## Latent Semantic Index ##
###########################

# Iterative LSI: two iterations
ArchRProject.filt <- addIterativeLSI(
  ArchRProj = ArchRProject.filt,
  useMatrix = args$matrix, 
  name = "IterativeLSI", 
  firstSelection = "Top",
  depthCol = "nFrags",
  iterations = opts$lsi.iterations, 
  # clusterParams = list(
  #  resolution = opts$lsi.cluster.resolution, 
  #  sampleCells = 10000, 
  #  n.start = 10
  # ), 
  dimsToUse=1:args$ndims, #uncomment to recreate ATAC UMAPs and ATAC DE analysis
  saveIterations = FALSE,
  varFeatures = args$nfeatures, 
  force = TRUE
)

############################
## LSI + Batch correction ##
############################

if (length(args$batch_variable)>0) {
  print(sprintf("Applying %s batch correction for variable: %s", args$batch_method, args$batch_variable))
  outfile <- sprintf("%s/lsi_nfeatures%d_ndims%d_%s.txt.gz",
                     args$clusterdir, args$nfeatures, args$ndims, args$incl_samples)
  
  # Harmony
  if (args$batch_method=="Harmony") {
    ArchRProject.filt <- addHarmony(
      ArchRProj = ArchRProject.filt,
      reducedDims = "IterativeLSI",
      name = "IterativeLSI_Harmony",
      groupBy = args$batch_variable,
      force = TRUE
    )
    lsi.dt <- getReducedDims(ArchRProject.filt, "IterativeLSI_Harmony") %>% round(3) %>% 
      as.data.table(keep.rownames = T) %>% setnames("rn","cell")
    
  } else {
    stop("Batch correction method not recognised")
  }
} else {
  lsi.dt <- getReducedDims(ArchRProject.filt, "IterativeLSI") %>% round(3) %>% 
    as.data.table(keep.rownames = T) %>% setnames("rn","cell")
  outfile <- sprintf("%s/lsi_nfeatures%d_ndims%d_%s.txt.gz",
                     args$clusterdir, args$nfeatures, args$ndims, args$incl_samples)
}

# Save LSI coordinates
fwrite(lsi.dt, outfile)
print(head(lsi.dt))

################
## Clustering ##
################

ArchRProject.filt <- addClusters(
  input = ArchRProject.filt,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = args$res,
  force=T
)

if (!is.null(args$reorder_cluster)){
  print("Reordering clusters...")
  print(length(args$reorder_clusters))
  print(args$reorder_clusters)
  print(length(unique(getCellColData(ArchRProject.filt)$Clusters)))
  print(unique(getCellColData(ArchRProject.filt)$Clusters))
  if (length(args$reorder_clusters)!=length(unique(getCellColData(ArchRProject.filt)$Clusters))){
    sample_metadata$Clusters <- getCellColData(ArchRProject.filt)$Clusters
    warning("Vector for cluster reordering of wrong length")
  } else {
    print(table(getCellColData(ArchRProject.filt)$Clusters,
                args$reorder_clusters[as.numeric(gsub("^C","",getCellColData(ArchRProject.filt)$Clusters))]))
    sample_metadata$Clusters <- paste0(rep("C",length(getCellColData(ArchRProject.filt)$Clusters)),
                                       args$reorder_clusters[as.numeric(gsub("^C","",getCellColData(ArchRProject.filt)$Clusters))])
  }
} else {
  sample_metadata$Clusters <- getCellColData(ArchRProject.filt)$Clusters
}

###########################
## Update ArchR metadata ##
###########################

sample_metadata.to.archr <- sample_metadata %>%
  .[cell%in%rownames(getCellColData(ArchRProject.filt))] %>% 
  as.data.frame() %>% tibble::column_to_rownames("cell")


ArchRProject.filt <- addCellColData(
  ArchRProject.filt,
  data = sample_metadata.to.archr$Clusters,
  name = "Clusters",
  cells = rownames(sample_metadata.to.archr),
  force = TRUE
)
print(table(getCellColData(ArchRProject.filt,"Clusters")[[1]]))


##########
## UMAP ##
##########

if (is.null(args$max_point_per_sample)){
  pt.size <- 0.8
} else {
  pt.size <- 1.2
}

i <- args$n_neighbors
j <- args$min_dist

# Run UMAP
ArchRProject.filt <- addUMAP(
  ArchRProj = ArchRProject.filt, 
  reducedDims = "IterativeLSI",
  name = "UMAP",
  metric = "cosine",
  nNeighbors = i, 
  minDist = j, 
  seed = args$seed,
  saveModel = FALSE,
  force = TRUE
)

# Fetch UMAP coordinates
umap.dt <- getEmbedding(ArchRProject.filt,"UMAP") %>%
  round(2) %>%
  as.data.table(keep.rownames = T) %>%
  setnames(c("cell","umap1","umap2"))

# Save UMAP coordinates
# outfile <- sprintf("%s/umap_%s_nfeatures%d_ndims%d_neigh%d_dist%s.txt.gz",args$outdir, args$matrix, args$nfeatures, args$ndims, i, j)
fwrite(umap.dt,  sprintf("%s/umap_nfeatures%d_ndims%d_%s.txt.gz",args$clusterdir, args$nfeatures, args$ndims, args$incl_samples))

# Plot
to.plot <- umap.dt %>%
  merge(sample_metadata,by="cell")

k <- "Clusters"

# log10 large numeric values
if (is.numeric(to.plot[[k]])) {
  if (max(to.plot[[k]],na.rm=T) - min(to.plot[[k]],na.rm=T) > 1000) {
    to.plot[[k]] <- log10(to.plot[[k]]+1)
    to.plot %>% setnames(k,paste0(k,"_log10")); k <- paste0(k,"_log10")
  }
}

if (args$sort_samples){
  to.plot$sample <- as.factor(to.plot$sample)
  to.plot$sample <- factor(to.plot$sample,levels=args$samples)
}

p <- ggplot(to.plot, aes_string(x="umap1", y="umap2",fill=k)) +
  geom_point(size=pt.size, shape=21, stroke=0.05) +
  # ggrastr::geom_point_rast(size=1.5, shape=21, stroke=0.05) +  # DOES NOT WORK IN THE CLUSTER
  theme_classic() +
  ggplot_theme_NoAxes()

# Save UMAP plot
output_plot(p, sprintf("%s/umap_%s_nfeatures%d_ndims%d_%s_res%s_%s",args$outdir, args$matrix, args$nfeatures, args$ndims, k, 
                       args$res, args$incl_samples), width=7, height=5, UMAP=TRUE)

###################
## Plot clusters ##
###################

to.plot$cluster <- gsub("^C","",to.plot$Clusters)
to.plot$cluster <- as.factor(to.plot$cluster)

# create output dir
dir.create(sprintf("%s/cluster_res%s_%s",args$clusterdir,args$res,args$incl_samples), showWarnings = F)
fwrite(sample_metadata, sprintf("%s/cluster_res%s_%s/sample_metadata_after_clustering.txt.gz",
                                args$clusterdir,args$res,args$incl_samples),
       sep="\t", na="NA", quote=F)

p <- ggplot(to.plot, aes_string(x="umap1", y="umap2", fill="cluster")) +
  geom_point(size=pt.size, shape=21, stroke=0.05) +
  theme_classic() +
  ggplot_theme_NoAxes() +
  #ggrepel::geom_label_repel(data = label.df_2, aes(label = label)) + 
  scale_fill_manual(values=opts$celltype.colors[1:length(levels(to.plot$cluster))]) 


output_plot(p, sprintf("%s/cluster_res%s_%s/umap_%s_nfeatures%d_ndims%s_%s_res%s", args$clusterdir, args$res, 
                         args$incl_samples, args$matrix, args$nfeatures, args$ndims, k, args$res), width=7, height=5, UMAP=TRUE)

###################
## Save metadata ##
###################

# Save
fwrite(sample_metadata, sprintf("%s/cluster_res%s_%s/sample_metadata_after_clustering.txt.gz",
                                  args$clusterdir, args$res, args$incl_samples), quote=F, na="NA", sep="\t")
                              
print(sessionInfo())