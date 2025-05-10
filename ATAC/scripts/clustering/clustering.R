####################
##                ##
##  clustering.R  ##
##                ##
####################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',    type="character",    help='ArchR directory')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--matrix',          type="character",  default="PeakMatrix",   help='Matrix to use')
p$add_argument('--samples',         type="character",                nargs='+',     help='Samples')
p$add_argument('--genome',          type="character", default="mm10",      help='Genome')
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
print(args)

#####################
## Define settings ##
#####################

# Options
opts$lsi.iterations = 2
opts$lsi.cluster.resolution = 2

args$filter_differentiated <- FALSE
if (args$batch_variable=="None"){
  args$batch_variable <- NULL
}

dir.create(args$clusterdir, showWarnings = F)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE] %>%
  .[,log_nFrags_atac:=log10(nFrags_atac)]

########################
## Load ArchR project ##
########################

setwd(args$archr_directory)

addArchRGenome(args$genome)
addArchRThreads(threads = args$threads)

ArchRProject <- loadArchRProject(args$archr_directory)
ArchRProject.filt <- ArchRProject[sample_metadata$cell]
ArchRProject.filt

###########################
## Update ArchR metadata ##
###########################

print(head(sample_metadata))
print(head(rownames(getCellColData(ArchRProject.filt))))

sample_metadata.to.archr <- sample_metadata %>%
  .[cell%in%rownames(getCellColData(ArchRProject.filt))] %>% 
  as.data.frame() %>% tibble::column_to_rownames("cell")

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
  dimsToUse=1:args$ndims, #uncomment to recreate ATAC UMAPs and ATAC DE analysis
  saveIterations = FALSE,
  varFeatures = args$nfeatures, 
  force = TRUE
)

#########
## LSI ##
#########

lsi.dt <- getReducedDims(ArchRProject.filt, "IterativeLSI") %>% round(3) %>% 
  as.data.table(keep.rownames = T) %>% setnames("rn","cell")
outfile <- sprintf("%s/lsi_nfeatures%d_ndims%d_%s.txt.gz",
                    args$clusterdir, args$nfeatures, args$ndims, args$incl_samples)

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

# Get current cluster assignments
sample_metadata$Clusters <- getCellColData(ArchRProject.filt)$Clusters

# Define cluster mapping (merging clusters)
# (C1 and C3) -> C1
# (C2 and C4) -> C3
# C5 -> C2
# (C6, C10 and C11) -> C4
# (C8 and C9) -> C6
# C7 -> C5
cluster_mapping <- list(
  "C1" = "C1", "C3" = "C1",
  "C2" = "C3", "C4" = "C3",
  "C5" = "C2",
  "C6" = "C4", "C10" = "C4", "C11" = "C4",
  "C8" = "C6", "C9" = "C6",
  "C7" = "C5"
)

sample_metadata$Clusters <- as.character(sapply(sample_metadata$Clusters, function(x) cluster_mapping[[x]]))

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

pt.size <- 0.8

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