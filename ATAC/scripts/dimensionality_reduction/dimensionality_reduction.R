########################################
##                                    ##
##  archR_dimensionality_reduction.R  ##
##                                    ##
########################################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',    type="character",    help='ArchR directory')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--matrix',          type="character",  default="PeakMatrix",   help='Matrix to use')
p$add_argument('--samples',         type="character",                nargs='+',     help='Samples')
p$add_argument('--seed',            type="integer",    default=42,                  help='Random seed')
p$add_argument('--seed3',     type="integer",    default=424242,    help='Random seed 3')
p$add_argument('--genome',          type="character", default="mm10",      help='Genome')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
p$add_argument('--nfeatures',       type="integer",    default=1000,               help='Number of features')
p$add_argument('--ndims',           type="integer",    default=30,                  help='Number of LSI dimensions')
p$add_argument('--batch_variable',  type="character",                               help='Metadata column to apply batch correction on')
p$add_argument('--batch_method',    type="character",  default="MNN",               help='Batch correctin method ("Harmony" or "MNN")')
p$add_argument('--n_neighbors',     type="integer",    default=30,   nargs='+',     help='(UMAP) Number of neighbours')
p$add_argument('--min_dist',        type="double",     default=0.3,  nargs='+',     help='(UMAP) Minimum distance')
p$add_argument('--colour_by',       type="character",  nargs='+',  help='Metadata columns to colour the UMAP by')
p$add_argument('--outdir',          type="character",                               help='Output directory')
p$add_argument("--sort_samples",  type="logical",default=TRUE, help="Should samples be sorted?")
p$add_argument("--incl_samples",  type="character",    help='Which samples should be included')
args <- p$parse_args(commandArgs(TRUE))

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

print(head(sample_metadata.to.archr))

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

# Correlation between latent dimensions and number of peaks
# cor(getReducedDims(ArchRProject.filt, "IterativeLSI"),ArchRProject.filt$nFrags)[,1] %>% abs %>% sort(decreasing = T)

#########
## LSI ##
#########

outfile <- sprintf("%s/lsi_nfeatures%d_ndims%d_%s.txt.gz",args$outdir, args$nfeatures, args$ndims, args$incl_samples)

lsi.dt <- getReducedDims(ArchRProject.filt, "IterativeLSI") %>% round(3) %>% 
  as.data.table(keep.rownames = T) %>% setnames("rn","cell")

# Save LSI coordinates
fwrite(lsi.dt, outfile)

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
fwrite(umap.dt,  sprintf("%s/umap_nfeatures%d_ndims%d_%s.txt.gz",args$outdir, args$nfeatures, args$ndims, args$incl_samples))

# Plot
to.plot <- umap.dt %>%
  merge(sample_metadata,by="cell")

for (k in args$colour_by) {
  # log10 large numeric values
  if (is.numeric(to.plot[[k]])&& k != "RNA_cluster") {
    if (max(to.plot[[k]],na.rm=T) - min(to.plot[[k]],na.rm=T) > 1000) {
      to.plot[[k]] <- log10(to.plot[[k]]+1)
      to.plot %>% setnames(k,paste0(k,"_log10")); k <- paste0(k,"_log10")
    }
  }
  
  if (args$sort_samples){
    to.plot$sample <- as.factor(to.plot$sample)
    to.plot$sample <- factor(to.plot$sample,levels=args$samples)
  }
  
  print(dim(to.plot))
  print(table(to.plot$sample))
  
  if (k == "RNA_cluster") {
    to.plot[[k]] <- factor(to.plot[[k]])
  }
  p <- ggplot(to.plot, aes_string(x="umap1", y="umap2",fill=k)) +
    geom_point(size=pt.size, shape=21, stroke=0.05, color="transparent") +
    theme_classic() +
    ggplot_theme_NoAxes() +
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  # Define colormap for RNA_cluster and other categorical variables
  if (k == "RNA_cluster" && !is.null(opts$color_scheme)) {
    # Define a custom color palette for RNA_cluster with distinct colors
    distinct_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#F781BF", "#A65628")  # Red, Blue, Green, Orange, Pink, Brown
    p <- p + scale_fill_manual(values = distinct_colors)
  } else if ((grepl("sample", k) || k == "RNA_cluster") && !is.null(opts$color_scheme)) {
    # Use user-provided color scheme for sample or RNA_cluster
    p <- p + scale_fill_manual(values = opts$color_scheme[1:length(args$samples)])
  }
  
  # Save UMAP plot
  output_plot(p, sprintf("%s/umap_nfeatures%d_ndims%d_%s_%s",args$outdir, args$nfeatures, args$ndims, k, args$incl_samples),
              width=7, height=5, UMAP=TRUE)
}
