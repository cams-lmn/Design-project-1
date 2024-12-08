########################################
##                                    ##
##  archR_dimensionality_reduction.R  ##
##                                    ##
########################################

source("/data/homes/louisc/Project_Babraham/ATAC/scripts/Settings.R")

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
p$add_argument('--genome',          type="character", default="hg38",      help='Genome')
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

## START TEST ##
# args <- list()
# args$metadata <- "/data/homes/louisc/Project_Babraham/ATAC/archR/qc/sample_metadata_after_qc.txt.gz"
# args$nfeatures <- 15000
# args$matrix <- "PeakMatrix"
# args$sort_samples <- TRUE
# args$ndims <- 25
# args$seed <- 42
# args$n_neighbors <- 25
# args$min_dist <- 0.25
# args$genome <- "hg38"
# args$batch_variable <- "SeqRun"
# args$batch_method <- "Harmony"
# args$threads <- 2
# args$max_point_per_sample <- NULL
# args$seed3 <- 424242
# args$filter_differentiated <- FALSE
# args$remove_ExE_cells <- FALSE
# args$samples <- c("d0","d1","d3","d5","d7","d10","d14","d18")
# args$colour_by <- c("sample","log_nFrags_atac","TSSEnrichment_atac","Phase")
# args$vars_to_regress <- c("nFeature_RNA","mitochondrial_percent_RNA")
# args$outdir <- "/data/louisc/Project_Babraham/ATAC/archR/dimensionality_reduction"
# args$archr_directory <- "/data/louisc/Project_Babraham/ATAC/archR"
# args$color_scheme <- c("#A06932","#DBB216","#EFE32A","#D7DB54","#A7B019","#9AA126",
#                        "#7E7721","#5D6821","#353D8A","#CA3639","#DCADCD","#7A327E")
# args$incl_samples <- "nodiff"
## END TEST ##


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

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  #.[ pass_atacQC==TRUE & doublet_call==FALSE & sample%in%args$samples] %>%
  # .[nFrags_atac>=3500 & TSSEnrichment_atac>=9,pass_atacQC:=TRUE] %>% .[stage=="E8.7",stage:="E8.75"] %>%    # temporary
  .[pass_atacQC==TRUE & doublet_call==FALSE] %>%
  .[,log_nFrags_atac:=log10(nFrags_atac)]

# stopifnot(args$colour_by %in% colnames(sample_metadata))

# if (args$remove_ExE_cells) {
#   print("Removing ExE cells...")
#   sample_metadata <- sample_metadata %>%
#     .[!celltype%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
# }

# table(sample_metadata$stage)
# table(sample_metadata$celltype)

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

print(head(sample_metadata.to.archr))
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

# Correlation between latent dimensions and number of peaks
# cor(getReducedDims(ArchRProject.filt, "IterativeLSI"),ArchRProject.filt$nFrags)[,1] %>% abs %>% sort(decreasing = T)

############################
## LSI + Batch correction ##
############################

if (length(args$batch_variable)>0) {
  print(sprintf("Applying %s batch correction for variable: %s", args$batch_method, args$batch_variable))
  outfile <- sprintf("%s/lsi_nfeatures%d_ndims%d_%s.txt.gz",
                     args$outdir, args$nfeatures, args$ndims, args$incl_samples)
  
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
  outfile <- sprintf("%s/lsi_nfeatures%d_ndims%d_%s.txt.gz",args$outdir, args$nfeatures, args$ndims, args$incl_samples)
  
  lsi.dt <- getReducedDims(ArchRProject.filt, "IterativeLSI") %>% round(3) %>% 
    as.data.table(keep.rownames = T) %>% setnames("rn","cell")
}

# Save LSI coordinates
fwrite(lsi.dt, outfile)

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
fwrite(umap.dt,  sprintf("%s/umap_nfeatures%d_ndims%d_%s.txt.gz",args$outdir, args$nfeatures, args$ndims, args$incl_samples))

# Plot
to.plot <- umap.dt %>%
  merge(sample_metadata,by="cell")

for (k in args$colour_by) {
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
  
  print(dim(to.plot))
  print(table(to.plot$sample))
  
  p <- ggplot(to.plot, aes_string(x="umap1", y="umap2",fill=k)) +
    geom_point(size=pt.size, shape=21, stroke=0.05, color="transparent") +
    # ggrastr::geom_point_rast(size=1.5, shape=21, stroke=0.05) +  # DOES NOT WORK IN THE CLUSTER
    theme_classic() +
    ggplot_theme_NoAxes() +
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  # Define colormap
  if (grepl("sample",k) & !is.null(opts$color_scheme)) {
    p <- p + scale_fill_manual(values=opts$color_scheme[1:length(args$samples)]) 
  }
  
  # Save UMAP plot
  output_plot(p, sprintf("%s/umap_nfeatures%d_ndims%d_%s_%s",args$outdir, args$nfeatures, args$ndims, k, args$incl_samples),
              width=7, height=5, UMAP=TRUE)
}

###################
## Plot animated ##
###################
# 
# to.plot <- umap.dt %>%
#   merge(sample_metadata,by="cell")
# 
# if (is.numeric(to.plot[[k]])) {
#   if (max(to.plot[[k]],na.rm=T) - min(to.plot[[k]],na.rm=T) > 1000) {
#     to.plot[[k]] <- log10(to.plot[[k]]+1)
#     to.plot %>% setnames(k,paste0(k,"_log10")); k <- paste0(k,"_log10")
#   }
# }
# 
# if (args$sort_sample){
#   to.plot$sample <- as.factor(to.plot$sample)
#   to.plot$sample <- factor(to.plot$sample,levels=args$samples)
# }
# 
# to.plot.elab <- NULL
# for (i in levels(to.plot$sample)){
#   to.plot.elab.tmp <- to.plot 
#   to.plot.elab.tmp$sample <- as.factor(c("other",i)[c(to.plot.elab.tmp$sample==i)+1])
#   to.plot.elab.tmp$plot <- i
#   to.plot.elab.tmp$sample <- factor(to.plot.elab.tmp$sample,levels=c("other",i))
#   to.plot.elab.tmp <- to.plot.elab.tmp[sort(as.character(to.plot.elab.tmp$sample),decreasing = T,index.return=T)$ix,]
#   to.plot.elab <- rbind(to.plot.elab,to.plot.elab.tmp)
# }
# if (args$sort_samples){
#   to.plot.elab$plot <- factor(to.plot.elab$plot,levels=args$samples)
# }
# p <- ggplot(to.plot.elab, aes_string(x="umap1", y="umap2",fill="sample")) +
#   geom_point(size=pt.size, shape=21, stroke=0.05) +
#   theme_classic() +
#   transition_states(plot,0.25,8) +
#   labs(title = 'Timepoint: {closest_state}') +
#   scale_fill_manual(values=c("black",rep("red",length(args$samples)))) +
#   theme(legend.position="none") +
#   ggplot_theme_NoAxes()
# 
# if (args$filter_differentiated){
#   anim_save(sprintf("%s/umap_%s_nfeatures%d_ndims%d_sample_nodiff.gif",args$outdir, args$matrix, args$nfeatures, args$ndims),
#             animation = (animate(p, renderer = gifski_renderer())))
# } else {
#   anim_save(sprintf("%s/umap_%s_nfeatures%d_ndims%d_sample.gif",args$outdir, args$matrix, args$nfeatures, args$ndims),
#             animation = (animate(p, renderer = gifski_renderer())))
# }
