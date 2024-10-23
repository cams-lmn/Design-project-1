####################
##                ##
##  clustering.R  ##
##                ##
####################

source("/data/homes/louisc/Project_Babraham/RNA_ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--mofa_model',        type="character",                               help='MOFA model')
p$add_argument('--matrix',          type="character",  default="PeakMatrix",   help='Matrix to use')
p$add_argument('--samples',         type="character",                nargs='+',     help='Samples')
p$add_argument('--sort_samples',  default=TRUE,  help='Sort samples?')
p$add_argument('--incl_samples',        type="character",                               help='Suffix')
p$add_argument('--batch_correction',  default=FALSE,  help='Correct for batch effect?')
p$add_argument('--n_factors',       type="integer",  help='Amount of factors to use for clustering')
p$add_argument('--am_clusters',       type="integer",  help='Amount of clusters to create')
p$add_argument('--reorder_clusters',       type="integer",  default= NULL, nargs='+',  help='New order for clusters')
p$add_argument('--seed1',            type="integer",    default=42,                  help='Random seed 1')
p$add_argument('--seed2',            type="integer",    default=4242,                  help='Random seed 2')
p$add_argument('--outdir',          type="character",                               help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

print(args)

## START TEST ##
# args <- list()
# args$filter_differentiated <- TRUE
# if (args$filter_differentiated){
#   args$metadata <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/mofa/sample_metadata_nodiff.txt.gz"
#   args$mofa_model <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/mofa/mofa_nodiff.rds"
# } else {
#   args$metadata <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/mofa/sample_metadata.txt.gz"
#   args$mofa_model <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/mofa/mofa.rds"
# }
# args$batch_correction <- FALSE
# args$outdir <-"/data/homes/louisc/Project_Babraham/RNA_ATAC/clustering/PeakMatrix/"
# args$samples <- c("d0","d1","d3","d5","d7","d10","d14","d18","DE","NE","AME_E","AME_L")
# args$sort_samples <- TRUE
# args$incl_samples <- "nodiff"
# args$reorder_clusters <- c(3, 5, 2, 6, 1, 4)
# args$am_clusters <- 6
# args$n_factors <- 30
# args$seed1 <- 42
# args$seed2 <- 4242
## END TEST ##

#####################
## Define settings ##
#####################

# I/O
dir.create(args$outdir, showWarnings = F, recursive = T)

if(grepl("nodiff",args$incl_samples)){
  args$filter_differentiated <- TRUE
} else {
  args$filter_differentiated <- FALSE
}

if (args$filter_differentiated){
  args$samples <- args$samples[grepl("d[0-9]+",args$samples)]
}

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) 
print(head(sample_metadata))

if (args$batch_correction){
  sample_metadata$SeqRun <- "Run2"
  sample_metadata$SeqRun[sample_metadata$batch=="d0"] <- "Jasmin"
  sample_metadata$SeqRun[sample_metadata$batch%in%c("d10","d14","d18","DE")] <- "Run1"
}

###############
## Load MOFA ##
###############

if (tools::file_ext(args$mofa_model)=="rds") {
  MOFAobject <- readRDS(args$mofa_model)
  
  # Add sample metadata
  samples_metadata(MOFAobject) <- sample_metadata
} else if (file_ext(args$mofa_model)=="hdf5") {
  MOFAobject <- load_model(args$mofa_model, load_data = F)
  
  # Add sample metadata
  cells <- as.character(unname(unlist(MOFA2::samples_names(MOFAobject))))
  sample_metadata_to_mofa <- copy(sample_metadata) %>%
    setnames("cell","sample") %>%
    .[sample%in%cells] %>% setkey(sample) %>% .[cells]
  stopifnot(all(cells==sample_metadata_to_mofa$cell))
  samples_metadata(MOFAobject) <- sample_metadata_to_mofa
}

##################
## Load factors ##
##################

# Select factors to use 
factors.to.use <- 1:args$n_factors
# factors.to.use <- factors.to.use[!factors.to.use%in%c("3")]
# factors.to.use <- c(1,2,3,6,7)

# Extract factors
Z <- get_factors(MOFAobject, factors=factors.to.use)[[1]]

print(factors.to.use)
print(head(Z))

if (args$batch_correction) {
  # Define stage and sample order
  timepoints <- MOFAobject@samples_metadata$SeqRun
  timepoint_order <- unique(MOFAobject@samples_metadata$SeqRun)
  samples <- MOFAobject@samples_metadata$batch
  sample_order <- args$samples
  
  Z_list    <- lapply(unique(timepoints), function(i){
    sub_pc   <- Z[timepoints == i, , drop = FALSE]
    sub_samp <- samples[timepoints == i]
    list     <- lapply(unique(sub_samp), function(j){ sub_pc[sub_samp == j, , drop = FALSE]})
    names(list) <- unique(sub_samp)
    return(list)
  })
  names(Z_list) <- unique(timepoints)
  
  # #arrange to match timepoint order
  Z_list <- Z_list[order(match(names(Z_list), timepoint_order))]
  Z_list <- lapply(Z_list, function(x){ x[order(match(names(x), sample_order))]})
  
  #perform corrections within stages
  correct_list <- lapply(Z_list, function(x){
    if(length(x) > 1){
      return(do.call(reducedMNN, x)$corrected)
    } else {
      return(x[[1]])
    }
  })
  
  #perform correction over stages
  Z_corrected <- reducedMNN(correct_list, merge.order=1:length(correct_list))$corrected 
  # colnames(Z_corrected) <- colnames(Z)
} else {
  Z_corrected <- Z
}

print(head(Z_corrected))

# Save factors
factors.dt <- Z_corrected %>% round(3) %>% as.data.table(keep.rownames = T) %>% setnames("rn","cell")
fwrite(factors.dt, sprintf("%s/factors_%s.txt.gz",args$outdir,args$incl_samples))

##########
## UMAP ##
##########

# Run
set.seed(args$seed1)
umap_embedding <- uwot::umap(Z_corrected, n_neighbors=25, min_dist=0.50, metric="cosine")

################
## Clustering ##
################

set.seed(args$seed2)
MOFA_clusters <- cluster_samples(MOFAobject, k = args$am_clusters, factors = 1:args$n_factors)

if(!is.null(args$reorder_clusters)){
  names_cluster <- names(MOFA_clusters$cluster)
  MOFA_clusters$cluster <- args$reorder_clusters[MOFA_clusters$cluster]
  names(MOFA_clusters$cluster) <- names_cluster
}

# fix different mofa
fix_data <- fread("/data/homes/louisc/Project_Babraham/RNA_ATAC/clustering/PeakMatrix/sample_metadata_nodiff_after_clustering.txt")
MOFA_clusters$cluster <- fix_data$cluster

##########
## Plot ##
##########

# Plot
to.plot <- umap_embedding %>% as.data.table %>%
  .[,sample:=rownames(Z_corrected)] %>%
  merge(MOFAobject@samples_metadata[,c("sample","batch")] %>% as.data.table)

sum(names(MOFA_clusters$cluster)==to.plot$sample)
dim(to.plot)

to.plot$cluster <- as.factor(MOFA_clusters$cluster)

if (args$sort_samples){
  to.plot$batch <- factor(to.plot$batch,levels=args$samples)
}

print(dim(to.plot))
print(table(to.plot$cluster))

print("UMAP")
print(head(to.plot))
fwrite(to.plot, sprintf("%s/umap_clusters_%s.txt.gz",args$outdir,args$incl_samples))

p <- ggplot(to.plot, aes(x=V1, y=V2, fill=cluster)) +
  geom_point(size=1, shape=21, stroke=0.1) +
  # ggrastr::geom_point_rast(size=1, shape=21, stroke=0.1) +
  scale_fill_manual(values=opts$celltype.colors) +
  guides(fill = guide_legend(override.aes = list(size=2))) +
  theme_classic() +
  # theme(legend.position="none") +
  ggplot_theme_NoAxes() 

# if (!is.null(args$color_scheme)){
#   p <- p + scale_fill_manual(values=args$color_scheme[1:length(args$samples)])
# }

pdf(sprintf("%s/mofa_umap_cluster_%s.pdf",args$outdir,args$incl_samples), width=7, height=5)
print(p)
dev.off()

###################
## Save metadata ##
###################

sum(names(MOFA_clusters$cluster)==sample_metadata$sample)
dim(sample_metadata)

sample_metadata$cluster <- MOFA_clusters$cluster

# Save
fwrite(sample_metadata, sprintf("%s/sample_metadata_%s_after_clustering.txt.gz",args$outdir,args$incl_samples), 
       quote=F, na="NA", sep="\t")

