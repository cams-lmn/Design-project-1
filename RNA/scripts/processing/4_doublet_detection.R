###########################
##                       ##
##  Doublet_detection.R  ##
##                       ##
###########################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',         type="character",                            help='SingleCellExperiment file')
p$add_argument('--metadata',    type="character",                            help='metadata file')
p$add_argument('--samples',                 type="character",   nargs='+',     help='Sample(s)')
p$add_argument('--outfile',                  type="character",                  help='Output file')
p$add_argument('--seed',                  type="integer",  default=42,                help='Seed')
p$add_argument("--sort_samples",  type="logical",default=TRUE, help="Should samples be sorted?")
args <- p$parse_args(commandArgs(TRUE))

# Parse arguments
dir.create(dirname(args$outfile))

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & sample%in%args$samples]
if (args$sort_samples){
  table(sample_metadata$sample)[args$samples]
} else {
  table(sample_metadata$sample)
}

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- load_SingleCellExperiment(args$sce, cells=sample_metadata$cell, normalise = TRUE)
dim(sce)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

#############################
## Calculate doublet score ##
#############################

set.seed(args$seed)

# sce <- cxds_bcds_hybrid(sce, estNdbl=TRUE) alternative
sce <- scDblFinder(sce, samples="sample")

dt <- colData(sce) %>%
  .[,c("sample","scDblFinder.score", "scDblFinder.class")] %>%
  as.data.frame %>% tibble::rownames_to_column("cell") %>% as.data.table %>%
  .[,c("scDblFinder.score", "scDblFinder.class"):=list(round(scDblFinder.score,2), scDblFinder.class)]

# Call doublets
dt[,doublet_call:= scDblFinder.class == "doublet"]
table(dt$doublet_call)

# Save
fwrite(dt, args$outfile, sep="\t", na="NA", quote=F)