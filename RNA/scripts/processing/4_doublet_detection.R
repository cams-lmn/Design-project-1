###########################
##                       ##
##  Doublet_detection.R  ##
##                       ##
###########################

source("/data/louisc/Project_Babraham/RNA/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',         type="character",                            help='SingleCellExperiment file')
p$add_argument('--metadata',    type="character",                            help='metadata file')
p$add_argument('--samples',                 type="character",   nargs='+',     help='Sample(s)')
p$add_argument('--doublet_score_threshold',  type="double",      default=1.25,   help='Doublet score threshold')
p$add_argument('--outfile',                  type="character",                  help='Output file')
p$add_argument('--seed',                  type="integer",  default=42,                help='Seed')
p$add_argument("--sort_samples",  type="logical",default=TRUE, help="Should samples be sorted?")
p$add_argument("--test",  type="logical",default=FALSE, help="Test?")
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# args$sce <-"/data/louisc/Project_Babraham/RNA/SingleCellExperiment.rds"
# args$metadata <- "/data/louisc/Project_Babraham/RNA/qc/sample_metadata_after_cellstate.txt.gz" # .io$metadata
# args$samples <- c("d0","d1","d3","d5","d7","d10","d14","d18","DE","NE","AME_E","AME_L")
# args$doublet_score_threshold <- 1.25
# args$test <- TRUE
# args$seed <- 42
# args$outfile <- sprintf("/data/louisc/Project_Babraham/RNA/mapping/doublets_all_%s.txt.gz",round(args$doublet_score_threshold,2))
# args$sort_samples <- TRUE
## END TEST ##

# Parse arguments
dir.create(dirname(args$outfile))
if (isTRUE(args$test)) print("Test mode activated...")

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

sce <- cxds_bcds_hybrid(sce, estNdbl=TRUE)

dt <- colData(sce) %>%
  .[,c("sample","cxds_score", "bcds_score", "hybrid_score")] %>%
  as.data.frame %>% tibble::rownames_to_column("cell") %>% as.data.table %>%
  .[,c("cxds_score","bcds_score","hybrid_score"):=list(round(cxds_score,2),round(bcds_score,2),round(hybrid_score,2))]

# Call doublets
dt[,doublet_call:=hybrid_score>args$doublet_score_threshold]
table(dt$doublet_call)

# Save
# io$outfile <- sprintf("%s/doublets_%s_%s.txt.gz",args$outdir, paste(args$samples,collapse="-"),round(args$doublet_score_threshold,2))
fwrite(dt, args$outfile, sep="\t", na="NA", quote=F)
