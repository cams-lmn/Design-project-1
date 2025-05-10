###############################
##                           ##
##  Convert_sce_to_anndata.R ##
##                           ##
###############################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA_ATAC/scripts/Settings.R")

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--python_path',   type="character",    help='Python path for reticulate')
p$add_argument('--metadata',   type="character",    help='Cell metadata')
p$add_argument('--sce',  type="character",              help='SingleCellExperiment input file')
p$add_argument('--samples',  type="character", nargs="+",      help='SingleCellExperiment input file')
p$add_argument('--outfile',          type="character",                help='Anndata output file')
p$add_argument('--sort_samples',  default=TRUE,  help='Sort samples?')
p$add_argument('--incl_samples',  type="character",     help='Which samples to include')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

args$filter_differentiated <- FALSE

################
## Reticulate ##
################

reticulate::use_python(args$python_path, required = TRUE)
ad <- reticulate::import("anndata")

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE & doublet_call==FALSE & sample%in%args$samples]

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(args$sce, cells=sample_metadata$cell, normalise = FALSE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

#############################################
## Convert SingleCellExperiment to AnnData ##
#############################################

adata_sce <- ad$AnnData(
  X   = t(counts(sce)),
  obs = as.data.frame(colData(sce)),
  var = data.frame(gene=rownames(sce), row.names=rownames(sce))
)

adata_sce

##########
## Save ##
##########

adata_sce$write_h5ad(args$outfile, compression="gzip")