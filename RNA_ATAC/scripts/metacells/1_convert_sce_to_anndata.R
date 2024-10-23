###############################
##                           ##
##  Convert_sce_to_anndata.R ##
##                           ##
###############################

source("/data/homes/louisc/Project_Babraham/RNA_ATAC/scripts/Settings.R")

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

# ## START TEST ##
# args <- list()
# args$python_path <- "/opt/conda/envs/ribocode/bin/python"
# args$metadata <- "/data/homes/louisc/Project_Babraham/ATAC/archR/qc/sample_metadata_after_qc.txt.gz"
# args$sce <- "/data/homes/louisc/Project_Babraham/RNA/SingleCellExperiment.rds"
# args$filter_differentiated <- TRUE
# args$samples <- c("d0","d1","d3","d5","d7","d10","d14","d18","DE","NE","AME_E","AME_L")
# args$sort_samples <- TRUE
# if (args$filter_differentiated){
#   args$outfile <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/anndata_nodiff.h5ad"
# } else {
#   args$outfile <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/anndata.h5ad"
# }
# ## END TEST ##

#####################
## Define settings ##
#####################

if(grepl("nodiff",args$incl_samples)){
  args$filter_differentiated <- TRUE
} else {
  args$filter_differentiated <- FALSE
}

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

# TO-DO: CAN WE MAKE THE MATRICES SPARSE??

##########
## Save ##
##########

adata_sce$write_h5ad(args$outfile, compression="gzip")

