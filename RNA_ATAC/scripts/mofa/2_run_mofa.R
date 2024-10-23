##################
##              ##
##  Run_mofa.R  ##
##              ##
##################

source("/data/homes/louisc/Project_Babraham/RNA_ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--atac_dimred',          type="character",  help='')
p$add_argument('--rna_dimred',          type="character",   help='')
p$add_argument('--incl_samples',        type="character",                               help='Suffix')
p$add_argument('--samples',       type="character",   nargs='+',  help='Samples to plot')
p$add_argument('--n_factors',           type="integer",    default=30,                  help='Number of MOFA factors')
p$add_argument('--outdir',          type="character",                               help='Output directory')
p$add_argument('--python_path',          type="character",                help='Python path')
p$add_argument('--seed',        type="integer", default=42,               help='Seed')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# # io$basedir <- file.path(io$basedir,"test")
# args$metadata <- "/data/homes/louisc/Project_Babraham/ATAC/archR/qc/sample_metadata_after_qc.txt.gz"
# args$atac_dimred <-"/data/homes/louisc/Project_Babraham/ATAC/archR/dimensionality_reduction/PeakMatrix/batch_correction_SeqRun_Harmony/lsi_nfeatures25000_ndims15_nodiff.txt.gz"
# args$rna_dimred <- "/data/homes/louisc/Project_Babraham/RNA/dimensionality_reduction/batch_correction_SeqRun/pca_features2500_pcs15_nodiff.txt.gz"
# args$samples <- c("d0","d1","d3","d5","d7","d10","d14","d18","DE","NE","AME_E","AME_L")
# args$factors <- 30
# args$filter_differentiated <- TRUE
# args$outdir <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/mofa/fast"
# args$python_path <- "/opt/conda/envs/ribocode/bin/python"
## END TEST ##

dir.create(args$outdir, showWarnings = F)
opts <- list()

if(grepl("nodiff",args$incl_samples)){
  args$filter_differentiated <- TRUE
} else {
  args$filter_differentiated <- FALSE
}

#####################
## Parse arguments ##
#####################

print(args)

########################
## Load cell metadata ##
########################

sample_metadata <- fread(args$metadata) %>%
  # .[(pass_atacQC==TRUE & pass_rnaQC==TRUE) & doublet_call==FALSE] %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE & doublet_call==FALSE]

if (args$filter_differentiated){
  sample_metadata <- sample_metadata[grepl("d[0-9]+",sample_metadata$sample),]
}

opts$rna.cells <- sample_metadata[pass_rnaQC==TRUE,cell]
opts$atac.cells <- sample_metadata[pass_atacQC==TRUE,cell]

###############################################
## Load precomputed dimensionality reduction ##
###############################################

# RNA (PCA)
tmp <- fread(args$rna_dimred)
opts$rna.cells <- intersect(tmp$cell,opts$rna.cells)
rna.mtx <- tmp %>% matrix.please %>% .[opts$rna.cells,] %>% t
rm(tmp)

# ATAC (LSI)
tmp <- fread(args$atac_dimred)
opts$atac.cells <- intersect(tmp$cell,opts$atac.cells)
atac.mtx <- tmp %>% matrix.please %>% .[opts$atac.cells,] %>% t
rm(tmp)

###########################
## Prepare data for MOFA ##
###########################

rna.mtx <- .augment_matrix(rna.mtx, unique(opts$rna.cells,opts$atac.cells))
atac.mtx <- .augment_matrix(atac.mtx, unique(opts$rna.cells,opts$atac.cells))

########################
## Create MOFA object ##
########################

MOFAobject <- create_mofa_from_matrix(list("RNA" = rna.mtx, "ATAC" = atac.mtx))
MOFAobject

####################
## Define options ##
####################

set.seed(args$seed)

# Data options
data_opts <- get_default_data_options(MOFAobject)
data_opts$use_float32 <- TRUE

# Model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- args$n_factors
model_opts$spikeslab_weights <- FALSE
model_opts$ard_weights <- FALSE

# Training options
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "fast"
# train_opts$maxiter <- 5

#########################
## Prepare MOFA object ##
#########################

MOFAobject <- prepare_mofa(
  MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

#####################
## Train the model ##
#####################

reticulate::use_python(args$python_path,require=T)
print(reticulate::py_config())
print(reticulate::py_module_available("mofapy2"))
set.seed(args$seed)
print(str(MOFAobject))
MOFAobject <- run_mofa(MOFAobject,use_basilisk = F)

#########################
# Add samples metadata ##
#########################

metadata.to.mofa <- sample_metadata %>% copy %>%
  setnames("sample","batch") %>% setnames("cell","sample") %>%
  .[sample%in%unlist(samples_names(MOFAobject))] %>%
  setkey(sample) %>% .[unlist(samples_names(MOFAobject))]
samples_metadata(MOFAobject) <- metadata.to.mofa

##########
## Save ##
##########

saveRDS(MOFAobject, file.path(args$outdir,sprintf("mofa_%s.rds",args$incl_samples)))
fwrite(metadata.to.mofa, file.path(args$outdir,sprintf("sample_metadata_%s.txt.gz",args$incl_samples)), quote=F, sep="\t", na="NA")

