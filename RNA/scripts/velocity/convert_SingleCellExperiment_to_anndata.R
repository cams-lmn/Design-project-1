###############################################
##                                           ##
##  convert_SingleCellExperiment_to_anndata  ##
##                                           ##
###############################################

source("/data/homes/louisc/Project_Babraham/RNA/scripts/Settings.R")

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--python_path',   type="character",    help='Python path for reticulate')
p$add_argument('--metadata',   type="character",    help='Cell metadata')
p$add_argument('--sce',  type="character",              help='SingleCellExperiment input file')
p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
p$add_argument('--incl_samples',        type="character",                               help='Suffix')
p$add_argument('--sort_samples',  type="logical",default=TRUE, help="Should samples be sorted?")
p$add_argument('--outdir',          type="character",                help='Anndata output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# args$python_path <- " /opt/conda/envs/Babraham/bin/python" # "/Users/argelagr/opt/anaconda3/envs/main/bin/python"
# args$metadata <- "/data/louisc/Project_Babraham/RNA/mapping/sample_metadata_after_doublets.txt.gz"
# args$sce <- "/data/louisc/Project_Babraham/RNA/SingleCellExperiment.rds"
# args$filter_differentiated <- TRUE
# args$samples <- c("d0","d1","d3","d5","d7","d10","d14","d18","DE","NE","AME_E","AME_L")
# args$sort_samples <- TRUE
## END TEST ##

#####################
## Define settings ##
#####################

if(grepl("nodiff",args$incl_samples)){
  args$filter_differentiated <- TRUE
} else {
  args$filter_differentiated <- FALSE
}

args$outfile <- paste0(args$outdir,"/anndata_",args$incl_samples,".h5ad")
print(args$outfile)

################
## Reticulate ##
################

print(args$python_path) # reticulate should be version 1.22 (https://ftp.belnet.be/mirror/CRAN/src/contrib/00Archive/reticulate/reticulate_1.22.tar.gz)
reticulate::use_python(args$python_path, required = TRUE)
ad <- import("anndata")
print("anndata loaded")

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>% 
  .[pass_rnaQC==TRUE & doublet_call==FALSE]# %>%
# .[,c("cell", "sample", "stage", "nFeature_RNA", "mitochondrial_percent_RNA", "ribosomal_percent_RNA", "celltype.mapped")] %>%
# setnames("celltype.mapped","celltype")

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

###########
## Parse ##
###########

# reducedDims
# reducedDimNames(sce)[1] <- "X_pca_precomputed"
# forceatlas.mtx <- colData(sce)[,c("cell","Dim1","Dim2")] %>% matrix.please
# reducedDims(sce)[["X_draw_graph_fa_precomputed"]] <- forceatlas.mtx

# colData
# colData(sce) <- colData(sce)[,c("cell","Dim1","Dim2","sample", "stage.mapped", "somite_count", "tube_name", "tube_name_corrected", "celltype.clustering")]
# colnames(colData(sce)) <- c("cell","x","y","sample", "stage", "somite_count", "tube_name", "tube_name_corrected", "celltype")


#############################################
## Convert SingleCellExperiment to AnnData ##
#############################################

# adata_sce <- sc$AnnData(
#     X   = Matrix::t(counts(sce)),
#     obs = as.data.frame(colData(sce)),
#     obsm = list(
#       "X_pca_precomputed" = as.matrix(reducedDim(sce, "X_pca_precomputed")),
#       "X_draw_graph_fa_precomputed" = as.matrix(reducedDim(sce, "X_draw_graph_fa_precomputed"))
#     ),
#     var = data.frame(gene=rownames(sce), row.names=rownames(sce))
# )

adata_sce <- ad$AnnData(
  X   = t(counts(sce)),
  obs = as.data.frame(colData(sce)),
  var = data.frame(gene=rownames(sce), row.names=rownames(sce))
)

adata_sce

# TO-DO: CAN WE MAKE THE MATRICES SPARSE??

##########################
## Parse anndata object ##
##########################

# # Add stage colors
# adata_sce$uns$update(stage_colors = opts$stage.colors[sort(unique(as.character(adata_sce$obs$stage)))])
# adata_sce$uns["stage_colors"]
# 
# # Add celltype colors
# adata_sce$uns$update(celltype_colors = opts$celltype.colors[sort(unique(as.character(adata_sce$obs$celltype)))])
# adata_sce$uns["celltype_colors"]

##########
## Save ##
##########

adata_sce$write_h5ad(args$outfile, compression="gzip")
