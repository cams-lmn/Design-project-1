######################
##                  ##
##  Prepare_mofa.R  ##
##                  ##
######################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA_ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",   help='Cell metadata file')
p$add_argument('--atac_matrix',          type="character",    help='ATAC Matrix to use')
p$add_argument('--atac_matrix_file',          type="character",   help='ATAC Matrix to use')
p$add_argument('--atac_feature_stats',          type="character",   help='ATAC feature stats')
p$add_argument('--sce',          type="character",   help='RNA SingleCellExperiment')
p$add_argument('--samples',       type="character",  default="all",  nargs='+',  help='Stages to plot')
p$add_argument('--nfeatures_atac',       type="integer",    default=25000,               help='Number of ATAC features')
p$add_argument('--nfeatures_rna',       type="integer",    default=5000,               help='Number of RNA features')
p$add_argument('--outdir',          type="character",                               help='Output directory')
p$add_argument('--sort_samples',  default=TRUE,  help='Sort samples?')
args <- p$parse_args(commandArgs(TRUE))

dir.create(args$outdir, showWarnings=F, recursive = T)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE & doublet_call==FALSE & sample%in%args$samples]

# print stats
if (args$sort_samples){
  table(sample_metadata$sample)[args$samples]
} else {
  table(sample_metadata$sample)
}

####################
## Load ATAC data ##
####################

print(sprintf("Fetching single-cell ATAC %s...",args$atac_matrix))

atac.se <- readRDS(args$atac_matrix_file)[,sample_metadata$cell]

############################
## Feature selection ATAC ##
############################

# Load feature stats
atac_feature_stats.dt <- data.table(
  feature = rownames(atac.se),
  var_cells = sparseMatrixStats::rowVars(assay(atac.se)) %>% round(5),
  mean_cells = Matrix::rowMeans(assay(atac.se)) %>% round(5)
) %>% setkey(feature)

# Define highly variable features
atac_features <- atac_feature_stats.dt %>% 
  setorder(-var_cells) %>% 
  head(n=args$nfeatures_atac) %>% .$Feature

#############################
## ATAC data normalisation ##
#############################

# TFIDF normalisation
atac_tfidf.mtx <- tfidf(assay(atac.se[atac_features,]), method=1, scale.factor=1e4) %>% round(2)
dim(atac_tfidf.mtx)

##############################
## Load RNA expression data ##
##############################

print("Fetching single-cell RNA %s...")

rna.sce <- load_SingleCellExperiment(args$sce, normalise = TRUE, cells = sample_metadata$cell)

# Add sample metadata to the colData of the SingleCellExperiment
colData(rna.sce) <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(rna.sce),] %>% DataFrame()

# Filter features manually
rna.sce <- rna.sce[grep("^MT-|^RPS|^RPL|^OR[0-9]",rownames(rna.sce), invert=T),]

# Feature selection
rna_features <- modelGeneVar(rna.sce) %>% as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  .[mean>=0.001] %>% setorder(-p.value) %>% head(args$nfeatures_rna) %>% .$gene

# Fetch matrix
rna.mtx <- assay(rna.sce[rna_features,],"logcounts") %>% round(2)
dim(rna.mtx)

###########################
## Prepare data for MOFA ##
###########################

cells <- intersect(colnames(rna.mtx),colnames(atac_tfidf.mtx))

rna.mtx <- rna.mtx[,cells]
atac_tfidf.mtx <- atac_tfidf.mtx[,cells]

###############
## Save data ##
###############

Matrix::writeMM(rna.mtx, file.path(args$outdir,"rna.mtx"))
Matrix::writeMM(atac_tfidf.mtx, file.path(args$outdir,"atac_tfidf.mtx"))
write.table(rownames(rna.mtx), file.path(args$outdir,"rna_features.txt"), quote=F, row.names=F, col.names=F)
write.table(rownames(atac_tfidf.mtx), file.path(args$outdir,"atac_features.txt"), quote=F, row.names=F, col.names=F)
write.table(cells, file.path(args$outdir,"cells.txt"), quote=F, row.names=F, col.names=F)
fwrite(sample_metadata, file.path(args$outdir,"/sample_metadata.txt.gz"))