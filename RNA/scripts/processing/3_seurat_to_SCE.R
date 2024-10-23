#######################
##                   ##
##  Seurat_to_SCE.R  ##
##                   ##
#######################

source("/data/louisc/Project_Babraham/RNA/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--test',          type="logical", default=FALSE,                 help='Testing mode')
p$add_argument('--normalise',        type="logical", default=TRUE,                   help='Log-Normalise?')
p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
p$add_argument('--seurat',         type="character", help='Seurat object (input)')
p$add_argument('--metadata',         type="character", help='Metadata file')
p$add_argument('--metadata_out',         type="character", help='Metadata output file')
p$add_argument('--outfile',         type="character", help='Output file')
p$add_argument('--sort_samples',       type="logical",  default=TRUE,       help='Sort samples?')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# args$outfile <- "/data/louisc/Project_Babraham/RNA/SingleCellExperiment.rds"
# args$samples <- c("d0","d1","d3","d5","d7","d10","d14","d18","DE","NE","AME_E","AME_L")
# args$metadata <- "/data/louisc/Project_Babraham/RNA/qc/sample_metadata_after_qc.txt.gz"
# args$metadata_out <- "/data/louisc/Project_Babraham/RNA/qc/sample_metadata_after_cellstate.txt.gz"
# args$seurat <- "/data/louisc/Project_Babraham/RNA/seurat.rds"
# args$test <- FALSE
# args$normalise <- TRUE
# args$sort_samples <- TRUE
## END TEST ##

#####################
## Define settings ##
#####################

# Sanity checks
# stopifnot(args$samples%in%opts$samples)
# if (args$test) args$samples <- head(args$samples,n=2)

###############
## Load data ##
###############

# Load sample metadata
sample_metadata <- fread(args$metadata) %>% .[pass_rnaQC==TRUE & sample%in%args$samples]
if (args$sort_samples){
  table(sample_metadata$sample)[args$samples]
} else {
  table(sample_metadata$sample)
}

# Load seurat
seurat <- readRDS(args$seurat)[,sample_metadata$cell]

#####################################
## Convert to SingleCellExperiment ##
#####################################

sce <- as.SingleCellExperiment(seurat)

# remove logcounts assays
sce@assays@data[["logcounts"]] <- NULL

# Add metadata
# stopifnot(sample_metadata$cell%in%colnames(sce))
# stopifnot(colnames(sce)%in%sample_metadata$cell)
sample_metadata <- sample_metadata %>% .[cell%in%colnames(sce)] %>% setkey(cell) %>% .[colnames(sce)]
stopifnot(sample_metadata$cell == colnames(sce))
colData(sce) <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(sce),] %>% DataFrame()

##########################
## Compute size factors ##
##########################

clusts <- as.numeric(quickCluster(sce, method = "igraph", min.size = 100, BPPARAM = mcparam))
# clusts <- as.numeric(quickCluster(sce))
min.clust <- min(table(clusts))/2
new_sizes <- c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
sce <- computeSumFactors(sce, clusters = clusts, sizes = new_sizes, max.cluster.size = 3000)

###################
## Log Normalise ##
###################

if (args$normalise) {
  sce <- logNormCounts(sce)
}

####################
## Add cell state ##
####################

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

seurat <- CellCycleScoring(seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
sample_metadata$Phase <- seurat[[]]$Phase

fwrite(sample_metadata,args$metadata_out, quote=F, na="NA", sep="\t")

##########
## Plot ##
##########

# to.plot <- data.frame(X = Matrix::colSums(counts(sce)), Y = sizeFactors(sce))
# ggplot(to.plot, mapping = aes(x = X, y = Y)) +
#   geom_point() +
#   labs(x = "Number of UMIs", y = "Size Factor") +
#   theme_classic()

##########
## Save ##
##########

saveRDS(sce, args$outfile)

