########################
##                    ##
##  Pseudobulk_rna.R  ##
##                    ##
########################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA_ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',         type="character",    help='SingleCellExperiment object')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--group_by',    type="character",    help='Metadata column to group cells by')
p$add_argument('--normalisation_method',    type="character",    help='Metadata column to group cells by')
p$add_argument('--outdir',      type="character",    help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

######################
## Define settings ##
######################

dir.create(args$outdir, showWarnings = F, recursive = T)

###############
## Load data ##
###############

# Load cell metadata
sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE & !is.na(eval(as.name(args$group_by)))] %>%
  dplyr::rename(cell=sample)

# Load SingleCellExperiment
sce <- load_SingleCellExperiment(args$sce, cells=sample_metadata$cell)
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

################
## Pseudobulk ##
################

sce_pseudobulk <- pseudobulk_sce_fn(
  x = sce,
  assay = "counts",
  by = args$group_by,
  fun = "sum",
  scale = FALSE # Should pseudo-bulks be scaled with the effective library size & multiplied by 1M?
)

assayNames(sce_pseudobulk) <- "counts"

###################
## Normalisation ##
###################

if (args$normalisation_method=="deseq2") {
  
  suppressPackageStartupMessages(library(DESeq2))
  dds <- DESeqDataSet(sce_pseudobulk, design=~1)
  dds <- varianceStabilizingTransformation(dds)
  logcounts(sce_pseudobulk) <- assay(dds)
  
} else if (args$normalisation_method=="cpm") {
  
  logcounts(sce_pseudobulk) <- log2(1e6*(sweep(counts(sce_pseudobulk),2,colSums(counts(sce_pseudobulk)),"/"))+1)
  
} else {
  stop("Normalisation method not recognised")
}

sce_pseudobulk$celltype <- colnames(sce_pseudobulk)

# Save
saveRDS(sce_pseudobulk, file.path(args$outdir,"SingleCellExperiment_pseudobulk.rds"))