#######################
##                   ##
##  Seurat_to_SCE.R  ##
##                   ##
#######################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--normalise',        type="logical", default=TRUE,                   help='Log-Normalise?')
p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
p$add_argument('--seurat',         type="character", help='Seurat object (input)')
p$add_argument('--metadata',         type="character", help='Metadata file')
p$add_argument('--metadata_out',         type="character", help='Metadata output file')
p$add_argument('--outfile',         type="character", help='Output file')
p$add_argument('--sort_samples',       type="logical",  default=TRUE,       help='Sort samples?')
args <- p$parse_args(commandArgs(TRUE))

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
sample_metadata <- sample_metadata %>% .[cell%in%colnames(sce)] %>% setkey(cell) %>% .[colnames(sce)]
stopifnot(sample_metadata$cell == colnames(sce))
colData(sce) <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(sce),] %>% DataFrame()

##########################
## Compute size factors ## problem non positive factors 
##########################

clusts <- as.numeric(quickCluster(sce, method = "igraph", min.size = 200, BPPARAM = mcparam))
min.clust <- min(table(clusts))/2
new_sizes <- c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
sce <- computeSumFactors(sce, clusters = clusts, sizes = new_sizes, max.cluster.size = 3000, positive = TRUE) # positive = TRUE for negative values

###################
## Log Normalise ##
###################

if (args$normalise) {
  sce <- logNormCounts(sce)
}

####################
## Add cell state ##
####################

#convertHumanGeneList <- function(human_genes) {
#    human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://asia.ensembl.org")
#    mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://asia.ensembl.org")
#    genesV2 <- getLDS(attributes = c("hgnc_symbol"),
#                      filters = "hgnc_symbol",
#                      values = human_genes,
#                      mart = human,
#                      attributesL = c("mgi_symbol"),
#                      martL = mouse,
#                      uniqueRows = TRUE)
#    mouse_genes <- unique(genesV2[, 2])
#    return(mouse_genes)
#}

#s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
#g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)

#seurat <- CellCycleScoring(seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#sample_metadata$Phase <- seurat[[]]$Phase

fwrite(sample_metadata,args$metadata_out, quote=F, na="NA", sep="\t")

##########
## Save ##
##########

saveRDS(sce, args$outfile)