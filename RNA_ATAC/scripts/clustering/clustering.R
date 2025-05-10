####################
##                ##
##  clustering.R  ##
##                ##
####################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA_ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--mofa_model',        type="character",                               help='MOFA model')
p$add_argument('--matrix',          type="character",  default="PeakMatrix",   help='Matrix to use')
p$add_argument('--samples',         type="character",                nargs='+',     help='Samples')
p$add_argument('--sort_samples',  default=TRUE,  help='Sort samples?')
p$add_argument('--incl_samples',        type="character",                               help='Suffix')
p$add_argument('--batch_correction',  default=FALSE,  help='Correct for batch effect?')
p$add_argument('--n_factors',       type="integer",  help='Amount of factors to use for clustering')
p$add_argument('--am_clusters',       type="integer",  help='Amount of clusters to create')
p$add_argument('--reorder_clusters',       type="integer",  default= NULL, nargs='+',  help='New order for clusters')
p$add_argument('--seed1',            type="integer",    default=42,                  help='Random seed 1')
p$add_argument('--seed2',            type="integer",    default=4242,                  help='Random seed 2')
p$add_argument('--outdir',          type="character",                               help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

print(args)

#####################
## Define settings ##
#####################

# I/O
dir.create(args$outdir, showWarnings = F, recursive = T)

args$filter_differentiated <- FALSE

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) 
print(head(sample_metadata))

###############
## Load MOFA ##
###############

if (tools::file_ext(args$mofa_model)=="rds") {
  MOFAobject <- readRDS(args$mofa_model)
  
  # Add sample metadata
  samples_metadata(MOFAobject) <- sample_metadata
} else if (file_ext(args$mofa_model)=="hdf5") {
  MOFAobject <- load_model(args$mofa_model, load_data = F)
  
  # Add sample metadata
  cells <- as.character(unname(unlist(MOFA2::samples_names(MOFAobject))))
  sample_metadata_to_mofa <- copy(sample_metadata) %>%
    setnames("cell","sample") %>%
    .[sample%in%cells] %>% setkey(sample) %>% .[cells]
  stopifnot(all(cells==sample_metadata_to_mofa$cell))
  samples_metadata(MOFAobject) <- sample_metadata_to_mofa
}
print(MOFAobject@dimensions$K)  # Number of trained factors
##################
## Load factors ##
##################

# Select factors to use 
factors.to.use <- 1:args$n_factors

# Extract factors
Z <- get_factors(MOFAobject, factors=factors.to.use)[[1]]

print(factors.to.use)
print(head(Z))

Z_corrected <- Z

print(head(Z_corrected))

# Save factors
factors.dt <- Z_corrected %>% round(3) %>% as.data.table(keep.rownames = T) %>% setnames("rn","cell")
fwrite(factors.dt, sprintf("%s/factors_%s.txt.gz",args$outdir,args$incl_samples))

##########
## UMAP ##
##########

# Run
set.seed(args$seed1)
umap_embedding <- uwot::umap(Z_corrected, n_neighbors=200, min_dist=0.05, metric="cosine")

################
## Clustering ##
################

set.seed(args$seed2)
MOFA_clusters <- cluster_samples(MOFAobject, k = args$am_clusters, factors = 1:args$n_factors)

if (!is.null(args$reorder_clusters)) {
  names_cluster <- names(MOFA_clusters$cluster)
  
  # Define the mapping
  cluster_mapping <- c(
    "1" = "4", "3" = "4", "9" = "4",
    "2" = "3", "7" = "3", "12" = "3",
    "4" = "2",
    "5" = "1", "8" = "1", "11" = "1",
    "6" = "5",
    "10" = "6"
  )
  
  # Convert clusters using the mapping
  MOFA_clusters$cluster <- as.character(cluster_mapping[as.character(MOFA_clusters$cluster)])
  
  # Preserve sample names
  names(MOFA_clusters$cluster) <- names_cluster
}

##########
## Plot ##
##########

# Plot
to.plot <- umap_embedding %>% as.data.table %>%
  .[,sample:=rownames(Z_corrected)] %>%
  merge(MOFAobject@samples_metadata[,c("sample","batch")] %>% as.data.table)

sum(names(MOFA_clusters$cluster)==to.plot$sample)
dim(to.plot)

to.plot$cluster <- as.factor(MOFA_clusters$cluster)

if (args$sort_samples){
  to.plot$batch <- factor(to.plot$batch,levels=args$samples)
}

print(dim(to.plot))
print(table(to.plot$cluster))

print("UMAP")
print(head(to.plot))
fwrite(to.plot, sprintf("%s/umap_clusters_%s.txt.gz",args$outdir,args$incl_samples))

legend.labels <- c(
  "1: Myeloid",
  "2: B cells",
  "3: T cells",
  "4: Fibroblasts",
  "5: Endothelial",
  "6: Tumor cells"
)

p <- ggplot(to.plot, aes(x=V1, y=V2, fill=cluster)) +
  geom_point(size=1, shape=21, stroke=0.1) +
  scale_fill_manual(
    name = "Cell Type",
    values = opts$celltype.colors,
    labels = c("1: Myeloid", "2: B cells", "3: T cells", "4: Fibroblasts", "5: Endothelial", "6: Tumor cells")
    ) +
  guides(fill = guide_legend(override.aes = list(size=5))) +
  theme_classic() +
  ggplot_theme_NoAxes() 

pdf(sprintf("%s/mofa_umap_cluster_%s.pdf",args$outdir,args$incl_samples), width=7, height=5)
print(p)
dev.off()

###################
## Save metadata ##
###################

sum(names(MOFA_clusters$cluster)==sample_metadata$sample)
dim(sample_metadata)

sample_metadata$cluster <- MOFA_clusters$cluster

# Save
fwrite(sample_metadata, sprintf("%s/sample_metadata_%s_after_clustering.txt.gz",args$outdir,args$incl_samples), 
       quote=F, na="NA", sep="\t")
