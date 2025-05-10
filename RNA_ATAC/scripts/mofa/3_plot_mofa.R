###################
##               ##
##  Plot_mofa.R  ##
##               ##
###################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA_ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--mofa_model',        type="character",                               help='MOFA model')
p$add_argument('--outdir',        type="character",                               help='Output directory')
p$add_argument('--samples',        type="character",      nargs="+",       help='Samples')
p$add_argument('--incl_samples',        type="character",                               help='Suffix')
p$add_argument('--n_factors',       type="integer",  help='Amount of factors to use for clustering')
p$add_argument('--sort_samples',  default=TRUE,  help='Sort samples?')
p$add_argument('--batch_correction',  default=FALSE,  help='Correct for batch effect?')
p$add_argument('--seed',        type="integer", default=42,               help='Seed')
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

print("Stat sample metadata")
print(colnames(sample_metadata))
print(dim(sample_metadata))
print(median(sample_metadata$nFeature_RNA))
print(median(sample_metadata$nFrags_atac))

#############################
## Plot variance explained ##
#############################

p <- plot_variance_explained(MOFAobject, plot_total = T)[[2]]

pdf(sprintf("%s/mofa_var_explained_total_%s.pdf",args$outdir,args$incl_samples), width=6, height=3)
print(p)
dev.off()

p <- plot_variance_explained(MOFAobject, factors = 1:MOFAobject@dimensions$K, x="view", y="factor", max_r2 = 5) +
  theme(legend.position = "top")

pdf(sprintf("%s/mofa_var_explained_%s.pdf",args$outdir,args$incl_samples), width=6, height=3)
print(p)
dev.off()

######################
## Batch correction ##
######################

# Select factors to use 
factors.to.use <- 1:get_dimensions(MOFAobject)[["K"]]
factors.to.use <- 1:args$n_factors

# Extract factors
Z <- get_factors(MOFAobject, factors=factors.to.use)[[1]]

print("factors")
print(factors.to.use)
print(dim(Z))
print(head(Z))

Z_corrected <- Z

print("Z corrected")
print(dim(Z_corrected))
print(head(Z_corrected))

# Save factors
factors.dt <- Z_corrected %>% round(3) %>% as.data.table(keep.rownames = T) %>% setnames("rn","cell")
fwrite(factors.dt, sprintf("%s/factors_%s.txt.gz",args$outdir,args$incl_samples))

##########
## UMAP ##
##########

# Run
set.seed(args$seed)
umap_embedding <- uwot::umap(Z_corrected, n_neighbors=200, min_dist=0.05, metric="cosine")

# Plot
to.plot <- umap_embedding %>% as.data.table %>%
  .[,sample:=rownames(Z_corrected)] %>%
  merge(MOFAobject@samples_metadata[,c("sample","batch", "RNA_cluster")] %>% as.data.table)

if (args$sort_samples){
  to.plot$batch <- factor(to.plot$batch,levels=args$samples)
}

print("UMAP")
print(head(to.plot))
fwrite(to.plot, sprintf("%s/umap_%s.txt.gz",args$outdir,args$incl_samples))

print(dim(to.plot))
print(table(to.plot$batch))

p <- ggplot(to.plot, aes(x=V1, y=V2, fill=batch)) +
  geom_point(size=1, shape=21, stroke=0.05) +
  guides(fill = guide_legend(override.aes = list(size=2))) +
  theme_classic() +
  ggplot_theme_NoAxes() 

if (!is.null(opts$color_scheme)){
  p <- p + scale_fill_manual(values=opts$color_scheme[1:length(args$samples)])
}

pdf(sprintf("%s/mofa_umap_sample_%s.pdf",args$outdir,args$incl_samples), width=7, height=5)
print(p)
dev.off()

to.plot$RNA_cluster <- as.factor(to.plot$RNA_cluster)

# Plot UMAP with RNA_cluster as color
p_cluster <- ggplot(to.plot, aes(x = V1, y = V2, fill = RNA_cluster)) +
  geom_point(size = 1, shape = 21, stroke = 0.05) +
  guides(fill = guide_legend(override.aes = list(size = 2))) +
  theme_classic() +
  ggplot_theme_NoAxes()

if (!is.null(opts$color_scheme)) {
  p_cluster <- p_cluster + scale_fill_manual(values = opts$color_scheme[1:length(unique(to.plot$RNA_cluster))])
}

# Save RNA cluster plot
pdf(sprintf("%s/mofa_umap_RNAcluster_%s.pdf", args$outdir, args$incl_samples), width = 7, height = 5)
print(p_cluster)
dev.off()