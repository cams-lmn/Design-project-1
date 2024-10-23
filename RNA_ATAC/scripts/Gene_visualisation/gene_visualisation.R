############################
##                        ##
##  Gene_visualisation.R  ##
##                        ##
############################

source("/data/homes/louisc/Project_Babraham/RNA_ATAC/scripts/Settings.R")

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

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args <- list()
# args$filter_differentiated <- TRUE
# if (args$filter_differentiated){
#   args$metadata <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/mofa/fast/sample_metadata_nodiff.txt.gz"
#   args$mofa_model <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/mofa/fast/mofa_nodiff.rds"
# } else {
#   args$metadata <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/mofa/fast/sample_metadata.txt.gz"
#   args$mofa_model <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/mofa/fast/mofa.rds"
# }
# args$batch_correction <- FALSE
# args$outdir <-"/data/homes/louisc/Project_Babraham/RNA_ATAC/mofa/fast/pdf"
# args$samples <- c("d0","d1","d3","d5","d7","d10","d14","d18","DE","NE","AME_E","AME_L")
# args$sort_samples <- TRUE
# args$color_scheme <- c("#A06932","#DBB216","#EFE32A","#D7DB54","#A7B019","#9AA126",
#                        "#7E7721","#5D6821","#353D8A","#CA3639","#DCADCD","#7A327E")
# args$seed <- 42
## END TEST ##


##########
## UMAP ##
##########

# Run
set.seed(args$seed)
umap_embedding <- uwot::umap(Z_corrected, n_neighbors=25, min_dist=0.50, metric="cosine")

# Plot
to.plot <- umap_embedding %>% as.data.table %>%
  .[,sample:=rownames(Z_corrected)] %>%
  merge(MOFAobject@samples_metadata[,c("sample","batch")] %>% as.data.table)

if (args$sort_samples){
  to.plot$batch <- factor(to.plot$batch,levels=args$samples)
}

print(dim(to.plot))
print(table(to.plot$batch))

p <- ggplot(to.plot, aes(x=V1, y=V2, fill=batch)) +
  geom_point(size=1, shape=21, stroke=0.05) +
  # ggrastr::geom_point_rast(size=1, shape=21, stroke=0.1) +
  guides(fill = guide_legend(override.aes = list(size=2))) +
  theme_classic() +
  # theme(legend.position="none") +
  ggplot_theme_NoAxes() 

if (!is.null(opts$color_scheme)){
  p <- p + scale_fill_manual(values=opts$color_scheme[1:length(args$samples)])
}

pdf(sprintf("%s/mofa_umap_sample_%s.pdf",args$outdir,args$incl_samples), width=7, height=5)
print(p)
dev.off()
