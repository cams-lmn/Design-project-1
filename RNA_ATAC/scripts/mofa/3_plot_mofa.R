###################
##               ##
##  Plot_mofa.R  ##
##               ##
###################

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

#####################
## Define settings ##
#####################

# I/O
dir.create(args$outdir, showWarnings = F, recursive = T)

if(grepl("nodiff",args$incl_samples)){
  args$filter_differentiated <- TRUE
} else {
  args$filter_differentiated <- FALSE
}

if (args$filter_differentiated){
  args$samples <- args$samples[grepl("d[0-9]+",args$samples)]
}

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) 
print(head(sample_metadata))

if (args$batch_correction){
  sample_metadata$SeqRun <- "Run2"
  sample_metadata$SeqRun[sample_metadata$batch=="d0"] <- "Jasmin"
  sample_metadata$SeqRun[sample_metadata$batch%in%c("d10","d14","d18","DE")] <- "Run1"
}

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

#################################
## Correlation between factors ##
#################################

# plot_factor_cor(MOFAobject)

####################
## Subset factors ##
####################

# r2 <- MOFAobject@cache$variance_explained$r2_per_factor
# factors <- sapply(r2, function(x) x[,"RNA"]>0.01)
# MOFAobject <- subset_factors(MOFAobject, which(apply(factors,1,sum)>=1))
# factors(MOFAobject) <- paste("Factor",1:get_dimensions(MOFAobject)[["K"]], sep=" ")

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

##################
## Plot factors ##
##################

# plot_factor(MOFAobject, factors = c(1), color_by = "celltype", group_by = "celltype", add_violin = T, add_boxplot = T, add_dots = F, dodge=T, legend=F) +
#   scale_fill_manual(values=opts$celltype.colors) +
#   theme(
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank()
#   )
# 
# 
# plot_factor(MOFAobject, factors = 3, color_by = "celltype", group_by = "celltype", 
#             add_boxplot = T, add_dots = F, add_violin = T, dodge=T, legend=F) +
#   scale_fill_manual(values=opts$celltype.colors) +
#   theme(
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank()
#   )
# 
# p <- plot_factors(MOFAobject, factors = c(1,2), color_by = "celltype", dot_size = 1) +
#   scale_fill_manual(values=opts$celltype.colors) +
#   theme(
#     legend.position = "none"
#   )

######################
## Sumarise factors ##
######################

# levels_df <- MOFAobject@samples_metadata[,c("sample","celltype")] %>% setnames("celltype","level")
# p <- summarise_factors(MOFAobject, levels_df, factors = 1:25, abs = F, return_data = F) +
#   guides(x = guide_axis(angle = 90)) +
#   theme(
#     legend.position = "top",
#     axis.text.x = element_text(color="black", size=rel(0.75)),
#     axis.text.y = element_text(color="black", size=rel(0.75)),
#   )
# 
# pdf(sprintf("%s/summarise_factor_celltype.pdf",args$outdir), width=7, height=5)
# print(p)
# dev.off()

##################
## Plot weights ##
##################

# plot_weights(MOFAobject, factor = 1, view="ATAC", nfeatures = 10, text_size = 3)
# plot_weights(MOFAobject, factor = 1, view="RNA", nfeatures = 15, text_size = 4)

#######################################
## Correlate factors with covariates ##
#######################################

# foo <- correlate_factors_with_covariates(MOFAobject, covariates = c("nFeature_RNA","nFrags_atac","ribosomal_percent_RNA","mitochondrial_percent_RNA"), return_data = T)
# pheatmap::pheatmap(foo)

######################
## Batch correction ##
######################

# Select factors to use 
factors.to.use <- 1:get_dimensions(MOFAobject)[["K"]]
#factors.to.use <- 1:args$n_factors
# factors.to.use <- factors.to.use[!factors.to.use%in%c("3")]
# factors.to.use <- c(1,2,3,6,7)

# Extract factors
Z <- get_factors(MOFAobject, factors=factors.to.use)[[1]]

print("factors")
print(factors.to.use)
print(dim(Z))
print(head(Z))

if (args$batch_correction) {
  # Define stage and sample order
  timepoints <- MOFAobject@samples_metadata$SeqRun
  timepoint_order <- unique(MOFAobject@samples_metadata$SeqRun)
  samples <- MOFAobject@samples_metadata$batch
  sample_order <- args$samples
  
  Z_list    <- lapply(unique(timepoints), function(i){
    sub_pc   <- Z[timepoints == i, , drop = FALSE]
    sub_samp <- samples[timepoints == i]
    list     <- lapply(unique(sub_samp), function(j){ sub_pc[sub_samp == j, , drop = FALSE]})
    names(list) <- unique(sub_samp)
    return(list)
  })
  names(Z_list) <- unique(timepoints)
  
  # #arrange to match timepoint order
  Z_list <- Z_list[order(match(names(Z_list), timepoint_order))]
  Z_list <- lapply(Z_list, function(x){ x[order(match(names(x), sample_order))]})
  
  #perform corrections within stages
  correct_list <- lapply(Z_list, function(x){
    if(length(x) > 1){
      return(do.call(reducedMNN, x)$corrected)
    } else {
      return(x[[1]])
    }
  })
  
  #perform correction over stages
  Z_corrected <- reducedMNN(correct_list, merge.order=1:length(correct_list))$corrected 
  # colnames(Z_corrected) <- colnames(Z)
} else {
  Z_corrected <- Z
}

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
umap_embedding <- uwot::umap(Z_corrected, n_neighbors=25, min_dist=0.50, metric="cosine")

# Plot
to.plot <- umap_embedding %>% as.data.table %>%
  .[,sample:=rownames(Z_corrected)] %>%
  merge(MOFAobject@samples_metadata[,c("sample","batch")] %>% as.data.table)

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

# p <- ggplot(to.plot[sample(.N,.N/3)], aes(x=V1, y=V2, fill=stage)) +
#   ggrastr::geom_point_rast(size=1, shape=21, stroke=0.1) +
#   scale_fill_manual(values=opts$stage.colors) +
#   guides(fill = guide_legend(override.aes = list(size=2))) +
#   theme_classic() +
#   ggplot_theme_NoAxes() +
#   theme(
#     legend.position="top",
#     legend.title = element_blank()
#   )
# 
# pdf(sprintf("%s/mofa_umap_stage.pdf",args$outdir), width=7, height=5)
# print(p)
# dev.off()

# plot_dimred(MOFAobject, method="UMAP", color_by="stage")

#########################
## Contribution scores ##
#########################

# factors.to.use <- "all"
# 
# r2.per.sample <- calculate_variance_explained_per_sample(MOFAobject)[[1]]
# # foo <- rowSums(r2.per.sample)
# # r2.per.sample["E8.0_rep1#GCGAAGTAGTACCGCA-1",]
# # foo["E8.0_rep1#GCGAAGTAGTACCGCA-1"]
# # sort(foo) %>% head
# MOFAobject <- calculate_contribution_scores(MOFAobject, factors = factors.to.use, scale = TRUE)
# 
# cells <- intersect(names(which(r2.per.sample[,"RNA"]>=3)), names(which(r2.per.sample[,"ATAC"]>=3)))
# 
# to.plot <- MOFAobject@samples_metadata[MOFAobject@samples_metadata$sample%in%cells,c("sample","RNA_contribution","ATAC_contribution","celltype")] %>% as.data.table
# 
# order.celltypes <- to.plot[,median(RNA_contribution),by="celltype"] %>% setorder(-V1) %>% .$celltype
# to.plot[,celltype:=factor(celltype, levels=order.celltypes)]
# 
# 
# p <- ggplot(to.plot, aes(x=celltype, y=RNA_contribution)) +
#   geom_boxplot(aes(fill = celltype), alpha=0.9, outlier.shape=NA, coef=1.5) +
#   coord_flip(ylim = c(0.10,0.90)) +
#   geom_hline(yintercept=0.5, linetype="dashed", size=0.5) +
#   scale_fill_manual(values=opts$celltype.colors, drop=F) +
#   theme_classic() +
#   labs(y="RNA contribution score", x="") +
#   theme(
#     legend.position = "none",
#     axis.title.y = element_blank(),
#     axis.text.y = element_text(color="black"),
#     axis.text.x = element_text(color="black")
#   )
# 
# 
# pdf(sprintf("%s/contribution_scores.pdf",args$outdir), width=5, height=7)
# print(p)
# dev.off()


##################################################################
## Plot cumulative variance explained per view vs factor number ##
##################################################################

# factors <- 1:25
# 
# r2.dt <- MOFAobject@cache$variance_explained$r2_per_factor[[1]][factors,] %>%
#   as.data.table %>% .[,factor:=as.factor(factors)] %>%
#   melt(id.vars="factor", variable.name="view", value.name = "r2") %>%
#   .[,cum_r2:=cumsum(r2), by="view"]
# 
# # threshold.var <- 5
# # max.factor <- max(which(apply(r2,1,sum) >= threshold.var))
# 
# p <- ggline(r2.dt, x="factor", y="cum_r2", color="view") +
#   # scale_color_manual(values=opts$colors.views) +
#   labs(x="Factor number", y="Cumulative variance explained (%)") +
#   # geom_vline(xintercept = max.factor, linetype="dashed") +
#   theme(
#     legend.title = element_blank(), 
#     legend.position = "top",
#     axis.text = element_text(size=rel(0.8))
#   )
# 
# pdf(paste0(args$outdir,"/r2_vs_factor.pdf"), width=8, height=5, useDingbats = F)
# print(p)
# dev.off()