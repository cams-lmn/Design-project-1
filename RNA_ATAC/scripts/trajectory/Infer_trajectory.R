##########################
##                      ##
##  Infer_trajectory.R  ##
##                      ##
##########################

source("/data/homes/louisc/Project_Babraham/RNA_ATAC/scripts/Settings.R")
library(ggbeeswarm)
library(ggthemes)

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--metadata',  type="character",              help='Metadata') 
p$add_argument('--samples',    type="character",  nargs="+", help='Samples')
p$add_argument('--sce',  type="character",              help='Single cell expression data (metacell)') 
p$add_argument('--sce.atac',  type="character",              help='Single cell accessibility data (metacell)') 
p$add_argument('--trajectory_name',  type="character",              help='Trajectory name') 
p$add_argument('--celltype_label',  type="character",              help='Celltype label') 
p$add_argument('--outdir',          type="character",                help='Output directory')
p$add_argument('--incl_samples',        type="character",                               help='Suffix')
args <- p$parse_args(commandArgs(TRUE))

# ## START TEST ##
# args <- list()
# args$sce <-"/data/homes/louisc/Project_Babraham/RNA_ATAC/metacells/SingleCellExperiment_metacells_nodiff.rds"
# args$sce.atac <-"/data/homes/louisc/Project_Babraham/RNA_ATAC/metacells/PeakMatrix_summarized_experiment_metacells_nodiff.rds"
# # args$metadata <- file.path(io$basedir,"results/atac/archR/celltype_assignment/sample_metadata_after_celltype_assignment.txt.gz")
# args$metadata <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/clustering/PeakMatrix/sample_metadata_nodiff_after_clustering.txt.gz"
# args$trajectory_name <- "N2P"
# args$celltype_label <- "cluster"
# args$outdir <- sprintf("/data/homes/louisc/Project_Babraham/RNA_ATAC/trajectories/%s",args$trajectory_name)
# args$samples <- c("d0","d1","d3","d5","d7","d10","d14","d18","DE","NE","AME_E","AME_L")
# args$filter_differentiated <- TRUE
# ## END TEST ##


#####################
## Define settings ##
#####################

# I/O
dir.create(args$outdir, showWarnings=F,recursive=T)

# incl samples
if(grepl("nodiff",args$incl_samples)){
  args$filter_differentiated <- TRUE
} else {
  args$filter_differentiated <- FALSE
}

# Options
opts$celltype_trajectory_dic <- list(
  "N2P" = c(1,2,3,4,5)
)

if (args$filter_differentiated){
  print("Removing differentiated cells...")
  args$samples <- args$samples[grepl("d[0-9]+",args$samples)]
}

stopifnot(args$trajectory_name%in%names(opts$celltype_trajectory_dic))
opts$celltypes <- opts$celltype_trajectory_dic[[args$trajectory_name]]

opts$genes2plot <- list(
  "N2P" = c("KLF4","NANOG","ETV4")
)

stopifnot(args$trajectory_name%in%names(opts$genes2plot))
opts$genes_to_plot <- opts$genes2plot[[args$trajectory_name]]

##########################
## Load sample metadata ##
##########################

print("Sample metadata")

sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE]

stopifnot(args$celltype_label%in%colnames(sample_metadata))
sample_metadata <- sample_metadata   %>%
  .[,celltype:=eval(as.name(args$celltype_label))] %>%
  .[celltype%in%opts$celltypes] %>%
  .[,celltype:=factor(celltype,levels=opts$celltypes)]

table(sample_metadata$celltype)
colnames(sample_metadata)[1] <- "cell"

# Save
fwrite(sample_metadata, file.path(args$outdir,sprintf("%s_sample_metadata.txt.gz",args$trajectory_name)))

#########################
## Load RNA expression ##
#########################

sce <- readRDS(args$sce)
sce.atac <- readRDS(args$sce.atac)

print(sce)
print(sce.atac)

print(sum(sample_metadata$cell[sample_metadata$cell%in%colnames(sce)]==colnames(sce)))

sce$celltype <- sample_metadata$celltype[sample_metadata$cell%in%colnames(sce)]
sce.atac$celltype <- sample_metadata$celltype[sample_metadata$cell%in%colnames(sce)]

assay(sce.atac,"logcounts") <- assay(sce.atac,"PeakMatrix")

#######################
## Feature selection ##
#######################

decomp <- modelGeneVar(sce)
decomp <- decomp[decomp$mean > 0.01,]
hvgs <- rownames(decomp)[decomp$p.value <= 0.01]

# Subset HVGs
sce_filt <- sce[hvgs,]
dim(sce_filt)

decomp.atac <- modelGeneVar(sce.atac)
decomp.atac <- decomp.atac[decomp.atac$mean > 0.01,]
hvgs.atac <- rownames(decomp.atac)[decomp.atac$p.value <= 0.01]

# Subset HVGs
sce.atac_filt <- sce.atac[hvgs.atac,]
dim(sce.atac_filt)

###############
## PCA : RNA ##
###############

sce_filt <- runPCA(sce_filt, ncomponents = 25, ntop=nrow(sce_filt))

# Define the starting point
foo <- mean(reducedDim(sce_filt,"PCA")[sce_filt$celltype==opts$celltypes[1],"PC1"])
if (foo>0) {
  reducedDim(sce_filt,"PCA")[,"PC1"] <- -reducedDim(sce_filt,"PCA")[,"PC1"]
}

# Plot PCA
p <- plotPCA(sce_filt, colour_by="celltype", ncomponents = c(1,2)) 
output_plot(p, file.path(args$outdir,sprintf("%s_pca_celltype.",args$trajectory_name)), width=8, height=5) 

# ##########
# ## MOFA ##
# ##########

# lsa_atac <- lsa(assay(sce.atac,"logcounts"),dims=25)$dk

# args$python_path <- "/opt/anaconda/Anaconda3/envs/python_Babraham/bin/python"
# reticulate::use_python(args$python_path,require=T)
# library(MOFA2)

# MOFAobject <- create_mofa_from_matrix(list("RNA" = t(reducedDim(sce_filt,"PCA")[,]), "ATAC" = t(lsa_atac)))
# # MOFAobject <- create_mofa_from_matrix(list("RNA" = assay(sce,"logcounts"), "ATAC" = assay(sce.atac,"logcounts")))
# MOFAobject

# data_opts <- get_default_data_options(MOFAobject)
# head(data_opts)

# model_opts <- get_default_model_options(MOFAobject)
# head(model_opts)

# train_opts <- get_default_training_options(MOFAobject)
# head(train_opts)

# MOFAobject <- prepare_mofa(
#   object = MOFAobject,
#   data_options = data_opts,
#   model_options = model_opts,
#   training_options = train_opts
# )

# outfile = "/data/homes/louisc/Project_Babraham/RNA_ATAC/trajectories/N2P/model.hdf5"
# MOFAobject.trained <- run_mofa(MOFAobject, outfile)

# Z <- get_factors(MOFAobject.trained, factors=1:15)[[1]]

# to.plot <- data.frame(MOFA=Z[,1],
#                       dos=factor(gsub("#.+","",rownames(Z)),levels=args$samples))
# colnames(to.plot) <- c("psdt","dos")
# p <- ggplot(to.plot, aes(x = psdt, y = dos, 
#                                colour = dos)) +
#   geom_quasirandom(groupOnX = FALSE) +
#   scale_color_manual() +  theme_classic() +
#   xlab("MOFA") + ylab("Timepoint") +
#   ggtitle(paste0("Cells ordered by MOFA"))
  
# pdf(file.path(args$outdir,sprintf("%s_MOFA_dayofsampling.pdf",args$trajectory_name)), width=6, height=3.5)
# print(p)
# dev.off()

# ################
# ## PCA : ATAC ##
# ################

# sce.atac_filt_pca <- calculatePCA(sce.atac_filt, ncomponents = 5, ntop=nrow(sce.atac_filt))
# sce.atac_filt_pca <- data.frame(sce.atac_filt_pca,stringsAsFactors = F)

# # Plot PCA
# hex <- hue_pal()(max(as.numeric(sce.atac_filt$celltype)))
# pdf(file.path(args$outdir,sprintf("%s_atac_pca_celltype.pdf",args$trajectory_name)), width=8, height=5) 
# plot(sce.atac_filt_pca$PC1,sce.atac_filt_pca$PC2, pch=18, col=hex[as.numeric(sce.atac_filt$celltype)]) 
# dev.off()

# #################
# ## PCA : Combo ##
# #################

# # Plot PCA
# hex <- hue_pal()(max(as.numeric(sce.atac_filt$celltype)))
# pdf(file.path(args$outdir,sprintf("%s_combo_pca_celltype.pdf",args$trajectory_name)), width=8, height=5) 
# plot(reducedDim(sce_filt,"PCA")[,"PC1"],sce.atac_filt_pca$PC1, pch=16, col=hex[as.numeric(sce.atac_filt$celltype)]) 
# legend("bottomright",as.character(1:5),pch=16,col=hex)
# dev.off()

# ################
# ## Combine PCA #
# ################

# summary(reducedDim(sce_filt,"PCA"))
# summary(sce.atac_filt_pca)

# PCA_combo <- scale(reducedDim(sce_filt,"PCA")[,"PC1"]) + scale(sce.atac_filt_pca$PC1)
# # PCA_combo <- reducedDim(sce_filt,"PCA")[,"PC1"] + c(sce.atac_filt_pca$PC1/10)

###################
## Diffusion map ##
###################

set.seed(42)
dm <- DiffusionMap(sce_filt, n_pcs=5)

reducedDim(sce_filt, "DiffusionMap") <- dm@eigenvectors[,1:5]

# Define the starting point
foo <- mean(reducedDim(sce_filt,"DiffusionMap")[sce_filt$celltype==opts$celltypes[1],"DC1"])
if (foo>0) {
  reducedDim(sce_filt,"DiffusionMap")[,"DC1"] <- -reducedDim(sce_filt,"DiffusionMap")[,"DC1"]
}

# Plot
p <- plotReducedDim(sce_filt, dimred = "DiffusionMap", colour_by="celltype", ncomponents = c(1,2)) 
output_plot(p, file.path(args$outdir,sprintf("%s_diffmap_celltype",args$trajectory_name)), width=8, height=5) 

####################
## Prepare output ##
####################

print("Creation pseudotime object")

pseudotime.dt <- data.table(
  cell = colnames(sce_filt),
  PC1 = reducedDim(sce_filt,"PCA")[,"PC1"],
  DC1 = reducedDim(sce_filt,"DiffusionMap")[,"DC1"],
  PC2 = reducedDim(sce_filt,"PCA")[,"PC2"],
  DC2 = reducedDim(sce_filt,"DiffusionMap")[,"DC2"],
  PC3 = reducedDim(sce_filt,"PCA")[,"PC3"],
  DC3 = reducedDim(sce_filt,"DiffusionMap")[,"DC3"],
  dos = factor(gsub("#.+$","",colnames(sce_filt)),levels=args$samples),
  cluster = as.factor(sample_metadata$cluster[sample_metadata$cell%in%colnames(sce_filt)])
) %>% setorder(PC1)

print(head(pseudotime.dt))

# Save
fwrite(pseudotime.dt, file.path(args$outdir,sprintf("/%s_trajectory_%s.txt.gz",args$trajectory_name,args$incl_samples)))

# Load precomputed pseudotime
# pseudotime.dt <- fread(paste0(args$outdir,"/blood_trajectory.txt.gz"))

################################################
## Boxplots of pseudotime values per celltype ## 
################################################

to.plot <- pseudotime.dt %>%
  merge(sample_metadata[,c("cell","celltype")], by="cell") %>%
  .[,pca:=minmax.normalisation(PC1)] %>% .[,diffmap:=minmax.normalisation(DC1)]

for (i in c("pca","diffmap")) {
  p <- ggboxplot(to.plot, x="celltype", y=i, fill="celltype", outlier.shape=NA) +
    stat_summary(fun.data = give.n, geom = "text", size=4) +
    labs(x="", y=sprintf("Pseudotime (%s)",i)) +
    theme_classic() 
  
  output_plot(p, file.path(args$outdir,sprintf("%s_%s_boxplot_celltype",args$trajectory_name,i)), width=6, height=3.5)
}

####################################################################
## pseudotime based on combo of PC1-3 & DC1-3 per day of sampling ## 
####################################################################

for (i in colnames(pseudotime.dt)[grepl("[DP]C",colnames(pseudotime.dt))]) {
  p <- ggplot(pseudotime.dt, aes(x = get(i), y = dos, 
                                colour = dos)) +
    geom_quasirandom(groupOnX = FALSE) +
    scale_color_manual(values=opts$color_scheme) +  theme_classic() +
    xlab(i) + ylab("Timepoint")

  output_plot(p, file.path(args$outdir,paste0(args$trajectory_name,"_",i,"_dayofsampling")), width=6, height=3.5)
}

# ######################################################
# ## pseudotime based on PCAcombo per day of sampling ## 
# ######################################################

# to.plot <- data.frame(PCA_combo=PCA_combo[,1],
#                       dos=factor(gsub("#.+","",rownames(PCA_combo)),levels=args$samples))
# colnames(to.plot) <- c("psdt","dos")
# p <- ggplot(to.plot, aes(x = psdt, y = dos, 
#                                colour = dos)) +
#   geom_quasirandom(groupOnX = FALSE) +
#   scale_color_manual() +  theme_classic() +
#   xlab("PCA combo") + ylab("Timepoint") +
#   ggtitle(paste0("Cells ordered by PCA combo"))
  
# pdf(file.path(args$outdir,sprintf("%s_PCAcombo_dayofsampling.pdf",args$trajectory_name)), width=6, height=3.5)
# print(p)
# dev.off()

# pseudotime.dt <- data.table(
#   cell = colnames(sce_filt),
#   PCA_combo = PCA_combo[,1],
#   dos = factor(gsub("#.+","",rownames(PCA_combo)),levels=args$samples)
# ) %>% setorder(PCA_combo)

# # Save
# fwrite(pseudotime.dt, file.path(args$outdir,sprintf("/%s_trajectory.txt.gz",args$trajectory_name)))

################################################################
## pseudotime based on combo of PC1 & PC2 per day of sampling ## 
################################################################

p <- ggplot(pseudotime.dt, aes(x = PC1+PC2/5, y = dos, 
                               colour = dos)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values=opts$color_scheme) +  theme_classic() +
  xlab("PC1+PC2/5") + ylab("Timepoint")

output_plot(p, file.path(args$outdir,paste0(args$trajectory_name,"_PC1-2combo_diffmap_dayofsampling")), width=6, height=3.5)

################################################################
## pseudotime based on combo of PC1 & DC1 per day of sampling ## 
################################################################

p <- ggplot(pseudotime.dt, aes(x = PC1+DC1*100, y = dos, 
                               colour = dos)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values=opts$color_scheme) +  theme_classic() +
  xlab("PC1+DC1*100") + ylab("Timepoint") 

output_plot(p, file.path(args$outdir,paste0(args$trajectory_name,"_PC1-DC1combo_diffmap_dayofsampling")), width=6, height=3.5)

p <- ggplot(pseudotime.dt, aes(x = PC1+DC1*100, y = cluster, 
                               colour = cluster)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values=opts$celltype.colors[1:length(levels(pseudotime.dt$cluster))]) +  theme_classic() +
  xlab("PC1+DC1*100") + ylab("Timepoint") 

output_plot(p, file.path(args$outdir,paste0(args$trajectory_name,"_PC1-DC1combo_diffmap_cluster")), width=6, height=2.5)


#########################################################################
## Scatterplots of pseudotime values versus expression of marker genes ## 
#########################################################################

# Denoise
# pca.rna <- fread(io$pca.rna) %>% matrix.please %>% .[colnames(sce),] # assay(sce_filt,"logcounts_denoised") <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(sce_filt)), D=pdist(pca.rna), k=25)

rna.dt <- data.table(as.matrix(logcounts(sce)[opts$genes_to_plot,]), keep.rownames = T) %>%
  setnames("rn","gene") %>%
  melt(id.vars="gene", variable.name="cell", value.name="expr") %>%
  merge(sample_metadata[,c("cell")])

to.plot2 <- to.plot %>%
  merge(rna.dt[gene%in%opts$genes_to_plot],by="cell") %>%
  .[,gene:=factor(gene,levels=opts$genes_to_plot)] %>%
  .[,expr:=minmax.normalisation(expr),by="gene"]

for (i in c("pca","diffmap")) {
  p <- ggplot(to.plot2, aes_string(x=i, y="expr")) +
    geom_point(aes(fill=celltype), size=1.75, shape=21, stroke=0.1, alpha=0.5) +
    # ggrastr::geom_point_rast(aes(fill=celltype), size=1.5, shape=21, stroke=0.1, alpha=0.75) +
    #stat_smooth(aes(fill=expr), method="loess", color="black", alpha=0.85, span=0.5) +
    geom_rug(aes(color=celltype), sides="b") +
    # viridis::scale_fill_viridis() +
    # scale_color_manual(values=opts$celltype.colors) +
    # scale_fill_manual(values=opts$celltype.colors) +
    facet_wrap(~gene, scales="free_y") +
    labs(x=sprintf("Pseudotime (%s)",i), y="Gene expression") +
    theme_classic() +
    guides(fill="none", color="none") +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      # axis.text.y = element_text(color="black"),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text = element_text(size=rel(1.25), color="black"),
      legend.position="right",
      legend.title = element_blank()
    )
  
  output_plot(p, file.path(args$outdir,sprintf("%s_%s_vs_expression",args$trajectory_name,i)), width=7.5, height=4)
}


###############################
## Save SingleCellExperiment ##
###############################

logcounts(sce) <- NULL
saveRDS(sce, file.path(args$outdir,sprintf("%s_SingleCellExperiment.rds",args$trajectory_name)))

##################
## Save Anndata ##
##################

# # sce <- readRDS(file.path(args$outdir,sprintf("%s_SingleCellExperiment.rds",args$trajectory_name)))
#
# # reticulate::use_condaenv("celloracle_env", required = T) # reticulate::use_condaenv("scanpy", required = T) # sc <- import("scanpy") # # # Create anndata # adata_sce <- sc$AnnData(
#     X   = t(counts(sce)),
#     obs = as.data.frame(colData(sce)),
#     var = data.frame(gene=rownames(sce), row.names=rownames(sce))
# )
#
# # Add cell type colors
# adata_sce$uns$update("celltype_colors" = opts$celltype.colors[sort(unique(as.character(adata_sce$obs[[args$celltype_label]])))])
#
# # Save
# adata_sce$write_h5ad(file.path(args$outdir,sprintf("%s_adata.h5ad",args$trajectory_name)))
