##########################################
##                                      ##
##  Build_GRN_metacells_trajectories.R  ##
##                                      ##
##########################################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA_ATAC/scripts/Settings.R")

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--trajectory_name',  type="character",              help='Trajectory name')
p$add_argument('--sce',       type="character",                help='SingleCellExperiment')
p$add_argument('--tf2gene_virtual_chip',       type="character",                help='Links between tfs and genes based on virtual chipseq')
p$add_argument('--max_distance',  type="integer",            default=1e5,      help='Maximum distance for a linkage between a peak and a gene')
p$add_argument('--ncores',  type="integer",            default=4,      help='Amount of cores to use')
p$add_argument('--min_chip_score',  type="double",            default=0.25,      help='Minimum in silico ChIP-seq score')
p$add_argument('--markers_TF',  type="character",      help='TF marker file')
p$add_argument('--plot_correlations',  type="logical", default=TRUE, help='Do you want to plot correlations plots')
p$add_argument('--outdir',       type="character",                help='Output directory')
p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
args <- p$parse_args(commandArgs(TRUE))

print(args)

##################
## GRN_function ##
##################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA_ATAC/scripts/GRN/build_GRN_functions.R")

#####################
## Define settings ##
#####################

if(grepl("full",args$outdir)){
  args$plot_correlations <- FALSE
}

if(grepl("supercluster",args$trajectory_name)){
  args$plot_correlations <- FALSE
}

dir.create(args$outdir, showWarnings = F, recursive = T)
dir.create(file.path(dirname(args$tf2gene_virtual_chip),sprintf("plots_cor_score%s",args$min_chip_score)))

myfun_used <- myfun_conv

if (args$plot_correlations){
  dir.create(sprintf("%s/plots_cor_score%s",args$outdir,args$min_chip_score))
}

####################################################
## Load TF2gene links based on in silico ChIP-seq ##
####################################################

tf2gene_chip.dt <- fread(args$tf2gene_virtual_chip) %>%
  .[chip_score>=args$min_chip_score & dist<=args$max_distance] %>% 
  .[,c("tf","gene")] %>% unique # Only keep TF-gene links

print(sprintf("Number of TFs: %s",length(unique(tf2gene_chip.dt$tf))))
print(sprintf("Number of genes: %s",length(unique(tf2gene_chip.dt$gene))))

print(head(tf2gene_chip.dt))

##############################
## Load RNA expression data ##
##############################

sce_mc <- readRDS(args$sce)

if(grepl("supercluster",args$trajectory_name)){
  metadata <- fread('/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA_ATAC/clustering/PeakMatrix/sample_metadata_all_after_clustering.txt.gz')
  metadata <- metadata[metadata$sample %in% colnames(sce_mc),]
  metadata <- metadata[match(colnames(sce_mc), metadata$sample),]

  if (grepl("A",args$trajectory_name)){
    clusters <- c(1, 2, 3)
  } else if (grepl("B",args$trajectory_name)) {
    clusters <- c(4, 5)
  } else if (grepl("C",args$trajectory_name)) {
    clusters <- 6
  }

  selected_cells <- metadata$sample[metadata$cluster %in% clusters]
  sce_mc <- sce_mc[, colnames(sce_mc) %in% selected_cells]
} else {
  clusters <- c(1, 2, 3, 4, 5, 6)
}

##########################
## Filter TFs and genes ##
##########################

TFs <- intersect(unique(tf2gene_chip.dt$tf),rownames(sce_mc))
print(head(TFs))
genes <- intersect(unique(tf2gene_chip.dt$gene),rownames(sce_mc))

if (grepl("full",args$outdir)){
  tf2gene_chip.dt <- tf2gene_chip.dt[tf%in%TFs & gene%in%genes,]
} else {
  tf2gene_chip.dt <- tf2gene_chip.dt[tf%in%TFs & gene%in%TFs,]
}

# Fetch RNA expression matrices
rna_tf.mtx <- logcounts(sce_mc)[unique(tf2gene_chip.dt$tf),]
rownames(rna_tf.mtx) <- rownames(rna_tf.mtx)
rna_targets.mtx <- logcounts(sce_mc)[unique(tf2gene_chip.dt$gene),]

print(dim(rna_tf.mtx))
print(dim(rna_targets.mtx))
print(head(rna_targets.mtx[,1:6]))

# Filter out lowly variable genes and TFs
print(dim(rna_tf.mtx))
print(head(rna_tf.mtx[,1:6]))

TFs <- intersect(unique(tf2gene_chip.dt$tf),rownames(rna_tf.mtx))
genes <- intersect(unique(tf2gene_chip.dt$gene),rownames(rna_targets.mtx))
tf2gene_chip.dt <- tf2gene_chip.dt[tf%in%TFs & gene%in%genes,]

print("Filtering for marker TF")
TF_markers <- read.table(args$markers_TF,header=T,sep="\t")$Gene
print(sum(TFs%in%TF_markers))
print(dim(tf2gene_chip.dt))
TFs <- intersect(unique(tf2gene_chip.dt$tf),rownames(rna_tf.mtx))
genes <- intersect(unique(tf2gene_chip.dt$gene),rownames(rna_targets.mtx))

print(TFs)
print(length(TFs))
print(sprintf("Number of TFs: %s",length(TFs)))
print(sprintf("Number of genes: %s",length(genes)))

####################
## run regression ##
####################

tmp <- tf2gene_chip.dt[tf2gene_chip.dt$gene%in%genes]
print(dim(tmp))
list_results <- mclapply(genes,myfun_conv, cutoff = 0.84, mc.cores=args$ncores,mc.preschedule=T)
GRN_coef.dt <- do.call(rbind,list_results)

system(sprintf("cat %s/x_y_*.txt > %s/x_y.txt", args$outdir, args$outdir))

# delete the temp files
file.remove(list.files(args$outdir, pattern="^x_y_\\d+\\.txt$", full.names=TRUE))

##############
## Save GRN ##
##############

print(head(GRN_coef.dt,20))
print(dim(GRN_coef.dt))
print(sum(unique(GRN_coef.dt$tf)%in%TF_markers))
print(sum(unique(GRN_coef.dt$gene)%in%TF_markers))
print(sum((GRN_coef.dt$tf%in%TF_markers) & (GRN_coef.dt$gene%in%TF_markers)))
print(sum((GRN_coef.dt$tf%in%TF_markers) & (GRN_coef.dt$gene%in%TF_markers) & 
            (abs(as.numeric(GRN_coef.dt$beta))>=0.25)))
print(sum((GRN_coef.dt$tf%in%TF_markers) & (GRN_coef.dt$gene%in%TF_markers) & 
            (abs(as.numeric(GRN_coef.dt$beta))>=0.25) & (as.numeric(GRN_coef.dt$pvalue)<0.10)))
print(sum((GRN_coef.dt$tf%in%TF_markers) & (GRN_coef.dt$gene%in%TF_markers) & 
            (abs(as.numeric(GRN_coef.dt$beta))>=0.25) & (as.numeric(GRN_coef.dt$pvalue)<0.10) &
            (GRN_coef.dt$tf!=GRN_coef.dt$gene)))

GRN_coef.dt <- data.table(GRN_coef.dt)

fwrite(GRN_coef.dt, file.path(args$outdir,sprintf('global_chip_GRN_coef_score%s.txt.gz', args$min_chip_score)), sep="\t")

##################
# Print warnings #
##################

warnings()