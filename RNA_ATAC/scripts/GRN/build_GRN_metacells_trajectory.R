##########################################
##                                      ##
##  Build_GRN_metacells_trajectories.R  ##
##                                      ##
##########################################

source("/data/homes/louisc/Project_Babraham/RNA_ATAC/scripts/Settings.R")

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
args <- p$parse_args(commandArgs(TRUE))

# # I/O
# args <- list()
# args$trajectory_name <- "N2P"
# args$sce <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/metacells/SingleCellExperiment_metacells_nodiff.rds"
# args$tf2gene_virtual_chip <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/gene_regulatory_network/metacells/CISBP/TF2gene_after_virtual_chip.txt.gz"
# args$outdir <-  "/data/homes/louisc/Project_Babraham/RNA_ATAC/gene_regulatory_network/metacells/trajectories/N2P"
# args$min_chip_score <- 0.06
# args$max_distance <- 5e4
# args$ncores <- 4

print(args)

##################
## GRN_function ##
##################

source("/data/homes/louisc/Project_Babraham/RNA_ATAC/scripts/GRN/build_GRN_functions.R")

#####################
## Define settings ##
#####################

# print(paste0("Method for GRN coef calculation: ",args$GRN_method))

if(grepl("full",args$outdir)){
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

# tf2gene_chip.dt[tf=="T"] %>% View
print(head(tf2gene_chip.dt))
# print(head(tf2gene_chip.dt[tf2gene_chip.dt$tf=="NANOG"]))
# print(head(tf2gene_chip.dt[tf2gene_chip.dt$tf=="SOX2"]))
# print(head(tf2gene_chip.dt[tf2gene_chip.dt$tf=="TFCP2L1"]))

##############################
## Load RNA expression data ##
##############################

sce_mc <- readRDS(args$sce)

# (Optional) restrict to marker genes
# marker_genes.dt <- fread(args$rna.atlas.marker_genes)
# sce_mc <- sce_mc[rownames(sce_mc)%in%unique(marker_genes.dt$gene),]

##########################
## Filter TFs and genes ##
##########################

TFs <- intersect(unique(tf2gene_chip.dt$tf),toupper(rownames(sce_mc)))
print(head(TFs))
genes <- intersect(unique(tf2gene_chip.dt$gene),rownames(sce_mc))

if (grepl("full",args$outdir)){
  tf2gene_chip.dt <- tf2gene_chip.dt[tf%in%TFs & gene%in%genes,]
} else {
  tf2gene_chip.dt <- tf2gene_chip.dt[tf%in%TFs & gene%in%TFs,]
}

# Fetch RNA expression matrices
rna_tf.mtx <- logcounts(sce_mc)[unique(tf2gene_chip.dt$tf),]
rownames(rna_tf.mtx) <- toupper(rownames(rna_tf.mtx))
rna_targets.mtx <- logcounts(sce_mc)[unique(tf2gene_chip.dt$gene),]

print(dim(rna_tf.mtx))
# print(head(rna_tf.mtx))
print(dim(rna_targets.mtx))
print(head(rna_targets.mtx[,1:6]))

# Filter out lowly variable genes and TFs
# rna_tf.mtx <- rna_tf.mtx[apply(rna_tf.mtx,1,var)>=1,]
# rna_targets.mtx <- rna_targets.mtx[apply(rna_targets.mtx,1,var)>=0.1,]
print(dim(rna_tf.mtx))
print(head(rna_tf.mtx[,1:6]))

TFs <- intersect(unique(tf2gene_chip.dt$tf),rownames(rna_tf.mtx))
genes <- intersect(unique(tf2gene_chip.dt$gene),rownames(rna_targets.mtx))
tf2gene_chip.dt <- tf2gene_chip.dt[tf%in%TFs & gene%in%genes,]

print("Filtering for marker TF")
TF_markers <- read.table(args$markers_TF,header=T,sep="\t")$Gene
print(sum(TFs%in%TF_markers))
# tf2gene_chip.dt <- tf2gene_chip.dt[tf%in%TF_markers & gene%in%TF_markers,]
print(dim(tf2gene_chip.dt))
# tf2gene_chip.dt <- tf2gene_chip.dt[tf%in%TF_markers,]
# print(dim(tf2gene_chip.dt))
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
list_results <- mclapply(genes,myfun_conv,mc.cores=args$ncores,mc.preschedule=T)
GRN_coef.dt <- do.call(rbind,list_results)

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

fwrite(GRN_coef.dt, file.path(args$outdir,sprintf('global_chip_GRN_coef_score%s.txt.gz',args$min_chip_score)), sep="\t")

##################
# Print warnings #
##################

warnings()