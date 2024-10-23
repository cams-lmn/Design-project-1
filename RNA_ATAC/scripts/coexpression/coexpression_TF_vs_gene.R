############################################
##                                        ##
##  coexpression_TF_vs_gene_pseudobulk.R  ##
##                                        ##
############################################

source("/data/homes/louisc/Project_Babraham/RNA_ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',             type="character",                               help='SingleCellExperiment file')
p$add_argument('--TFs_file',             type="character",                               help='txt file with a list of TFs')
p$add_argument('--cor_test',         type="character",        help='Pearson or spearman')
p$add_argument('--outdir',          type="character",                               help='Output file')
args <- p$parse_args(commandArgs(TRUE))

# # START TEST ##
# args <- list()
# args$sce <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/RNA/pseudobulk_with_replicates.rds"
# args$TFs_file <- "/data/homes/louisc/Project_Babraham/TFs.txt" # "/Users/argelagr/data/mm10_regulation/TFs/TFs.txt"
# args$cor_test <- "pearson"
# args$outdir <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/coexpression/pseudobulk/"
# # END TEST ##

#####################
## Define settings ##
#####################

print(args)

# Sanity checks
stopifnot(args$cor_test%in%c("pearson","spearman"))

# I/O
dir.create(args$outdir, showWarnings=F, recursive=T)

#########################
## Load pseudobulk RNA ##
#########################

# Load SingleCellExperiment
rna.sce <- readRDS(args$sce)

##################################################
## Split RNA expression matrix into TF vs genes ##
##################################################

TFs <- fread(args$TFs_file)[[1]]# %>% unique

print("The following TFs are not found in the SingleCellExperiment:")
print(TFs[!TFs%in%toupper(rownames(rna.sce))])

rna_tfs.sce <- rna.sce[toupper(rownames(rna.sce))%in%TFs,]
rna_genes.sce <- rna.sce
rownames(rna_tfs.sce) <- toupper(rownames(rna_tfs.sce))

print(sprintf("Number of TFs: %s",nrow(rna_tfs.sce)))
print(sprintf("Number of genes: %s",nrow(rna_genes.sce)))

##########################
## Correlation analysis ##
##########################

if (args$cor_test=="pearson") {
  tf2gene_cor.mtx <- cor(t(logcounts(rna_tfs.sce)),t(logcounts(rna_genes.sce))) %>% round(2)
  tf2tf_cor.mtx <- cor(t(logcounts(rna_tfs.sce)),t(logcounts(rna_tfs.sce))) %>% round(2)
} else if (args$cor_test=="spearman") {
  tf2gene_cor.mtx <- cor(t(logcounts(rna_tfs.sce)), t(logcounts(rna_genes.sce)), method = "spearman")
  tf2tf_cor.mtx <- cor(t(logcounts(rna_tfs.sce)), t(logcounts(rna_tfs.sce)), method = "spearman")
}

##########
## Save ##
##########

outfile <- sprintf("correlation_matrix_%s_tf2gene",args$cor_test)
outfile <- paste0(outfile,".rds")
saveRDS(tf2gene_cor.mtx, file.path(args$outdir,outfile))

outfile <- sprintf("correlation_matrix_%s_tf2tf",args$cor_test)
outfile <- paste0(outfile,".rds")
saveRDS(tf2tf_cor.mtx, file.path(args$outdir,outfile))

##########
## TEST ##
##########

# tf2gene_cor.mtx <- readRDS(file.path(args$outdir,"correlation_matrix_tf2gene.rds"))

# i <- "FOXA2"
# j <- "Cab39l"

# to.plot <- data.table(
#   TF = logcounts(rna_tfs.sce[i,])[1,],
#   target_gene = logcounts(rna_genes.sce[j,])[1,],
#   celltype = colnames(rna_tfs.sce)
# )


# ggscatter(to.plot, x="TF", y="target_gene", fill="celltype", size=4, shape=21, 
#           add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
#   stat_cor(method = "pearson") +
#   scale_fill_manual(values=opts$celltype.colors) +
#   labs(x=sprintf("%s expression",i), y=sprintf("%s expression",j)) +
#   guides(fill=F) +
#   theme(
#     plot.title = element_text(hjust = 0.5, size=rel(0.85)),
#     axis.text = element_text(size=rel(0.7))
#   )
