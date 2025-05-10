###########################################
##                                       ##
##   Cor_TFexpr_vs_peakAcc_pseudobulk.R  ##
##                                       ##
###########################################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA_ATAC/scripts/Settings.R")

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--sce',  type="character",              help='RNA SingleCellExperiment (pseudobulk)')
p$add_argument('--atac_peak_matrix',  type="character",              help='ATAC Peak matrix (pseudobulk)')
p$add_argument('--motif_annotation',  type="character",              help='Motif annotation')
p$add_argument('--motif2gene',  type="character",              help='Motif2gene')
p$add_argument('--motifmatcher',  type="character",              help='Motif annotation')
p$add_argument('--force_rerun',  type="logical",            default=FALSE,      help='Should a rerun be forced if correlation file already exists')
p$add_argument('--which_clusters',  type="character",            default=NULL,      help='Should a rerun be forced if correlation file already exists')
p$add_argument('--outdir',          type="character",                help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

args$outfile <- paste0(args$outdir,args$which_clusters,"/",args$motif_annotation,"_cor_TFexpr_vs_peakAcc_pseudobulk.txt.gz")
print(args$outfile)
args$outdir <- dirname(args$outfile)

# I/O
dir.create(args$outdir, showWarnings = F, recursive = T)

if (args$which_clusters!="all"){
  args$remove_clusters <- gsub("-","",gsub("^no","",args$which_clusters))
  print(paste0("Removing the following clusters: ",paste(args$remove_clusters,collapse=", ")))
}

cols3 <- rev(brewer.pal(10, "Spectral"))

#######################
## Load RNA and ATAC ##
#######################

# Load SingleCellExperiment
rna_pseudobulk.sce <- readRDS(args$sce)
rna_pseudobulk.sce

# Load ATAC SummarizedExperiment
atac_peakMatrix_pseudobulk.se <- readRDS(args$atac_peak_matrix)
atac_peakMatrix_pseudobulk.se
# head(assay(atac_peakMatrix_pseudobulk.se,"logcounts"))
# head(logcounts(atac_peakMatrix_pseudobulk.se))


# Normalise ATAC data
# assayNames(atac_peakMatrix_pseudobulk.se) <- "counts"
# assay(atac_peakMatrix_pseudobulk.se,"logcounts") <- log(1e6*(sweep(assay(atac_peakMatrix_pseudobulk.se),2,colSums(assay(atac_peakMatrix_pseudobulk.se),na.rm=T),"/"))+1)

# Make sure that samples are consistent
samples <- intersect(colnames(rna_pseudobulk.sce),colnames(atac_peakMatrix_pseudobulk.se))
if (!is.null(args$remove_clusters)){
  samples <- samples[!(samples%in%args$remove_clusters)]
}
print(samples)
rna_pseudobulk.sce <- rna_pseudobulk.sce[,samples]
atac_peakMatrix_pseudobulk.se <- atac_peakMatrix_pseudobulk.se[,samples]

###############################
## Load motifmatcher results ##
###############################

motifmatcher.se <- readRDS(args$motifmatcher)

# Subset peaks
stopifnot(sort(rownames(motifmatcher.se))==sort(rownames(atac_peakMatrix_pseudobulk.se)))
motifmatcher.se <- motifmatcher.se[rownames(atac_peakMatrix_pseudobulk.se),]

################################
## Load motif2gene annotation ##
################################

motif2gene.dt <- fread(args$motif2gene)

################
## Filter TFs ##
################

motifs <- intersect(colnames(motifmatcher.se),motif2gene.dt$motif)
motifmatcher.se <- motifmatcher.se[,motifs]
motif2gene.dt <- motif2gene.dt[motif%in%motifs]

genes <- intersect(rownames(rna_pseudobulk.sce),motif2gene.dt$gene)
# rna_tf_pseudobulk.sce <- rna_pseudobulk.sce[str_to_title(genes),] this is when pseudobulk genes are in lower case?
# rownames(rna_tf_pseudobulk.sce) <- toupper(rownames(rna_tf_pseudobulk.sce))
print("TFs")
print(length(genes))
print(head(sort(genes)))
print(tail(sort(genes)))
rna_tf_pseudobulk.sce <- rna_pseudobulk.sce[genes,]
motif2gene.dt <- motif2gene.dt[gene%in%genes]

# Manually remove some motifs
# motif2gene.dt <- motif2gene.dt[!motif%in%c("T_789")] 

# Remove duplicated gene-motif pairs
genes.to.remove <- names(which(table(motif2gene.dt$gene)>1))
cat(sprintf("Removing %d TFs that have duplicated gene-motif pairs:\n%s", length(genes.to.remove), paste(genes.to.remove, collapse=", ")))
motif2gene.dt <- motif2gene.dt[!gene%in%genes.to.remove]
rna_tf_pseudobulk.sce <- rna_tf_pseudobulk.sce[rownames(rna_tf_pseudobulk.sce)%in%motif2gene.dt$gene]
motifmatcher.se <- motifmatcher.se[,colnames(motifmatcher.se)%in%motif2gene.dt$motif]
stopifnot(table(motif2gene.dt$gene)==1)

#########################################################
## Correlate peak accessibility with TF RNA expression ##
#########################################################

# Sanity checks
stopifnot(colnames(rna_tf_pseudobulk.sce)==colnames(atac_peakMatrix_pseudobulk.se))

TFs <- rownames(rna_tf_pseudobulk.sce)

# Prepare output data objects
cor.mtx <- matrix(as.numeric(NA), nrow=nrow(atac_peakMatrix_pseudobulk.se), ncol=length(TFs))
pvalue.mtx <- matrix(as.numeric(NA), nrow=nrow(atac_peakMatrix_pseudobulk.se), ncol=length(TFs))
rownames(cor.mtx) <- rownames(atac_peakMatrix_pseudobulk.se); colnames(cor.mtx) <- TFs
dimnames(pvalue.mtx) <- dimnames(cor.mtx)

if (args$force_rerun | !file.exists(args$outfile)){
  # i <- "GATA1"
  counter <- 1
  for (i in TFs) {
    motif_i <- motif2gene.dt[gene==i,motif]
    
    print(paste0(i,", (",counter,"/",length(TFs),") motifs done."))
    all_peaks_i <- rownames(motifmatcher.se)[which(assay(motifmatcher.se[,motif_i],"motifMatches")==1)]
    # print(head(all_peaks_i))
    
    # calculate correlations
    corr_output <- psych::corr.test(
      x = t(logcounts(rna_tf_pseudobulk.sce[i,])), 
      y = t(assay(atac_peakMatrix_pseudobulk.se,"logcounts")[all_peaks_i,]), 
      ci = FALSE
    )
    
    # Fill matrices
    cor.mtx[all_peaks_i,i] <- round(corr_output$r[1,],3)
    pvalue.mtx[all_peaks_i,i] <- round(corr_output$p[1,],10)
    
    counter <- counter + 1
  }
  
  to.save <- SummarizedExperiment(
    assays = SimpleList("cor" = dropNA(cor.mtx), "pvalue" = dropNA(pvalue.mtx)),
    rowData = rowData(atac_peakMatrix_pseudobulk.se)
  )
  saveRDS(to.save, args$outfile)
  print("Data saved")
} else {
  to.save <- readRDS(args$outfile)
  print("Data loaded")
}

head(to.save)