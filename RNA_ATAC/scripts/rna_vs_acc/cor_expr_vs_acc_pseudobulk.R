####################################
##                                ##
##  cor_expr_vs_acc_pseudobulk.R  ##
##                                ##
####################################

source("/data/homes/louisc/Project_Babraham/RNA/scripts/Settings.R")

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--sce',  type="character",              help='RNA SingleCellExperiment (pseudobulk)')
p$add_argument('--atac_peak_matrix',  type="character",              help='ATAC Peak matrix (pseudobulk)')
p$add_argument('--peak2gene',  type="character",              help='Peak2gene file')
p$add_argument('--atac_peak_metadata',  type="character",              help='Peak metadata')
p$add_argument('--distance',  type="integer",            default=1e5,      help='Maximum distance for a linkage between a peak and a gene')
p$add_argument('--outdir',          type="character",                help='Output directory')
p$add_argument('--min_cor',  type="double",            default=0.25,      help='Minimal correlation for correlations of interest')
p$add_argument('--max_pval',  type="double",            default=0.01,      help='Maximal pvalue for correlations of interest')
p$add_argument('--force_rerun',  type="logical",            default=FALSE,      help='Should a rerun be forced if correlation file already exists')
p$add_argument('--which_clusters',  type="character",            default=NULL,      help='Which clusters should be used for analysis')
p$add_argument('--n_clusters',  type="integer",            default=7,      help='Amount of gene clusters')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

# ## START TEST ##
# args <- list()
# args$sce <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/RNA/pseudobulk_with_replicates.rds"
# args$atac_peak_matrix <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/PeakMatrix/pseudobulk_with_replicates.rds"
# args$peak2gene <- "/data/homes/louisc/Project_Babraham/ATAC/archR/PeakCalls/peaks2genes/peaks2genes_all.txt.gz"
# args$atac_peak_metadata <- "/data/homes/louisc/Project_Babraham/ATAC/archR/peak_metadata.tsv.gz"
# args$distance <- 5e4
# args$max_pval <- 0.01
# args$min.cor <- c(0.25,0.5,0.75)
# args$force_rerun <- FALSE
# args$remove_clusters <- NULL
# # args$remove_clusters <- 6
# # args$remove_clusters <- c(2;3)
# # args$remove_clusters <- c(1,4,5,6)
# ## END TEST ##

#####################
## Define settings ##
#####################

print("Pseudobulk")

# I/O
args$outfile <- paste0(args$outdir,args$which_clusters,"/cor_gene_expr_vs_peak_acc_pseudobulk.txt.gz")
args$outdir <- dirname(args$outfile)

unlink(sprintf("%s/completed*",args$outdir))

dir.create(args$outdir, showWarnings = F, recursive = TRUE)

print(sprintf("Pvalue cutoff: %s",args$max_pval))
print(sprintf("Correlation cutoff: %s",args$min_cor))

if (args$which_clusters!="all"){
  args$remove_clusters <- gsub("-","",gsub("^no","",args$which_clusters))
  print(paste0("Removing the following clusters: ",paste(args$remove_clusters,collapse=", ")))
}

cols3 <- rev(brewer.pal(10, "Spectral"))

n_clusters <- args$n_clusters

##################################
## Load pseudobulk RNA and ATAC ##
##################################

# Load SingleCellExperiment
rna_pseudobulk.sce <- readRDS(args$sce)
rna_pseudobulk.sce

# Load ATAC SummarizedExperiment
atac_peakMatrix_pseudobulk.se <- readRDS(args$atac_peak_matrix)
rowData(atac_peakMatrix_pseudobulk.se) <- NULL
atac_peakMatrix_pseudobulk.se

# Normalise ATAC data
assayNames(atac_peakMatrix_pseudobulk.se) <- "counts"
assay(atac_peakMatrix_pseudobulk.se,"logcounts") <- log(1e6*(sweep(assay(atac_peakMatrix_pseudobulk.se),2,colSums(assay(atac_peakMatrix_pseudobulk.se),na.rm=T),"/"))+1)

# Make sure that samples are consistent
celltypes <- intersect(colnames(rna_pseudobulk.sce),colnames(atac_peakMatrix_pseudobulk.se))
if (!is.null(args$remove_clusters)){
  celltypes <- celltypes[!(celltypes%in%args$remove_clusters)]
}
print(celltypes)
rna_pseudobulk.sce <- rna_pseudobulk.sce[,celltypes]
atac_peakMatrix_pseudobulk.se <- atac_peakMatrix_pseudobulk.se[,celltypes]

#############################
## Load peak2gene linkages ##
#############################

peak2gene.dt <- fread(args$peak2gene) %>%
  .[dist<=args$distance]

# peak2gene.dt <- fread(io$archR.peak2gene.nearest) %>%
#   .[dist<=args$distance]

###################
## Density plots ##
###################

print("Density")

to.plot.density_rna <- cbind(logcounts(rna_pseudobulk.sce)%>% melt(),rep("RNA",ncol(rna_pseudobulk.sce)*nrow(rna_pseudobulk.sce)))         
to.plot.density_atac <- cbind(assay(atac_peakMatrix_pseudobulk.se,"logcounts")[peak2gene.dt$peak,]%>% melt(),
                         rep("ATAC",ncol(assay(atac_peakMatrix_pseudobulk.se,"logcounts")[peak2gene.dt$peak,])*nrow(atac_peakMatrix_pseudobulk.se[peak2gene.dt$peak,])))
colnames(to.plot.density_rna) <- c("feature","sample","logcount","modality")
colnames(to.plot.density_atac) <- c("feature","sample","logcount","modality")
to.plot.density <- rbind(to.plot.density_rna,to.plot.density_atac)
to.plot.density$sample <- as.factor(to.plot.density$sample)
to.plot.density$modality <- as.factor(to.plot.density$modality)
print(head(to.plot.density))
print(summary(to.plot.density))

p <- ggplot(data=to.plot.density) +
  geom_density(mapping = aes(x=logcount, color=sample)) + 
  facet_wrap(~ modality)

pdf(file = sprintf("%s/density_plot.pdf", args$outdir),
    height = 6, width = 8)
print(p)
dev.off()

###########################################################
## Correlate peak accessibility with gene RNA expression ##
###########################################################

# TFs <- intersect(colnames(motifmatcher.se),rownames(rna_pseudobulk.sce))# %>% head(n=5)
# TFs <- c("GATA1","TAL1")
genes <- intersect(rownames(rna_pseudobulk.sce),peak2gene.dt$gene)

if (args$force_rerun | !file.exists(args$outfile)){
  cor.dt <- genes %>% map(function(i) {
    tmp <- peak2gene.dt[gene==i,c("gene","peak")]
    print(sprintf("%s (%d/%d)",i,match(i,genes),length(genes)))
    # calculate correlations
    corr_output <- psych::corr.test(
      x = t(logcounts(rna_pseudobulk.sce[i,])), 
      y = t(logcounts(atac_peakMatrix_pseudobulk.se[tmp$peak,])), 
      ci = FALSE
    )
    
    data.table(
      gene = i,
      peak = tmp$peak,
      cor = corr_output$r[1,] %>% round(3),
      pvalue = corr_output$p[1,] %>% format(digits=3)
    )
  }) %>% rbindlist
  
  fwrite(cor.dt, args$outfile, sep="\t")
  cor.dt <- fread(args$outfile) # read in again to be sure the pval column is numerical
  print(head(cor.dt))
} else {
  cor.dt <- fread(args$outfile)
  print(head(cor.dt))
}

####################
## Transform data ##
####################

gene_expr_vs_peak_acc.dt <- cor.dt %>%
  .[is.na(cor), c("cor","pvalue"):=list(0,1)] %>% .[,dist:=NULL]

# Filter genes
genes.to.use <- grep("*^MT-|^RPS|^RPL|^OR[0-9]",unique(gene_expr_vs_peak_acc.dt$gene), invert=T,value=T)
head(genes.to.use)
length(genes.to.use)
gene_expr_vs_peak_acc.dt <- gene_expr_vs_peak_acc.dt[gene%in%genes.to.use]

########################
## Load peak metadata ##
########################

peak_metadata.dt <- fread(args$atac_peak_metadata) %>%
  .[,peak:=sprintf("%s:%s-%s",chr,start,end)] %>%
  .[,c("peak","peakType","distToGeneStart","nearestGene")]

############################
## Select peak2gene links ##
############################

peak2gene.dt <- peak2gene.dt %>% 
  .[gene%in%genes.to.use] %>%
  .[,c("peak","gene","dist")] %>%
  merge(peak_metadata.dt[,c("peak","peakType")])

# Note that this does not contain peaks with no links to genes?
length(unique(peak2gene.dt$peak))

stopifnot(unique(peak2gene.dt$peak)%in%unique(peak_metadata.dt$peak))

###############################################
## Calculate number of peaks linked to genes ##
###############################################

peak_metadata.dt[distToGeneStart<=args$distance,.N] / nrow(peak_metadata.dt)

###########
## Merge ##
###########

gene_expr_vs_peak_acc_filt.dt <- gene_expr_vs_peak_acc.dt %>% 
  merge(peak2gene.dt,by=c("peak","gene")) %>%
  .[,sign:=c("-","+")[as.numeric(cor>0)+1]]

print(dim(gene_expr_vs_peak_acc_filt.dt))
print(str(gene_expr_vs_peak_acc_filt.dt))
print(head(gene_expr_vs_peak_acc_filt.dt))

length(unique(gene_expr_vs_peak_acc_filt.dt$peak))

sum(gene_expr_vs_peak_acc_filt.dt$pvalue<=args$max_pval & abs(gene_expr_vs_peak_acc_filt.dt$cor)>=args$min_cor)

###########################################################################################
## Calculate number of peaks that are correlated with RNA expression of the nearest gene ##
###########################################################################################

tmp <- gene_expr_vs_peak_acc_filt.dt %>% .[,sum(pvalue<=args$max_pval & abs(cor)>=args$min_cor), by=c("peak")] %>% .[,mean(V1>=1)]

to.plot <- gene_expr_vs_peak_acc_filt.dt %>%
  .[,sum(pvalue<=args$max_pval & abs(cor)>=args$min_cor), by=c("peak","sign","peakType")] %>%
  .[,100*mean(V1>=1),by=c("sign","peakType")]

p <- ggbarplot(to.plot, x="peakType", fill="sign", y="V1", position=position_dodge(width = 0.75)) +
  labs(x="", y="% of peaks correlated with proximal gene expr.") +
  # scale_fill_brewer(palette="Dark2") +
  # guides(color=F) + 
  theme(
    axis.title = element_text(size=rel(0.80)),
    legend.title = element_blank(),
    axis.text.x = element_text(size=rel(0.80), color="black"),
    axis.text.y = element_text(size=rel(0.80), color="black")
  )

pdf(file.path(args$outdir,"barplots_number_gene_expr_vs_peak_acc_correlations.pdf"), width=7, height=5)
print(p)
dev.off()

####################################################
## Histogram of ATAC peak acc vs RNA correlations ##
####################################################

to.plot <- gene_expr_vs_peak_acc_filt.dt# %>% .[pvalue<=0.05]

p <- gghistogram(to.plot, x="cor", y="..density..", bins=50, fill="gray70") +
  labs(x="Pearson correlation coefficient\n(ATAC peak vs gene expr.)", y="Density") +
  theme(
    axis.title.x = element_text(size=rel(0.9)),
    axis.title.y = element_text(size=rel(1.0)),
    legend.title = element_blank(),
    axis.text.x = element_text(size=rel(0.8), color="black"),
    axis.text.y = element_text(size=rel(0.8), color="black")
  )

pdf(file.path(args$outdir,"histogram_gene_expr_vs_peak_acc_corr.pdf"), width=7, height=5)
print(p)
dev.off()

to.plot <- gene_expr_vs_peak_acc_filt.dt %>% .[pvalue<=0.05]

p <- gghistogram(to.plot, x="cor", y="..density..", bins=50, fill="gray70") +
  labs(x="Pearson correlation coefficient\n(ATAC peak vs gene expr.)", y="Density") +
  theme(
    axis.title.x = element_text(size=rel(0.9)),
    axis.title.y = element_text(size=rel(1.0)),
    legend.title = element_blank(),
    axis.text.x = element_text(size=rel(0.8), color="black"),
    axis.text.y = element_text(size=rel(0.8), color="black")
  )

pdf(file.path(args$outdir,"histogram_sign_gene_expr_vs_peak_acc_corr.pdf"), width=7, height=5)
print(p)
dev.off()

print(head(gene_expr_vs_peak_acc_filt.dt))

sum((as.numeric(gene_expr_vs_peak_acc_filt.dt$pvalue)<=args$max_pval) & (gene_expr_vs_peak_acc_filt.dt$sign=="+"))
sum((as.numeric(gene_expr_vs_peak_acc_filt.dt$pvalue)<=args$max_pval) & (gene_expr_vs_peak_acc_filt.dt$sign=="+")) /nrow(gene_expr_vs_peak_acc_filt.dt)

sum((as.numeric(gene_expr_vs_peak_acc_filt.dt$pvalue)<=args$max_pval) & (gene_expr_vs_peak_acc_filt.dt$sign=="-"))
sum((as.numeric(gene_expr_vs_peak_acc_filt.dt$pvalue)<=args$max_pval) & (gene_expr_vs_peak_acc_filt.dt$sign=="-")) /nrow(gene_expr_vs_peak_acc_filt.dt)

##############################
## Number of peaks per gene ##
##############################

tmp <- peak2gene.dt %>% .[dist<=args$distance] %>% .[,.(npeaks=.N),by=c("gene")]

to.plot <- peak2gene.dt %>% .[dist<=args$distance] %>%
  .[,.(npeaks=.N),by=c("gene")] %>%
  .[,.N,by="npeaks"] %>%
  .[npeaks>=2] %>%
  .[npeaks>=200,npeaks:=200] # for viz purposes

p <- ggline(to.plot, x="npeaks", y="N") +
  labs(x="Number of peaks per gene", y="Number of genes") +
  theme(
    axis.title = element_text(size=rel(1.0), color="black"),
    axis.text = element_text(size=rel(0.80), color="black")
  )

pdf(file.path(args$outdir,"lineplot_npeaks_per_gene.pdf"), width=7, height=5)
print(p)
dev.off()

to.plot <- peak2gene.dt %>% .[dist<=args$distance] %>% .[,.(npeaks=.N),by=c("gene")]
p <- gghistogram(to.plot, x="npeaks", y="..density..", bins=50, fill="gray70") +
  scale_x_continuous(limits=c(0,100)) +
  labs(x="Number of peaks per gene", y="Density") +
  theme(
    axis.title = element_text(size=rel(1.0), color="black"),
    axis.text = element_text(size=rel(1.0), color="black")
  )

pdf(file.path(args$outdir,"histogram_npeaks_per_gene.pdf"), width=7, height=5)
print(p)
dev.off()

mean(to.plot$npeaks)

##############################
## Number of genes per peak ##
##############################

tmp <- peak2gene.dt %>% .[,.(ngenes=.N),by=c("peak")]
mean(tmp$ngenes)

to.plot <- peak2gene.dt %>%
  .[,.(ngenes=.N),by=c("peak")] %>%
  .[,.N,by="ngenes"] %>%
  .[ngenes<=30]

p <- ggline(to.plot, x="ngenes", y="N") +
  labs(x="Number of genes per peak", y="Number of peaks") +
  theme(
    axis.title = element_text(size=rel(1.0), color="black"),
    axis.text = element_text(size=rel(0.80), color="black")
  )

pdf(file.path(args$outdir,"lineplot_ngenes_per_peak.pdf"), width=7, height=5)
print(p)
dev.off()

to.plot <- peak2gene.dt %>% .[,.(ngenes=.N),by=c("peak")]
p <- gghistogram(to.plot, x="ngenes", y="..density..", bins=50, fill="gray70") +
  scale_x_continuous(limits=c(0,100)) +
  labs(x="Number of peaks per gene", y="Density") +
  theme(
    axis.title = element_text(size=rel(1.0), color="black"),
    axis.text = element_text(size=rel(1.0), color="black")
  )

pdf(file.path(args$outdir,"histogram_ngenes_per_peak.pdf"), width=7, height=5)
print(p)
dev.off()

mean(to.plot$ngenes)

##########################################################
## Expression profile of significantly correlated genes ##
##########################################################

genes_neg_corr <- unique(gene_expr_vs_peak_acc_filt.dt$gene[abs(gene_expr_vs_peak_acc_filt.dt$cor)>=args$min_cor &
                                                              gene_expr_vs_peak_acc_filt.dt$pval<args$max_pval &
                                                              gene_expr_vs_peak_acc_filt.dt$sign=="-"])
print(length(genes_neg_corr))

genes_pos_corr <- unique(gene_expr_vs_peak_acc_filt.dt$gene[abs(gene_expr_vs_peak_acc_filt.dt$cor)>=args$min_cor &
                                                              gene_expr_vs_peak_acc_filt.dt$pval<args$max_pval &
                                                              gene_expr_vs_peak_acc_filt.dt$sign=="+"])
print(length(genes_pos_corr))

genes_purely_neg_corr <- genes_neg_corr[!(genes_neg_corr%in%genes_pos_corr)]
print(length(genes_purely_neg_corr))
genes_purely_pos_corr <- genes_pos_corr[!(genes_pos_corr%in%genes_neg_corr)]
print(length(genes_purely_pos_corr))

# Negative correlation
print("Negative correlation")
if(length(genes_neg_corr)>1){
  sce <- readRDS(args$sce)
  print(sce)
  exp.values <- NULL
  for (i in 1:5){
    if (sum(sce$celltype==i)>1) {
      exp.values <- cbind(exp.values,rowMeans(assays(sce)[["logcounts"]][rownames(sce)%in%genes_neg_corr,sce$celltype==i]))
    }
    else {
      exp.values <- cbind(exp.values,assays(sce)[["logcounts"]][rownames(sce)%in%genes_neg_corr,sce$celltype==i])
    }
  }
  colnames(exp.values) <- paste0(rep("cluster",5),1:5)
  head(exp.values)
  
  t(apply(exp.values, 1, z_score)) %>%
    as_tibble(rownames = "Gene") %>%
    column_to_rownames(var="Gene") -> exp.values.scaled
  exp.values.scaled <- exp.values.scaled[!(rowSums(is.na(exp.values.scaled))>0),]
  head(exp.values.scaled)
  dim(exp.values.scaled)
  summary(exp.values.scaled)
  
  set.seed(42)
  Rclusterpp.hclust(exp.values.scaled) -> clustering.genes
  
  exp.values.scaled_multi <- exp.values.scaled %>% rownames_to_column(var="gene") %>% 
    merge(cor.dt)
  print(dim(exp.values.scaled_multi))
  exp.values.scaled_multi <- exp.values.scaled_multi[exp.values.scaled_multi$pvalue<=args$max_pval & 
                                                       exp.values.scaled_multi$cor<=(-args$min_cor),] %>% 
    select(-c("cor","pvalue")) 
  print(dim(exp.values.scaled_multi))
  
  atac.se <- readRDS(args$atac_peak_matrix) 
  exp.values_atac <- NULL
  for (i in 1:5){
    if (sum(sce$celltype==i)>1){
      exp.values_atac <- cbind(exp.values_atac,rowMeans(assays(atac.se)[["logcounts"]][rownames(atac.se)%in%exp.values.scaled_multi$peak,sce$celltype==i]))
    } else {
      exp.values_atac <- cbind(exp.values_atac,assays(atac.se)[["logcounts"]][rownames(atac.se)%in%exp.values.scaled_multi$peak,sce$celltype==i])
    }
  }
  colnames(exp.values_atac) <- paste0(rep("atac_cluster",5),1:5)
  
  t(apply(exp.values_atac, 1, z_score)) %>%
    as_tibble(rownames = "Gene") %>%
    column_to_rownames(var="Gene") -> exp.values_atac.scaled
  exp.values_atac.scaled <- exp.values_atac.scaled[!(rowSums(is.na(exp.values.scaled))>0),]
  head(exp.values_atac.scaled)
  dim(exp.values_atac.scaled)
  summary(exp.values_atac.scaled)
  
  exp.values.scaled_multi <- exp.values_atac.scaled %>% rownames_to_column(var="peak") %>% 
    merge(exp.values.scaled_multi) %>% select(-peak) %>% dplyr::rename(Gene = gene)
  
  count.gaps <- function(x,clusters){
    clusters %>%
      filter(cluster<=x) %>%
      nrow -> y
    return(y)
  }
  
  cutree(clustering.genes, k=n_clusters) %>%
    as_tibble(rownames = "Gene") %>%
    dplyr::rename(cluster = value) -> clusters.n
  
  lapply(c(1:n_clusters), count.gaps, clusters.n %>% merge(exp.values.scaled_multi)) %>%
    do.call("rbind",.) %>%
    as.vector -> clusters.gaps.n
  
  annotation.colors <- list(cluster = c("1"=cols3[1], "2"=cols3[2], "3"=cols3[3], 
                                        "4"=cols3[4], "5"=cols3[5], "6"=cols3[6], 
                                        "7"=cols3[7], "8"=cols3[8], "9"=cols3[9], 
                                        "10"=cols3[10])[1:n_clusters])
  
  pdf(file = sprintf("%s/heatmap_%sclusters_negcorr_corr%s.pdf",
                     args$outdir,n_clusters,args$min_cor), height = 6, width = 8)
  exp.values.scaled_multi %>%
    arrange(Gene) %>%
    full_join(clusters.n) %>%
    arrange(cluster) %>%
    select(-c("cluster","Gene")) %>%
    pheatmap(cluster_rows = F,
             cluster_cols = F,
             show_rownames = F,
             gaps_row = clusters.gaps.n,
             gaps_col = 5,
             # annotation_row = (clusters.n %>% column_to_rownames(var="Gene") %>% arrange(cluster)),
             annotation_colors = annotation.colors,
             color = colorRampPalette(cols3)(100),breaks=seq(from=-2,to=2,length.out=101))
  dev.off()
  
  file.create(sprintf("%s/completed_negcorr_corr%s_%spairs.pdf",args$outdir,args$min_cor,nrow(exp.values.scaled_multi)))
  file.create(sprintf("%s/neg_corr_completed_%sclusters_corr%s.pdf",args$outdir,n_clusters,args$min_cor))
} else {
  print("No significant correlations with current cutoffs.")
  file.create(sprintf("%s/neg_corr_completed_%sclusters_corr%s.pdf",args$outdir,n_clusters,args$min_cor))
}                        

# Positive correlation
print("Positive correlation")
if(length(genes_pos_corr)>1){
  sce <- readRDS(args$sce)
  print(sce)
  exp.values <- NULL
  for (i in 1:5){
    if (sum(sce$celltype==i)>1) {
      exp.values <- cbind(exp.values,rowMeans(assays(sce)[["logcounts"]][rownames(sce)%in%genes_pos_corr,sce$celltype==i]))
    }
    else {
      exp.values <- cbind(exp.values,assays(sce)[["logcounts"]][rownames(sce)%in%genes_pos_corr,sce$celltype==i])
    }
  }
  colnames(exp.values) <- paste0(rep("cluster",5),1:5)
  head(exp.values)
  
  t(apply(exp.values, 1, z_score)) %>%
    as_tibble(rownames = "Gene") %>%
    column_to_rownames(var="Gene") -> exp.values.scaled
  exp.values.scaled <- exp.values.scaled[!(rowSums(is.na(exp.values.scaled))>0),]
  head(exp.values.scaled)
  dim(exp.values.scaled)
  summary(exp.values.scaled)
  
  set.seed(42)
  Rclusterpp.hclust(exp.values.scaled) -> clustering.genes
  
  exp.values.scaled_multi <- exp.values.scaled %>% rownames_to_column(var="gene") %>% 
    merge(cor.dt)
  print(dim(exp.values.scaled_multi))
  exp.values.scaled_multi <- exp.values.scaled_multi[exp.values.scaled_multi$pvalue<=args$max_pval & 
                                                       exp.values.scaled_multi$cor>=args$min_cor,] %>% 
    select(-c("cor","pvalue")) 
  print(dim(exp.values.scaled_multi))
  
  atac.se <- readRDS(args$atac_peak_matrix) 
  exp.values_atac <- NULL
  for (i in 1:5){
    if (sum(sce$celltype==i)>1){
      exp.values_atac <- cbind(exp.values_atac,rowMeans(assays(atac.se)[["logcounts"]][rownames(atac.se)%in%exp.values.scaled_multi$peak,sce$celltype==i]))
    } else {
      exp.values_atac <- cbind(exp.values_atac,assays(atac.se)[["logcounts"]][rownames(atac.se)%in%exp.values.scaled_multi$peak,sce$celltype==i])
    }
  }
  head(exp.values)
  colnames(exp.values_atac) <- paste0(rep("atac_cluster",5),1:5)
  
  t(apply(exp.values_atac, 1, z_score)) %>%
    as_tibble(rownames = "Gene") %>%
    column_to_rownames(var="Gene") -> exp.values_atac.scaled
  exp.values_atac.scaled <- exp.values_atac.scaled[!(rowSums(is.na(exp.values.scaled))>0),]
  head(exp.values_atac.scaled)
  dim(exp.values_atac.scaled)
  summary(exp.values_atac.scaled)
  
  exp.values.scaled_multi <- exp.values_atac.scaled %>% rownames_to_column(var="peak") %>% 
    merge(exp.values.scaled_multi) %>% select(-peak) %>% dplyr::rename(Gene = gene)
  
  count.gaps <- function(x,clusters){
    clusters %>%
      filter(cluster<=x) %>%
      nrow -> y
    return(y)
  }
  
  cutree(clustering.genes, k=n_clusters) %>%
    as_tibble(rownames = "Gene") %>%
    dplyr::rename(cluster = value) -> clusters.n
  
  lapply(c(1:n_clusters), count.gaps, clusters.n %>% merge(exp.values.scaled_multi)) %>%
    do.call("rbind",.) %>%
    as.vector -> clusters.gaps.n
  
  annotation.colors <- list(cluster = c("1"=cols3[1], "2"=cols3[2], "3"=cols3[3], 
                                        "4"=cols3[4], "5"=cols3[5], "6"=cols3[6], 
                                        "7"=cols3[7], "8"=cols3[8], "9"=cols3[9], 
                                        "10"=cols3[10])[1:n_clusters])
  
  pdf(file = sprintf("%s/heatmap_%sclusters_poscorr_corr%s.pdf",
                     args$outdir,n_clusters,args$min_cor), height = 6, width = 8)
  exp.values.scaled_multi %>%
    arrange(Gene) %>%
    full_join(clusters.n) %>%
    arrange(cluster) %>%
    select(-c("cluster","Gene")) %>%
    pheatmap(cluster_rows = F,
             cluster_cols = F,
             show_rownames = F,
             gaps_row = clusters.gaps.n,
             gaps_col = 5,
             # annotation_row = (clusters.n %>% column_to_rownames(var="Gene") %>% arrange(cluster)),
             annotation_colors = annotation.colors,
             color = colorRampPalette(cols3)(100),breaks=seq(from=-2,to=2,length.out=101))
  dev.off()

  file.create(sprintf("%s/completed_poscorr_corr%s_%spairs.pdf",args$outdir,args$min_cor,nrow(exp.values.scaled_multi)))
  file.create(sprintf("%s/pos_corr_completed_%sclusters_corr%s.pdf",args$outdir,n_clusters,args$min_cor))
} else {
  print("No significant correlations with current cutoffs.")
  file.create(sprintf("%s/pos_corr_completed_%sclusters_corr%s.pdf",args$outdir,n_clusters,args$min_cor))
}    
