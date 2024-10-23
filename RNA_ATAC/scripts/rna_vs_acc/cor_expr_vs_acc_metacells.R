###################################
##                               ##
##  cor_expr_vs_acc_metacells.R  ##
##                               ##
###################################

source("/data/homes/louisc/Project_Babraham/RNA/scripts/Settings.R")

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--sce',  type="character",              help='RNA SingleCellExperiment (metacells)')
p$add_argument('--sce_pseudobulk',  type="character",              help='RNA SingleCellExperiment (pseudobulk)')
p$add_argument('--atac_peak_matrix',  type="character",              help='ATAC Peak matrix (metacells)')
p$add_argument('--atac_peak_matrix_pseudobulk',  type="character",              help='ATAC Peak matrix (pseudobulk)')
p$add_argument('--peak2gene',  type="character",              help='Peak2gene file')
p$add_argument('--atac_peak_metadata',  type="character",              help='Peak metadata')
p$add_argument('--distance',  type="integer",            default=1e5,      help='Maximum distance for a linkage between a peak and a gene')
p$add_argument('--outdir',          type="character",                help='Output directory')
p$add_argument('--min_cor',  type="double",            default=0.25,      help='Minimal correlation for correlations of interest')
p$add_argument('--max_pval',  type="double",            default=0.01,      help='Maximal pvalue for correlations of interest')
p$add_argument('--force_rerun',  type="logical",            default=FALSE,      help='Should a rerun be forced if correlation file already exists')
p$add_argument('--sort_samples',  type="logical",            default=TRUE,      help='Should samples be sorted')
p$add_argument('--filter_differentiated',  type="logical",            default=TRUE,      help='Filter out differentiated cells')
p$add_argument('--samples',  type="character",  nargs="+",    help='Samples')
p$add_argument('--which_clusters',  type="character",            default=NULL,      help='Should a rerun be forced if correlation file already exists')
p$add_argument('--n_clusters',  type="integer",            default=7,      help='Amount of gene clusters')
args <- p$parse_args(commandArgs(TRUE))

# ## START TEST ##
# args <- list()
# args$sce <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/metacells/SingleCellExperiment_metacells_nodiff.rds"
# args$atac_peak_matrix <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/metacells/PeakMatrix_summarized_experiment_metacells_nodiff.rds"
# args$peak2gene <- "/data/homes/louisc/Project_Babraham/ATAC/archR/PeakCalls/peaks2genes/peaks2genes_all.txt.gz"
# args$atac_peak_metadata <- "/data/homes/louisc/Project_Babraham/ATAC/archR/PeakCalls/peak_metadata.tsv.gz"
# args$distance <- 5e4
# args$max_pval <- 0.01
# args$min.cor <- c(0.25,0.5,0.75)
# args$force_rerun <- FALSE
# args$remove_clusters <- NULL
# # args$remove_clusters <- 6
# # args$remove_clusters <- c(2;3)
# # args$remove_clusters <- c(1,4,5,6)
# args$outdir <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/rna_vs_acc/metacells/gene_expr_vs_peak_acc/all"
# ## END TEST ##

#####################
## Define settings ##
#####################

print("Metacells")

# I/O
args$outfile <- paste0(args$outdir,args$which_clusters,"/cor_gene_expr_vs_peak_acc_metacells.txt.gz")
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
rna_metacells.sce <- readRDS(args$sce)
rna_metacells.sce

# Load ATAC SummarizedExperiment
atac_peakMatrix_metacells.se <- readRDS(args$atac_peak_matrix)
atac_peakMatrix_metacells.se

ncol(rna_metacells.sce)
sum(colnames(rna_metacells.sce)==colnames(atac_peakMatrix_metacells.se))

# Normalise ATAC data
assayNames(atac_peakMatrix_metacells.se) <- "counts"
assay(atac_peakMatrix_metacells.se,"logcounts") <- log(1e6*(sweep(assay(atac_peakMatrix_metacells.se),2,colSums(assay(atac_peakMatrix_metacells.se),na.rm=T),"/"))+1)
atac_peakMatrix_metacells.se

# Make sure that samples are consistent
celltypes <- intersect(colnames(rna_metacells.sce),colnames(atac_peakMatrix_metacells.se))

if (!is.null(args$remove_clusters)){
  for (i in args$remove_clusters){
    celltypes <- celltypes[!grepl(paste0("^",i),celltypes)]
  }
}
print(table(gsub("#.*$","",celltypes)))
rna_metacells.sce <- rna_metacells.sce[,celltypes]
atac_peakMatrix_metacells.se <- atac_peakMatrix_metacells.se[,celltypes]

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

to.plot.density_rna <- cbind(logcounts(rna_metacells.sce)%>% melt(),rep("RNA",ncol(rna_metacells.sce)*nrow(rna_metacells.sce)))         
to.plot.density_atac <- cbind(assay(atac_peakMatrix_metacells.se,"logcounts")[peak2gene.dt$peak,]%>% melt(),
                         rep("ATAC",ncol(assay(atac_peakMatrix_metacells.se,"logcounts")[peak2gene.dt$peak,])*nrow(atac_peakMatrix_metacells.se[peak2gene.dt$peak,])))
colnames(to.plot.density_rna) <- c("feature","sample","logcount","modality")
colnames(to.plot.density_atac) <- c("feature","sample","logcount","modality")
to.plot.density <- rbind(to.plot.density_rna,to.plot.density_atac)
to.plot.density$sample <- as.factor(to.plot.density$sample)
to.plot.density$modality <- as.factor(to.plot.density$modality)
print(head(to.plot.density))
print(summary(to.plot.density))

p <- ggplot(data=to.plot.density) +
  geom_density(mapping = aes(x=logcount, color=sample)) + 
  theme(legend.position = "none") +
  facet_wrap(~ modality)

pdf(file = sprintf("%s/density_plot.pdf", args$outdir),
    height = 6, width = 8)
print(p)
dev.off()

metacells_dos <- unique(gsub("#.+$","",colnames(logcounts(rna_metacells.sce))))
to.plot.density.av <- NULL
for (i in 1:length(metacells_dos)){
  to.plot.density.av <- rbind(to.plot.density.av,cbind(rownames(rna_metacells.sce),
                                                       rep(metacells_dos[i],nrow(rna_metacells.sce)),
                                                       rowMeans(logcounts(rna_metacells.sce)[,grepl(metacells_dos[i],colnames(logcounts(rna_metacells.sce)))],na.rm=T),                                                       
                                                       rep("RNA",nrow(rna_metacells.sce))))
  to.plot.density.av <- rbind(to.plot.density.av,cbind(rownames(assay(atac_peakMatrix_metacells.se,"logcounts")[peak2gene.dt$peak,]),
                                                       rep(metacells_dos[i],nrow(assay(atac_peakMatrix_metacells.se,"logcounts")[peak2gene.dt$peak,])),  
                                                       rowMeans(assay(atac_peakMatrix_metacells.se,"logcounts")[peak2gene.dt$peak,grepl(metacells_dos[i],colnames(assay(atac_peakMatrix_metacells.se,"logcounts")))],na.rm=T),                                                       
                                                       rep("ATAC",nrow(assay(atac_peakMatrix_metacells.se,"logcounts")[peak2gene.dt$peak,]))))                                              
}
to.plot.density.av <- data.frame(to.plot.density.av)
colnames(to.plot.density.av) <- c("feature","sample","logcount","modality")
to.plot.density.av$logcount <- as.numeric(to.plot.density.av$logcount)
to.plot.density.av$sample <- factor(to.plot.density.av$sample,levels=args$samples)
to.plot.density.av$modality <- as.factor(to.plot.density.av$modality)

print(head(to.plot.density.av))
print(summary(to.plot.density.av))

p <- ggplot(data=to.plot.density.av) +
  geom_density(mapping = aes(x=logcount, color=sample)) + 
  facet_wrap(~ modality)

pdf(file = sprintf("%s/density_plot_av.pdf", args$outdir),
    height = 6, width = 8)
print(p)
dev.off()

###########################################################
## Correlate peak accessibility with gene RNA expression ##
###########################################################

# TFs <- intersect(colnames(motifmatcher.se),rownames(rna_metacells.sce))# %>% head(n=5)
# TFs <- c("GATA1","TAL1")
genes <- intersect(rownames(rna_metacells.sce),peak2gene.dt$gene)

if (args$force_rerun | !file.exists(args$outfile)){
  cor.dt <- genes %>% map(function(i) {
    tmp <- peak2gene.dt[gene==i,c("gene","peak")]
    print(sprintf("%s (%d/%d)",i,match(i,genes),length(genes)))
    # calculate correlations
    corr_output <- psych::corr.test(
      x = t(logcounts(rna_metacells.sce[i,])), 
      y = t(assay(atac_peakMatrix_metacells.se[tmp$peak,],"logcounts")), 
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

sum((gene_expr_vs_peak_acc_filt.dt$pvalue<=args$max_pval) & (abs(gene_expr_vs_peak_acc_filt.dt$cor)>=args$min_cor))

###########################################################################################
## Calculate number of peaks that are correlated with RNA expression of the nearest gene ##
###########################################################################################

tmp <- gene_expr_vs_peak_acc_filt.dt %>% .[,sum(pvalue<=args$max_pval), by=c("peak")] %>% .[,mean(V1>=1)]

to.plot <- gene_expr_vs_peak_acc_filt.dt %>%
  .[,sum(pvalue<=args$max_pval & abs(cor)>=args$min_cor), by=c("peak","sign","peakType")] %>%
  .[,100*mean(V1>=1),by=c("sign","peakType")]

print(to.plot)

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

to.plot <- gene_expr_vs_peak_acc_filt.dt %>% .[pvalue<=args$max_pval]

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

to.plot <- gene_expr_vs_peak_acc_filt.dt %>% .[pvalue<=args$max_pval  & abs(cor)>=args$min_cor]

p <- gghistogram(to.plot, x="cor", y="..density..", bins=50, fill="gray70") +
  labs(x="Pearson correlation coefficient\n(ATAC peak vs gene expr.)", y="Density") +
  theme(
    axis.title.x = element_text(size=rel(0.9)),
    axis.title.y = element_text(size=rel(1.0)),
    legend.title = element_blank(),
    axis.text.x = element_text(size=rel(0.8), color="black"),
    axis.text.y = element_text(size=rel(0.8), color="black")
  )

pdf(file.path(args$outdir,paste0("histogram_sign_cor",args$min_cor,"gene_expr_vs_peak_acc_corr.pdf")), width=7, height=5)
print(p)
dev.off()

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
  sce <- readRDS(args$sce_pseudobulk)
  exp.values <- NULL
  for (i in 1:5){
    exp.values <- cbind(exp.values,rowMeans(assays(sce)[["logcounts"]][rownames(sce)%in%genes_neg_corr,sce$celltype==i]))
  }
  colnames(exp.values) <- paste0(rep("cluster",5),1:5)
  
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
  
  atac.se <- readRDS(args$atac_peak_matrix_pseudobulk)
  exp.values_atac <- NULL
  for (i in 1:5){
    exp.values_atac <- cbind(exp.values_atac,rowMeans(assays(atac.se)[["logcounts"]][rownames(atac.se)%in%exp.values.scaled_multi$peak,sce$celltype==i]))
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
  sce <- readRDS(args$sce_pseudobulk)
  exp.values <- NULL
  for (i in 1:5){
    exp.values <- cbind(exp.values,rowMeans(assays(sce)[["logcounts"]][rownames(sce)%in%genes_pos_corr,sce$celltype==i]))
  }
  colnames(exp.values) <- paste0(rep("cluster",5),1:5)
  
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
  
  atac.se <- readRDS(args$atac_peak_matrix_pseudobulk)
  exp.values_atac <- NULL
  for (i in 1:5){
    exp.values_atac <- cbind(exp.values_atac,rowMeans(assays(atac.se)[["logcounts"]][rownames(atac.se)%in%exp.values.scaled_multi$peak,sce$celltype==i]))
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


##########################################################
## Expression profile of significantly correlated genes ##
##########################################################

unlink(paste0(args$outdir,"/gene_plots"), recursive = T, force = T)
dir.create(paste0(args$outdir,"/gene_plots"))

genelist <- c(opts$naieve_pluri,opts$general_pluri,opts$common_post_impl,opts$housekeeping)
for (i in genelist){
  peaklist <- cor.dt$peak[(cor.dt$gene==i) & (cor.dt$pvalue<args$max_pval)]
  
  if (length(peaklist)==0){
    pdf(file.path(args$outdir,sprintf("/gene_plots/%s_no_sign_corr.pdf",i)), width = 8, height = 5)
    dev.off()
    next()
  }
  
  to.plot <- NULL
  for (j in peaklist){
    to.plot <- rbind(to.plot,data.table(
      atac = assay(atac_peakMatrix_metacells.se[j,],"logcounts")[1,],
      rna = logcounts(rna_metacells.sce[i,])[1,],
      peak = rep(j,length(logcounts(rna_metacells.sce[i,])[1,])),
      sample = gsub("#.*$","",colnames(atac_peakMatrix_metacells.se)))
    )
  }
  
  if(args$sort_samples){
    to.plot$sample <- factor(to.plot$sample,levels=args$samples)
  }
  
  p <- ggplot(to.plot, aes(x=rna, y=atac)) +
    geom_point(aes(fill=sample), color="black", size=2, shape=21) +
    geom_smooth(method=lm, se=FALSE) +
    stat_cor(method = "pearson") +
    scale_fill_manual(values=opts$color_scheme) +
    facet_wrap(~peak, scales="fixed") +
    labs(x="RNA expression", y="Peak accessibility", title=sprintf("%s expression vs accessibility",i)) +
    theme_classic() + 
    theme(
      plot.title = element_text(hjust=0.5, size=rel(0.8)),
      axis.text = element_text(color="black"),
    )
  
  pdf(file.path(args$outdir,sprintf("/gene_plots/%s_rna_vs_acc_metacells.pdf",i)), width = 8, height = 5)
  print(p)
  dev.off()
}

