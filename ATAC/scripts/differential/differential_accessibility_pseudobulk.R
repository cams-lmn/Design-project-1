###############################################
##                                           ##
##  Differential_accessibility_pseudobulk.R  ##
##                                           ##
###############################################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--matrix',        type="character",     help='')
p$add_argument('--cluster_ids',            type="character",     nargs='+',             help='Which cluster ids should be compared')
p$add_argument('--seed1',            type="integer",     default=42,             help='Random seed 1')
p$add_argument('--peak_metadata',  type="character",                help='Peak metadata')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
p$add_argument('--fdr_cutoff',            type="double",     default=0.01,             help='FDR cutoff')
p$add_argument('--logfc_cutoff',            type="double",     default=1.00,             help='LogFC cutoff')
p$add_argument("--rm_prev_res",  type="logical",default=T, help="Remove previous differential results?")
p$add_argument('--min_cdr',            type="double",     default=0.30,             help='LogFC cutoff')
p$add_argument('--group_by',    type="character",    help='Metadata column to group cells by')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

args$atac_matrix_file <- sprintf("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/ATAC/archR/pseudobulk/%s/%s/pseudobulk_with_replicates.rds",
                                 args$group_by,args$matrix)

args$groups <- combn(args$cluster_ids,2)
print(args$groups)

args$outdir <- sprintf("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/ATAC/archR/pseudobulk/%s/%s/DE_res", args$group_by,args$matrix)

###############
## Load data ##
###############

print(sprintf("Fetching ATAC matrix: '%s'...",args$atac_matrix_file))

atac.se <- readRDS(args$atac_matrix_file)

# # parse
# if (!"celltype" %in% colnames(colData(atac.se))) {
#   atac.se$celltype <- colnames(atac.se) %>% strsplit("-") %>% map_chr(1)
# }
# atac.se <- atac.se[,atac.se$celltype%in%opts$groups]
# atac.se$celltype <- factor(atac.se$celltype, levels=opts$groups)

# Normalise with log2 counts (for consistentcy with edgeR)
assayNames(atac.se)[1] <- "counts"
logcounts(atac.se) <- log2(1e6*(sweep(counts(atac.se),2,colSums(counts(atac.se)),"/"))+1)

###########################
## Differential analysis ##
###########################

for (i in 1:ncol(args$groups)){
  args$groupA <- as.character(args$groups[1,i])
  args$groupB <- as.character(args$groups[2,i])
  args$outfile <-  sprintf("%s/res_cluster%svsclusters%s.txt",args$outdir,args$groupA,args$groupB)

  if(args$rm_prev_res & i==1){
    unlink(dirname(args$outfile), recursive = T, force = T)
  }
  
  dir.create(dirname(args$outfile), showWarnings=F, recursive = T)
  
  #####################
  ## Define settings ##
  #####################
  
  # Define groups
  opts$groups <- c(args$groupA,args$groupB)

  # stupid stuff but otherwise the snakemake pipeline doesn't work
  if (args$groupA==args$groupB) {
    out <- data.table(feature=NA, logFC=NA, padj_fdr=NA, mean_groupA=NA, mean_groupB=NA)
    fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
    warning("groupA and groupB are the same, saving an empty file...")
    quit(status=0)
  }
  
  # check that we have pseudobulk replicates for both cell types
  if (any(!opts$groups%in%unique(atac.se$celltype))) {
    out <- data.table(feature=NA, logFC=NA, padj_fdr=NA, mean_groupA=NA, mean_groupB=NA)
    fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
    warning("groups not found, saving an empty file...")
    quit(status=0)
  }
  
  ##################################
  ## calculate feature statistics ##
  ##################################
  
  acc.dt <- data.table(
    feature = rownames(atac.se),
    detection_rate_groupA = rowMeans(assay(atac.se,"counts")[,atac.se$celltype==opts$groups[1]]>0) %>% round(2),
    detection_rate_groupB = rowMeans(assay(atac.se,"counts")[,atac.se$celltype==opts$groups[2]]>0) %>% round(2),
    mean_groupA = rowMeans(logcounts(atac.se[,atac.se$celltype==args$groupA])) %>% round(2),
    mean_groupB = rowMeans(logcounts(atac.se[,atac.se$celltype==args$groupB])) %>% round(2)
  )
  
  if(i==1){
    DE_genes <- matrix(rep(F,nrow(atac.se)*ncol(args$groups)),ncol=ncol(args$groups))
    rownames(DE_genes) <- sort(rownames(atac.se))
    colnames(DE_genes) <- paste0(args$groups[1,],rep("vs",ncol(args$groups)),args$groups[2,])
  }
  
  #####################################
  ## Differential testing with edgeR ##
  #####################################
  
  # Consider only features that have a minimum detection rate
  features.to.use <- acc.dt[detection_rate_groupA>=args$min_cdr | detection_rate_groupB>=args$min_cdr,feature]
  
  # Convert SCE to DGEList
  atac_metacells.dge <- DGEList(assay(atac.se[features.to.use,],"counts"))
  
  # Define design matrix (with intercept)
  design <- model.matrix(~atac.se$celltype)
  
  # Estimate dispersions
  atac_metacells.dge <- estimateDisp(atac_metacells.dge,design)
  
  # Fit GLM
  fit <- glmQLFit(atac_metacells.dge,design)
  
  # Likelihood ratio test
  lrt <- glmQLFTest(fit)
  
  # Construct output data.frame
  out <- topTags(lrt, n=nrow(lrt))$table %>% as.data.table(keep.rownames=T) %>%
    setnames(c("feature","logFC","logCPM","LR","p.value","padj_fdr")) %>%
    .[,c("logCPM","LR","p.value"):=NULL] %>%
    merge(acc.dt[,c("feature","mean_groupA","mean_groupB")], by="feature", all.y=TRUE) %>%
    .[,c("padj_fdr","logFC"):=list(signif(padj_fdr,digits=3),round(logFC,3))] %>%
    .[is.na(logFC),c("logFC","padj_fdr"):=list(0,1)] %>%
    setorder(padj_fdr, na.last=T)
  
  mask_out <- (out$padj_fdr<args$fdr_cutoff) & (abs(out$logFC)>args$logfc_cutoff)
  DE_genes[rownames(DE_genes)%in%out$feature[mask_out],i] <- T
  
  ##################
  ## Save results ##
  ##################
  
  fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
}

print(dim(DE_genes))
print(sum(rowSums(DE_genes)>0))
print(head(DE_genes))
saveRDS(DE_genes,sprintf("%s/DE_genes.rds",args$outdir))

###################
# Gene clustering #
###################

exp.values <- NULL
for (i in 1:6){
  exp.values <- cbind(exp.values,rowMeans(assay(atac.se,"logcounts")[rownames(atac.se)%in%rownames(DE_genes)[rowSums(DE_genes)>0],atac.se$celltype==paste0("C",i)]))
}
colnames(exp.values) <- args$cluster_ids

z_score <- function(x){
  (x - mean(x)) / sd(x)
}

t(apply(exp.values, 1, z_score)) %>%
  as_tibble(rownames = "Gene") %>%
  column_to_rownames(var="Gene") -> exp.values.scaled


# Elbow plot to estimate the number of clusters
wss <- (nrow(exp.values.scaled)-1)*sum(apply(exp.values.scaled,2,var))
for (i in 2:50) wss[i] <- sum(kmeans(exp.values.scaled,
                                     centers=i)$withinss)

pdf(file = sprintf("%s/elbow.plot.pdf",args$outdir),
    height = 6, width = 8)
plot(1:50, wss,
     type="b",
     xlab="Number of Clusters",
     ylab="Within groups sum of squares")
dev.off()


# Then we perform clustering and select various numbers of clusters
set.seed(args$seed1)

dist_matrix <- dist(exp.values.scaled, method = "euclidean")

clustering.genes <- hclust(dist_matrix, method = "ward.D2")

saveRDS(clustering.genes,sprintf("%s/clustering.genes.rds",args$outdir))

### This is to visulalise the tree
# plot(clustering.genes)

### Then we look at various options, how many clusters to extract

count.gaps <- function(x,clusters){
  clusters %>%
    filter(cluster<=x) %>%
    nrow -> y
  return(y)
}

for (i in 5:10){
  cutree(clustering.genes, k=i) %>%
    as_tibble(rownames = "Gene") %>%
    dplyr::rename(cluster = value) -> clusters.n
  
  lapply(c(1:i), count.gaps, clusters.n) %>%
    do.call("rbind",.) %>%
    as.vector -> clusters.gaps.n
  
  cols3 <- rev(brewer.pal(10, "Spectral"))
  
  annotation.colors <- list(cluster = c("1"=cols3[1], "2"=cols3[2], "3"=cols3[3], "4"=cols3[4], "5"=cols3[5],
                                        "6"=cols3[6], "7"=cols3[7], "8"=cols3[8], "9"=cols3[9], "10"=cols3[10]
  ))
  annotation.colors <- annotation.colors[1:i]
  
  pdf(file = sprintf("%s/heatmap_%sgeneclusters.pdf",args$outdir,i),
      height = 6, width = 8)
  exp.values.scaled %>%
    rownames_to_column(var="Gene") %>%
    full_join(clusters.n) %>%
    arrange(cluster) %>%
    column_to_rownames(var = "Gene") %>%
    select(-cluster) %>%
    pheatmap(cluster_rows = F,
             cluster_cols = F,
             show_rownames = F,
             gaps_row = clusters.gaps.n,
             annotation_row = (clusters.n %>% column_to_rownames(var="Gene") %>% arrange(cluster)),
             annotation_colors = annotation.colors,
             color = colorRampPalette(cols3)(100),breaks=seq(from=-2,to=2,length.out=101))
  dev.off()
}

#############
# DEG table #
#############

DEG_table <- NULL
c_names <- NULL

letters <- c("A","B","C","D","E","F","G","H","I","J")

if (args$matrix=="PeakMatrix") {
  order_gene_clusters <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
} else if (grepl("GeneScoreMatrix_TSS",args$matrix)) {
  order_gene_clusters <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
} else  if (grepl("GeneScoreMatrix_distal",args$matrix)){
  order_gene_clusters <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
}

pairw_comp <- paste0(combn(args$cluster_ids,2)[1,],
                     rep("vs",ncol(combn(args$cluster_ids,2))),
                     combn(args$cluster_ids,2)[2,])

# Average expression
for (i in args$cluster_ids){
  DEG_table <- cbind(DEG_table,rowMeans(assay(atac.se,"logcounts")[,atac.se$celltype==i]))
  c_names <- c(c_names,paste0("log_av_",i))
}

DEG_table_sorted <- DEG_table[sort(rownames(DEG_table),index.return=T)$ix,]

# DE results
for (i in 1:length(pairw_comp)){
  res_i <- read.table(sprintf("%s/res_cluster%svsclusters%s.txt",args$outdir,substr(pairw_comp[i],1,2),substr(pairw_comp[i],5,6)),
                      header=T, sep="\t")
  res_i_sorted <- res_i[sort(res_i$feature,index.return=T)$ix,]
  DEG_table_sorted <- cbind(DEG_table_sorted,res_i_sorted[,2:3],
                            (abs(res_i_sorted$logFC)>1) & (res_i_sorted$padj_fdr<0.01))
  c_names <- c(c_names,paste0(c("logFC_","FDR","Sign"),rep(pairw_comp[i],3)))
}

if (args$matrix=="PeakMatrix") {
  n_clusters <- 9
} else if (grepl("GeneScoreMatrix",args$matrix)) {
  n_clusters <- 9
} 

# Gene clustering
cutree(clustering.genes, k=n_clusters) %>%
  as_tibble(rownames = "Gene") %>%
  dplyr::rename(cluster = value) -> clusters.n

saveRDS(clusters.n,sprintf("%s/clusters.%s.rds",args$outdir,n_clusters))

clusters.n_sorted <- clusters.n[sort(clusters.n$Gene,index.return=T)$ix,]
clusters.n_sorted$cluster <- letters[order_gene_clusters][clusters.n_sorted$cluster]

DEG_table_sorted <- cbind(DEG_table_sorted,rep(NA,nrow(DEG_table_sorted)))
DEG_table_sorted[rownames(DEG_table_sorted)%in%clusters.n_sorted$Gene,ncol(DEG_table_sorted)] <- clusters.n_sorted$cluster

c_names <- c(c_names,"gene_cluster")

# Cluster spec
for (i in paste0(rep("C",7),1:7)){
  cols_logFC <- which(grepl(i,c_names) & grepl("vs",c_names) & grepl("logFC",c_names))
  cols_FDR <-  which(grepl(i,c_names) & grepl("vs",c_names) & grepl("FDR",c_names))
  sign <- c(1,-1)[c(substr(gsub("[A-Z,a-z]+_","",c_names[cols_logFC]),1,1)==i)+1]
  
  DEG_table_sorted <- cbind(DEG_table_sorted,rowSums((t(t(DEG_table_sorted[,cols_logFC])*sign)>1) & (DEG_table_sorted[,cols_FDR]<0.01))>=3)
  c_names <- c(c_names,paste0("clust",i,"_spec"))
}

# Write table
DEG_table_sorted <- cbind(rownames(DEG_table_sorted),DEG_table_sorted)
colnames(DEG_table_sorted) <- c("Gene",c_names)
write.table(DEG_table_sorted,file=sprintf("%s/DEG_overview.txt",args$outdir),
            col.names=T,row.names = F,sep="\t",quote=F)

####################
# Final clustering #
####################

clusters.n_ordered <- clusters.n
clusters.n_ordered$cluster <- as.numeric(clusters.n_ordered$cluster)
clusters.n_ordered$cluster_new <- order_gene_clusters[as.numeric(clusters.n_ordered$cluster)]
table(clusters.n_ordered$cluster,clusters.n_ordered$cluster_new)
clusters.n_ordered$cluster <- clusters.n_ordered$cluster_new
clusters.n_ordered <- clusters.n_ordered[,-3]

lapply(c(1:n_clusters), count.gaps, clusters.n_ordered) %>%
  do.call("rbind",.) %>%
  as.vector -> clusters.gaps.n

cols3 <- rev(brewer.pal(10, "Spectral"))

annotation.colors <- list(cluster = c("1"=cols3[1], "2"=cols3[2], "3"=cols3[3], "4"=cols3[4], "5"=cols3[5],
                                      "6"=cols3[6], "7"=cols3[7], "8"=cols3[8], "9"=cols3[9], "10"=cols3[10]
))
annotation.colors <- annotation.colors[1:n_clusters]

head(exp.values.scaled)

# pdf 
pdf(file = sprintf("%s/heatmap_%sgeneclusters_final.pdf",args$outdir,n_clusters),
    height = 6, width = 8)
exp.values.scaled %>%
  rownames_to_column(var="Gene") %>%
  full_join(clusters.n_ordered) %>%
  arrange(cluster) %>%
  column_to_rownames(var = "Gene") %>%
  select(-cluster) %>%
  pheatmap(cluster_rows = F,
           cluster_cols = F,
           show_rownames = F,
           gaps_row = clusters.gaps.n,
           annotation_row = (clusters.n_ordered %>% column_to_rownames(var="Gene") %>% arrange(cluster)),
           annotation_colors = annotation.colors,
           color = colorRampPalette(cols3)(100),breaks=seq(from=-2,to=2,length.out=101))
dev.off()

# pdf - nolabs
pdf(file = sprintf("%s/heatmap_%sgeneclusters_final_nolabs.pdf",args$outdir,n_clusters),
    height = 6, width = 8)
exp.values.scaled %>%
  rownames_to_column(var="Gene") %>%
  full_join(clusters.n_ordered) %>%
  arrange(cluster) %>%
  column_to_rownames(var = "Gene") %>%
  select(-cluster) %>%
  pheatmap(cluster_rows = F,
           cluster_cols = F,
           show_rownames = F,           
           show_colnames = F,
           gaps_row = clusters.gaps.n,
           legend=FALSE,
           annotation_legend = FALSE,
           annotation_names_row	= FALSE,
           annotation_names_col = FALSE,
           color = colorRampPalette(cols3)(100),breaks=seq(from=-2,to=2,length.out=101))
dev.off()

#################################
# Violin plots per gene cluster #
#################################

exp.values.scaled <- as.data.table(exp.values.scaled)
exp.values.scaled.melted <- melt(exp.values.scaled, variable.name = "variable", value.name = "value")

exp.values.scaled.melted <- cbind(exp.values.scaled.melted,
                                  rep(c("A","B","C","D","E","F","G","H","I")[clusters.n_ordered$cluster],6))
colnames(exp.values.scaled.melted)[3] <- "gene_cluster"
head(exp.values.scaled.melted)

pdf(file = sprintf("%s/violin_%sgeneclusters_final.pdf",args$outdir,n_clusters),
    height = 6, width = 8)
ggplot(exp.values.scaled.melted, aes(x=variable, y=value, fill=variable)) + 
  geom_violin(trim=FALSE) + 
  facet_wrap(~gene_cluster, scales="fixed") +
  ggplot_theme_NoAxes()
dev.off()


# ##########################################
# ## Peak distribution in genomic context ##
# ##########################################

#  if (grepl("PeakMatrix",args$matrix)){

#   print("Peak distribution")

#   peak_metadata.dt <- fread(args$peak_metadata) %>%
#     .[,peak:=sprintf("%s:%s-%s",chr,start,end)] %>%
#     .[,c("peak","peakType","distToGeneStart","nearestGene")]
  
#   colnames(peak_metadata.dt)[1] <- "Gene"

#   print(dim(DEG_table_sorted))
#   DEG_table_sorted <- merge(DEG_table_sorted,peak_metadata.dt,all.x=TRUE)
#   print(dim(DEG_table_sorted))
#   # print(colnames(DEG_table_sorted))

#   to.plot <- table(DEG_table_sorted$gene_cluster,DEG_table_sorted$peakType)
#   to.plot <- melt(to.plot)
#   colnames(to.plot) <- c("gene_cluster","peakType","V1") 
#   print(head(to.plot))

#   for (j in unique(to.plot$gene_cluster)){
#     to.plot$V1[to.plot$gene_cluster==j] <- to.plot$V1[to.plot$gene_cluster==j]/sum(to.plot$V1[to.plot$gene_cluster==j])
#   }

#   to.plot$V1 <- to.plot$V1*100

#   p <- ggbarplot(to.plot, x="gene_cluster", fill="peakType", y="V1", position=position_dodge(width = 0.75)) +
#     labs(x="", y="% of peaks present in genomic context") +
#     # scale_fill_brewer(palette="Dark2") +
#     # guides(color=F) + 
#     theme(
#       axis.title = element_text(size=rel(0.80)),
#       legend.title = element_blank(),
#       axis.text.x = element_text(size=rel(0.80), color="black"),
#       axis.text.y = element_text(size=rel(0.80), color="black")
#     )

#   pdf(file.path(args$outdir,"barplots_number_DAP_vs_genomic_context.pdf"), width=7, height=5)
#   print(p)
#   dev.off()
#  }


###################
# Peak annotation #
###################

if (args$matrix=="PeakMatrix"){
  print("Peak distribution per cluster")

  # Load peak metadata
  peakSet.dt <- fread(args$peak_metadata) %>%
    .[,chr:=as.factor(sub("chr","",chr))] %>%
    .[,peak:=sprintf("chr%s:%s-%s",chr,start,end)] 

  print(head(peakSet.dt))

  print("Barplot peak type attribution - per cluster")
  if(!dir.exists(sprintf("%s/barplots_per_cluster",args$outdir))){dir.create(sprintf("%s/barplots_per_cluster",args$outdir))}

  for (i in unique(clusters.n_ordered$cluster)){
    print(i)

    peaks_i <- clusters.n_ordered$Gene[clusters.n_ordered$cluster==i]
    # Barplot of the fraction of cells that pass QC for each sample
    peakSet.dt_i <- peakSet.dt[peakSet.dt$peak%in%peaks_i]
    to.plot <- data.frame(peakType=names(table(peakSet.dt_i$peakType)),freq=as.numeric(table(peakSet.dt_i$peakType)/nrow(peakSet.dt_i)),stringsAsFactors=F)
    print(to.plot)

    p <- ggbarplot(to.plot, x="peakType", y="freq",fill="peakType",label=as.character(round(to.plot$freq,digits=3)),lab.vjust=-4) +
      coord_cartesian(ylim=c(0,1)) +
      theme(
        legend.position = "none",
        axis.text.y = element_text(colour="black",size=rel(0.8)),
        axis.text.x = element_text(colour="black",size=rel(0.65), angle=20, hjust=1, vjust=1),
      ) + 
      scale_fill_manual(values=c("#E6E7E9", "#BCBDBF", "#828387", "#3E3E40")) 

    output_plot(p, sprintf("%s/barplots_per_cluster/barplot_cluster%s_peakTypeAnnotation",args$outdir,i), width=6, height=5)
  }
}