###############################
##                           ##
##  Differential_analysis.R  ##
##                           ##
###############################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA_ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--fdr_cutoff',            type="double",     default=0.01,             help='FDR cutoff')
p$add_argument('--logfc_cutoff',            type="double",     default=1.00,             help='LogFC cutoff')
p$add_argument("--rm_prev_res",  type="logical",default=T, help="Remove previous differential results?")
p$add_argument('--matrix',          type="character",  nargs="+",    help='Matrix to use')
p$add_argument('--min_expr',            type="double",     default=4,             help='LogFC cutoff')
p$add_argument('--min_cdr',            type="double",     default=0.30,             help='LogFC cutoff')
p$add_argument('--cluster_ids',            type="character",     nargs='+',             help='Which cluster ids should be compared')
p$add_argument('--seed',            type="integer",     default=42,             help='Random seed')
p$add_argument('--atac_peak_metadata',  type="character",        default=NULL,       help='Peak metadata')
args <- p$parse_args(commandArgs(TRUE))


#####################
## Define settings ##
#####################
args$matrix <- "RNA"
args$sce <- sprintf("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA_ATAC/pseudobulk/cluster/%s/pseudobulk_with_replicates.rds",
                    args$matrix)

args$groups <- combn(args$cluster_ids,2)

args$QC_dir <- sprintf("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA_ATAC/pseudobulk/cluster/%s/QC",
                           args$matrix)

dir.create(args$QC_dir, showWarnings=F)

if(args$rm_prev_res){
  args$outdir <- sprintf("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA_ATAC/pseudobulk/cluster/%s/DE_res",
                    args$matrix)
  unlink(args$outdir, recursive = T, force = T)
}
  
dir.create(args$outdir, showWarnings = F)

###################################
## Load expression/accessibility ##
###################################

# Load SingleCellExperiment object
sce <- readRDS(args$sce)

###################
## Density plots ##
###################

sce_logcount <- assay(sce,"logcounts")
sce_logcount_filt <- sce_logcount[rowMeans(sce_logcount)>args$min_expr,]

# Ensure sce_logcount_filt is a data.table
sce_logcount_filt_dt <- as.data.table(sce_logcount_filt, keep.rownames = "gene")

# Melt the data.table
to.plot <- melt(sce_logcount_filt_dt, id.vars = "gene", variable.name = "pseudorep", value.name = "logcount")

# to.plot <- to.plot[!(grepl("6_",to.plot$pseudorep)),]
to.plot$pseudorep <- as.factor(as.character(to.plot$pseudorep))

to.plot$cluster <- gsub("_.+$","",to.plot$pseudorep)
to.plot$rep <- gsub("^.+_","",to.plot$pseudorep)
print(head(to.plot))


p <- ggplot(data=to.plot) +
  geom_density(mapping = aes(x=logcount, color=rep)) +
  facet_wrap(~ cluster)

pdf(file = sprintf("%s/%s_density_pseudorep_percluster.pdf", args$outdir, args$matrix),
    height = 6, width = 8)
print(p)
dev.off()

p <- ggplot(data=to.plot) +
  geom_density(mapping = aes(x=logcount, color=pseudorep))

pdf(file = sprintf("%s/%s_density_pseudorep.pdf", args$outdir, args$matrix),
    height = 6, width = 8)
print(p)
dev.off()
  
############################
## Option mutliple groups ##
############################

for (i in 1:ncol(args$groups)){
  args$groupA <- as.character(args$groups[1,i])
  args$groupB <- as.character(args$groups[2,i])
  args$outfile <-  sprintf("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA_ATAC/pseudobulk/cluster/%s/DE_res/res_cluster%svsclusters%s.txt",
                           args$matrix,args$groupA,args$groupB)

  print(sprintf("%s: Differential analysis for cluster %s versus cluster %s",args$matrix,args$groupA,args$groupB))
  
  #####################
  ## Define settings ##
  #####################
  
  # Define groups
  opts$groups <- c(args$groupA,args$groupB)
  
  # otherwise the snakemake pipeline doesn't work
  if (args$groupA==args$groupB) {
    out <- data.table(feature=NA, logFC=NA, padj_fdr=NA)
    fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
    warning("groupA and groupB are the same, saving an empty file...")
    quit(status=0)
  }

  ################
  ## Subset sce ##
  ################

  sce_groups <- sce[,sce$celltype %in% opts$groups]
  
  sce_groups$celltype <- factor(sce_groups$celltype, levels=opts$groups)
  table(sce_groups$celltype)
    
  ###########################
  ## Differential_analysis ##
  ###########################

  if(i==1){
    DE_genes <- matrix(rep(F,nrow(sce)*ncol(args$groups)),ncol=ncol(args$groups))
    rownames(DE_genes) <- sort(rownames(sce))
    colnames(DE_genes) <- paste0(args$groups[1,],rep("vs",ncol(args$groups)),args$groups[2,])
  }
  
  if (grepl("RNA",args$matrix)){
    
    ## Calculate average expression levels
    #########################################
    
    expr.dt <- data.table(
      gene = rownames(sce_groups),
      mean_groupA = rowMeans(logcounts(sce_groups[,sce_groups$celltype==args$groupA])) %>% round(2),
      mean_groupB = rowMeans(logcounts(sce_groups[,sce_groups$celltype==args$groupB])) %>% round(2)
    )
    
    ## Feature selection 
    #######################
    
    features.to.use <- expr.dt[mean_groupA>=args$min_expr | mean_groupB>=args$min_expr,gene]
    
  } else {
    
    ## Calculate average expression levels
    #########################################
    
    acc.dt <- data.table(
      feature = rownames(sce_groups),
      detection_rate_groupA = rowMeans(assay(sce_groups,"counts")[,sce_groups$celltype==opts$groups[1]]>0) %>% round(2),
      detection_rate_groupB = rowMeans(assay(sce_groups,"counts")[,sce_groups$celltype==opts$groups[2]]>0) %>% round(2),
      mean_groupA = rowMeans(logcounts(sce_groups[,sce_groups$celltype==args$groupA])) %>% round(2),
      mean_groupB = rowMeans(logcounts(sce_groups[,sce_groups$celltype==args$groupB])) %>% round(2))
    
    ## Feature selection 
    #######################
    
    features.to.use <- acc.dt[detection_rate_groupA>=args$min_cdr | detection_rate_groupB>=args$min_cdr,feature]
    
  }
  
  ################################################
  ## Differential expression testing with edgeR ##
  ################################################
  
  # Convert SCE to DGEList
  sce_edger <- scran::convertTo(sce_groups[features.to.use,], type="edgeR")
  
  # Define design matrix (with intercept)
  design <- model.matrix(~sce_groups$celltype)
  
  # Estimate dispersions
  sce_edger  <- estimateDisp(sce_edger,design)
  
  # Fit GLM
  fit <- glmQLFit(sce_edger,design)
  
  # Likelihood ratio test
  lrt <- glmQLFTest(fit)
  
  # Construct output data.frame
  if (grepl("RNA",args$matrix)){
    out <- topTags(lrt, n=nrow(lrt))$table %>% as.data.table(keep.rownames=T) %>%
      setnames(c("gene","logFC","logCPM","LR","p.value","padj_fdr")) %>%
      .[,c("logCPM","LR","p.value"):=NULL] %>%
      .[,c("padj_fdr","logFC"):=list(signif(padj_fdr,digits=3), round(logFC,3))] %>%
      merge(expr.dt, by="gene", all.y=TRUE) %>%
      setorder(padj_fdr, na.last=T)
    
    mask_out <- (out$padj_fdr<args$fdr_cutoff) & (abs(out$logFC)>args$logfc_cutoff)
    DE_genes[rownames(DE_genes)%in%out$gene[mask_out],i] <- T
  } else {
    out <- topTags(lrt, n=nrow(lrt))$table %>% as.data.table(keep.rownames=T) %>%
      setnames(c("feature","logFC","logCPM","LR","p.value","padj_fdr")) %>%
      .[,c("logCPM","LR","p.value"):=NULL] %>%
      merge(acc.dt[,c("feature","mean_groupA","mean_groupB")], by="feature", all.y=TRUE) %>%
      .[,c("padj_fdr","logFC"):=list(signif(padj_fdr,digits=3),round(logFC,3))] %>%
      .[is.na(logFC),c("logFC","padj_fdr"):=list(0,1)] %>%
      setorder(padj_fdr, na.last=T)
    
    mask_out <- (out$padj_fdr<args$fdr_cutoff) & (abs(out$logFC)>args$logfc_cutoff)
    DE_genes[rownames(DE_genes)%in%out$feature[mask_out],i] <- T
  }
  print(head(out))
  
  ##################
  ## Save results ##
  ##################
  
  fwrite(out, args$outfile, sep="\t", na="NA", quote=F)

  ########
  ## QC ##
  ########

  out_all <- topTags(lrt, n=nrow(lrt))$table %>% as.data.table(keep.rownames=T) %>%
      setnames(c("gene","logFC","logCPM","LR","p.value","padj_fdr")) 
  print(head(out_all))

  # BCV plot
  jpeg(sprintf("%s/BCVplot_cluster%svsclusters%s.jpg",args$QC_dir,args$groupA,args$groupB), type = "cairo")
  plotBCV(sce_edger) 
  dev.off()

  ## MA plot
  jpeg(sprintf("%s/MAplot_cluster%svsclusters%s.jpg",args$QC_dir,args$groupA,args$groupB), type = "cairo")
  with(out_all,plot(logCPM,logFC,pch=16,cex=0.2)) 
  with(out_all,points(logCPM[padj_fdr<0.05],logFC[padj_fdr<0.05],pch=16,col="red",cex=0.6)) 
  abline(0,0) 
  dev.off()

  ## Pvalue distribution
  jpeg(sprintf("%s/Pval_dist_cluster%svsclusters%s.jpg",args$QC_dir,args$groupA,args$groupB), type = "cairo")
  hist(out_all$p.value)
  dev.off()
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
  exp.values <- cbind(exp.values,rowMeans(assays(sce)[["logcounts"]][rownames(sce)%in%rownames(DE_genes)[rowSums(DE_genes)>0],sce$celltype==i]))
}
colnames(exp.values) <- paste0(rep("cluster",6),1:6)

t(apply(exp.values, 1, z_score)) %>%
  as_tibble(rownames = "Gene") %>%
  column_to_rownames(var="Gene") -> exp.values.scaled

# Elbow plot to estimate the number of clusters
wss <- (nrow(exp.values.scaled)-1)*sum(apply(exp.values.scaled,2,var))
for (i in 2:50) wss[i] <- sum(kmeans(exp.values.scaled,
                                     centers=i, iter.max=20)$withinss)

pdf(file = sprintf("%s/elbow.plot.pdf", args$outdir),
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

saveRDS(clustering.genes,sprintf("%s/clustering.genes.rds", 
                                 args$outdir))

for (i in 5:10){
  cutree(clustering.genes, k=i) %>% 
    as_tibble(rownames = "Gene") %>%
    dplyr::rename(cluster = value) -> clusters.n
  
  count.gaps <- function(x,clusters){
    clusters %>%
      filter(cluster<=x) %>%
      nrow -> y
    return(y)
  }
   
  lapply(c(1:10), count.gaps, clusters.n) %>%
    do.call("rbind",.) %>% 
    as.vector -> clusters.gaps.n
  
  
  cols3 <- rev(brewer.pal(10, "Spectral"))
  
  annotation.colors <- list(cluster = c("1"=cols3[1], "2"=cols3[2], "3"=cols3[3], "4"=cols3[4], "5"=cols3[5],
                                        "6"=cols3[6], "7"=cols3[7], "8"=cols3[8], "9"=cols3[9], "10"=cols3[10])[1:i])
  annotation_row <- clusters.n %>% 
    column_to_rownames(var = "Gene") %>%
    arrange(cluster)
  annotation_row$cluster <- as.factor(annotation_row$cluster)
  pdf(file = sprintf("%s/heatmap_%sgeneclusters.pdf", args$outdir, i), 
      height = 6, width = 8)
  exp.values.scaled %>%
    rownames_to_column(var="Gene") %>% 
    full_join(clusters.n, by = "Gene") %>%
    arrange(cluster) %>%
    column_to_rownames(var = "Gene") %>%
    select(-cluster) %>%
    pheatmap(cluster_rows = F,
             cluster_cols = F,
             show_rownames = F,
             gaps_row = clusters.gaps.n,
             annotation_row = annotation_row,
             annotation_colors = annotation.colors,
             color = colorRampPalette(cols3)(100),
             breaks=seq(from=-2,to=2,length.out=101))
  dev.off()
}

#############
# DEG table #
#############

DEG_table <- NULL
c_names <- NULL
letters <- c("A","B","C","D","E","F","G","H","I","J")

if (args$matrix=="RNA"){
  order_gene_clusters <- c(1,2,3,4,5,6,7,8,9)
  n_clusters <- 9
} else if(grepl("Peak",args$matrix)){
  order_gene_clusters <- c(1,2,3,4,5,6,7,8,9)
  n_clusters <- 9
} else if(grepl("TSS",args$matrix)){ 
  order_gene_clusters <- c(1,2,3,4,5,6,7,8,9)
  n_clusters <- 9
} else if(grepl("distal",args$matrix)){ 
  order_gene_clusters <- c(1,2,3,4,5,6,7,8,9)
  n_clusters <- 9
}
pairw_comp <- paste0(combn(1:6,2)[1,],rep("vs",ncol(combn(1:6,2))),combn(1:6,2)[2,])

# Average expression
for (i in 1:6){
  DEG_table <- cbind(DEG_table,rowMeans(assays(sce)[["logcounts"]][,sce$celltype==i]))
  c_names <- c(c_names,paste0("log_av_",i))
}

DEG_table_sorted <- DEG_table[sort(rownames(DEG_table),index.return=T)$ix,]

# DE results
for (i in 1:length(pairw_comp)){
  res_i <- read.table(sprintf("%s/res_cluster%svsclusters%s.txt",args$outdir, 
                              substr(pairw_comp[i],1,1), substr(pairw_comp[i],4,4)),
                      header=T, sep="\t")
  res_i_sorted <- res_i[sort(res_i[,1],index.return=T)$ix,]
  DEG_table_sorted <- cbind(DEG_table_sorted,res_i_sorted[,2:3],
                            (abs(res_i_sorted$logFC)>1) & (res_i_sorted$padj_fdr<0.01))
  c_names <- c(c_names,paste0(c("logFC_","FDR","Sign"),rep(pairw_comp[i],3)))
}

# Gene clustering
cutree(clustering.genes, k=n_clusters) %>% 
  as_tibble(rownames = "Gene") %>%
  dplyr::rename(cluster = value) -> clusters.n

saveRDS(clusters.n,paste0(args$outdir,"/clusters.",
                          n_clusters,".rds"))

clusters.n_sorted <- clusters.n[sort(clusters.n$Gene,index.return=T)$ix,]
clusters.n_sorted$cluster <- letters[order_gene_clusters][clusters.n_sorted$cluster]

DEG_table_sorted <- cbind(DEG_table_sorted,rep(NA,nrow(DEG_table_sorted)))
DEG_table_sorted[rownames(DEG_table_sorted)%in%clusters.n_sorted$Gene,ncol(DEG_table_sorted)] <- clusters.n_sorted$cluster 

c_names <- c(c_names,"gene_cluster")

# cluster specificity when OE to at least 3 other clusters
for (i in 1:6){
  cols_logFC <- which(grepl(i,c_names) & grepl("vs",c_names) & grepl("logFC",c_names))
  cols_FDR <-  which(grepl(i,c_names) & grepl("vs",c_names) & grepl("FDR",c_names))
  sign <- c(1,-1)[c(substr(gsub("[A-Z,a-z]+_","",c_names[cols_logFC]),1,1)==i)+1]

  DEG_table_sorted <- cbind(DEG_table_sorted,rowSums((t(t(DEG_table_sorted[,cols_logFC])*sign)>1) & (DEG_table_sorted[,cols_FDR]<0.01))>=3)
  c_names <- c(c_names,paste0("clust",i,"_spec"))
}

# Write table
DEG_table_sorted <- cbind(rownames(DEG_table_sorted),DEG_table_sorted)
colnames(DEG_table_sorted) <- c("Gene",c_names)
write.table(DEG_table_sorted,file=paste0(args$outdir,"/DEG_overview.txt"),
            col.names=T,row.names = F,sep="\t",quote=F)

####################
# Final clustering #
####################

clusters.n_ordered <- clusters.n
clusters.n_ordered$cluster <- as.numeric(clusters.n_ordered$cluster)
clusters.n_ordered$cluster_new <- order_gene_clusters[as.numeric(clusters.n_ordered$cluster)]
clusters.n_ordered_sorted <- clusters.n_ordered[sort(clusters.n_ordered$Gene,index.return=T)$ix,]
table(clusters.n_ordered$cluster,clusters.n_ordered$cluster_new)
clusters.n_ordered$cluster <- clusters.n_ordered$cluster_new
clusters.n_ordered <- clusters.n_ordered[,-3]

lapply(c(1:n_clusters), count.gaps, clusters.n_ordered) %>%
  do.call("rbind",.) %>% 
  as.vector -> clusters.gaps.n

cols3 <- rev(brewer.pal(10, "Spectral"))

annotation.colors <- list(cluster = c("1"=cols3[1], "2"=cols3[2], "3"=cols3[3], "4"=cols3[4], "5"=cols3[5],
                                      "6"=cols3[6], "7"=cols3[7], "8"=cols3[8], "9"=cols3[9], "10"=cols3[10])[1:n_clusters])

# pdf
pdf(file = paste0(args$outdir,"/heatmap_",
                  n_clusters,"geneclusters_final.pdf"),
    height = 6, width = 8)
annotation_row <- clusters.n_ordered %>% 
    column_to_rownames(var = "Gene") %>%
    arrange(cluster)
annotation_row$cluster <- as.factor(annotation_row$cluster)
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
           annotation_row = annotation_row,
           annotation_colors = annotation.colors,
           color = colorRampPalette(cols3)(100),breaks=seq(from=-2,to=2,length.out=101))
dev.off()

# pdf - nolabs
pdf(file = paste0(args$outdir,"/heatmap_",
                  n_clusters,"geneclusters_final_nolabs.pdf"),
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
           legend = FALSE,
           annotation_legend = FALSE,
           annotation_names_row	= FALSE,
           annotation_names_col = FALSE,
           color = colorRampPalette(cols3)(100),breaks=seq(from=-2,to=2,length.out=101))
dev.off()

print(clusters.gaps.n)

genes_ordered <- exp.values.scaled %>%
  rownames_to_column(var="Gene") %>% 
  full_join(clusters.n_ordered) %>%
  arrange(cluster) %>%
  column_to_rownames(var = "Gene") %>% rownames
write.table(genes_ordered,file=paste0(args$outdir,"/genes_gene_clustering_ordered.txt"),col.names = F,row.names = F,quote=F,sep="\n")

clusters.n_ordered_sorted$cluster_letter <- clusters.n_sorted$cluster
print(clusters.n_sorted)
table(clusters.n_ordered_sorted$cluster,clusters.n_ordered_sorted$cluster_letter)
table(clusters.n_ordered_sorted$cluster_new,clusters.n_ordered_sorted$cluster_letter)

#################################
# Violin plots per gene cluster #
#################################

print("Violin plots")
exp.values.scaled <- as.data.table(exp.values.scaled)
exp.values.scaled.melted <- melt(exp.values.scaled, variable.name = "variable", value.name = "value")

exp.values.scaled.melted <- cbind(exp.values.scaled.melted,
                                  rep(c("A","B","C","D","E","F","G", "I", "J")[clusters.n_ordered$cluster],6))
colnames(exp.values.scaled.melted)[3] <- "gene_cluster"
print(head(exp.values.scaled.melted))

g <- ggplot(exp.values.scaled.melted, aes(x=variable, y=value, fill=variable)) + 
  geom_violin(trim=FALSE) + 
  facet_wrap(~gene_cluster, scales="fixed") +
  scale_fill_manual(values=opts$celltype.colors) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  )

output_plot(g, sprintf("%s/violin_%sgeneclusters_final", args$outdir,n_clusters), height = 6, width = 8)

if(!dir.exists(sprintf("%s/violin_percluster",args$outdir))){dir.create(sprintf("%s/violin_percluster",args$outdir))}
for (cluster in unique(exp.values.scaled.melted$gene_cluster)){
  g <- ggplot(exp.values.scaled.melted[exp.values.scaled.melted$gene_cluster==cluster,], aes(x=variable, y=value, fill=variable)) + 
    geom_violin(trim=FALSE) + 
    scale_fill_manual(values=opts$celltype.colors) +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank()
    )
  output_plot(g, sprintf("%s/violin_percluster/violin_%sgeneclusters_cluster%s_final", args$outdir,n_clusters,cluster), height = 6, width = 8)
}


if (grepl("RNA",args$matrix)){
  correl_obj <- NULL
  for (t in c("GeneScoreMatrix_TSS","GeneScoreMatrix_distal")){
    print(t)
    args$sce_atac <- sprintf("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA_ATAC/pseudobulk/cluster/%s/pseudobulk_with_replicates.rds",
                             t)
    
    sce <- readRDS(args$sce_atac)
    exp.values_atac <- NULL
    cor_RNA_seq <- NULL
    for (i in 1:6){
      exp.values_atac <- cbind(exp.values_atac,rowMeans(assays(sce)[["logcounts"]][rownames(sce)%in%rownames(DE_genes)[rowSums(DE_genes)>0],sce$celltype==i]))
      obj_i <- merge(rownames_to_column(data.frame(val=exp.values[,i]),"gene"),rownames_to_column(data.frame(val_atac=exp.values_atac[,i]),"gene"),all=T)
      cor_RNA_seq <- cbind(cor_RNA_seq,c(cor.test(obj_i$val,obj_i$val_atac)$estimate,cor.test(obj_i$val,obj_i$val_atac)$p.value))
    }
    colnames(exp.values_atac) <- paste0(rep("atac_cluster",6),1:6)
    colnames(cor_RNA_seq) <- paste0(rep("cluster",6),1:6)
    rownames(cor_RNA_seq) <- c("rho","pval")
    print(summary(exp.values_atac))
    print(cor_RNA_seq)
    correl_obj <- rbind(correl_obj,cbind(colnames(cor_RNA_seq),t(cor_RNA_seq),rep(t,ncol(cor_RNA_seq))))

    t(apply(exp.values_atac, 1, z_score)) %>%
      as_tibble(rownames = "Gene") %>%
      column_to_rownames(var="Gene") -> exp.values_atac.scaled
    exp.values_atac.scaled <- exp.values_atac.scaled[!(rowSums(is.na(exp.values_atac.scaled))>0),]
    print(head(exp.values_atac.scaled))
    print(dim(exp.values_atac.scaled))
    print(summary(exp.values_atac.scaled))
    
    exp.values.scaled_multi <- exp.values.scaled %>% rownames_to_column(var="gene")
    dim(exp.values.scaled_multi)
    
    exp.values_atac.scaled <- exp.values_atac.scaled %>% rownames_to_column(var="gene")
    dim(exp.values_atac.scaled)
    
    exp.values.scaled_multi <- exp.values.scaled_multi %>% 
      merge(exp.values_atac.scaled,all.x=T) %>% dplyr::rename(Gene = gene)
    dim(exp.values.scaled_multi)

    print("Combined object")
    print(dim(exp.values.scaled_multi))
    print(head(exp.values.scaled_multi))

    count.gaps <- function(x,clusters){
      clusters %>%
        filter(cluster<=x) %>%
        nrow -> y
      return(y)
    }
    
    n_clusters <- 9
    
    cutree(clustering.genes, k=n_clusters) %>%
      as_tibble(rownames = "Gene") %>%
      dplyr::rename(cluster = value) -> clusters.n
    
    lapply(c(1:n_clusters), count.gaps, clusters.n_ordered %>% merge(exp.values.scaled_multi)) %>%
      do.call("rbind",.) %>%
      as.vector -> clusters.gaps.n
    
    annotation.colors <- list(cluster = c("1"=cols3[1], "2"=cols3[2], "3"=cols3[3], 
                                          "4"=cols3[4], "5"=cols3[5], "6"=cols3[6], 
                                          "7"=cols3[7], "8"=cols3[8], "9"=cols3[9], 
                                          "10"=cols3[10])[1:n_clusters])
    
    pdf(file = sprintf("%s/heatmap_%sgeneclusters_%s.pdf",
                       args$outdir,n_clusters,t), height = 6, width = 8)
    exp.values.scaled_multi %>%
      arrange(Gene) %>%
      full_join(clusters.n_ordered) %>%
      arrange(cluster) %>%
      select(-c("cluster","Gene")) %>%
      pheatmap(cluster_rows = F,
               cluster_cols = F,
               show_rownames = F,
               gaps_row = clusters.gaps.n,
               gaps_col = 5,
               annotation_colors = annotation.colors,
               color = colorRampPalette(cols3)(100),breaks=seq(from=-2,to=2,length.out=101))
    dev.off()
    if (FALSE) {
        exp.values.scaled_multi <- exp.values.scaled_multi %>% select(-Gene)
        exp.values.scaled_multi.melted <- melt(exp.values.scaled_multi)
        exp.values.scaled_multi.melted <- cbind(exp.values.scaled_multi.melted,
                                                rep(c("A","B","C","D","E","F","G", "I", "J")[clusters.n_ordered$cluster[sort(clusters.n_ordered$Gene,index.return=T)$ix]],6))
        colnames(exp.values.scaled_multi.melted)[3] <- "gene_cluster"
        print("Combined object (melted)")
        print(dim(exp.values.scaled_multi.melted))
        print(head(exp.values.scaled_multi.melted))
        print(str(exp.values.scaled_multi.melted))
    
        to.plot <- exp.values.scaled_multi.melted[grepl("atac", exp.values.scaled_multi.melted$variable), ]
        to.plot$variable <- as.factor(as.character(to.plot$variable))
    
        print(str(to.plot))
    
        g <- ggplot(to.plot, 
                    aes(x=variable, y=value, fill=variable)) + 
          geom_violin(trim=FALSE) + 
          facet_wrap(~gene_cluster, scales="fixed") +
          scale_fill_manual(values=opts$celltype.colors) +
          theme_classic()  +
          theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank()
          )
    
        output_plot(g, sprintf("%s/violin_%sgeneclusters_%s_final", args$outdir,n_clusters,t), height = 6, width = 8)
    
        for (cluster in unique(to.plot$gene_cluster)){
          g <- ggplot(to.plot[to.plot$gene_cluster==cluster,], aes(x=variable, y=value, fill=variable)) + 
            geom_violin(trim=FALSE) + 
            scale_fill_manual(values=opts$celltype.colors) +
            theme_classic() +
            theme(
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank()
            )
          output_plot(g, sprintf("%s/violin_percluster/violin_%sgeneclusters_%s_cluster%s_final", args$outdir,n_clusters,t,cluster), height = 6, width = 8)
        }
            
        g <- ggplot(exp.values.scaled_multi.melted, aes(x=variable, y=value, fill=variable)) + 
          geom_violin(trim=FALSE) + 
          facet_wrap(~gene_cluster, scales="fixed") +
          ggplot_theme_NoAxes()
        
        pdf(file = sprintf("%s/violin_%sgeneclusters_%s.pdf",args$outdir,n_clusters,t),
            height = 6, width = 8)
        print(g)
        dev.off()
    }
  }
  colnames(correl_obj) <- c("cluster","rho","pval","modality")
  correl_obj <- data.frame(correl_obj,stringsAsFactors = FALSE)
  correl_obj$rho <- as.numeric(correl_obj$rho)
  correl_obj$pval <- as.numeric(correl_obj$pval)
  correl_obj$modality <- as.factor(correl_obj$modality)
  print(correl_obj)
  g <- ggplot(data=correl_obj, aes(x=cluster, y=rho, fill=modality)) +
        geom_bar(stat="identity", position=position_dodge())
  pdf(file = sprintf("%s/correlation_barplot.pdf",args$outdir),
      height = 6, width = 8)
  print(g)
  dev.off()
} else if (grepl("GeneScore",args$matrix)){
  args$sce_RNA <- sprintf("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA_ATAC/pseudobulk/cluster/%s/pseudobulk_with_replicates.rds",
                          "RNA")
  
  sce <- readRDS(args$sce_RNA)
  exp.values_expr <- NULL
  for (i in 1:6){
    exp.values_expr <- cbind(exp.values_expr,rowMeans(assays(sce)[["logcounts"]][rownames(sce)%in%rownames(DE_genes)[rowSums(DE_genes)>0],sce$celltype==i]))
  }
  colnames(exp.values_expr) <- paste0(rep("expr_cluster",6),1:6)
  
  t(apply(exp.values_expr, 1, z_score)) %>%
    as_tibble(rownames = "Gene") %>%
    column_to_rownames(var="Gene") -> exp.values_expr.scaled
  exp.values_expr.scaled <- exp.values_expr.scaled[!(rowSums(is.na(exp.values_expr.scaled))>0),]
  head(exp.values_expr.scaled)
  dim(exp.values_expr.scaled)
  summary(exp.values_expr.scaled)
  
  exp.values.scaled_multi <- exp.values.scaled %>% rownames_to_column(var="gene")
  dim(exp.values.scaled_multi)
  
  exp.values_expr.scaled <- exp.values_expr.scaled %>% rownames_to_column(var="gene")
  dim(exp.values_expr.scaled)
  
  exp.values.scaled_multi <- exp.values.scaled_multi %>% 
    merge(exp.values_expr.scaled,all.x=T) %>% dplyr::rename(Gene = gene)
  dim(exp.values.scaled_multi)
  
  count.gaps <- function(x,clusters){
    clusters %>%
      filter(cluster<=x) %>%
      nrow -> y
    return(y)
  }
  
  n_clusters <- 9
  
  cutree(clustering.genes, k=n_clusters) %>%
    as_tibble(rownames = "Gene") %>%
    dplyr::rename(cluster = value) -> clusters.n
  
  lapply(c(1:n_clusters), count.gaps, clusters.n_ordered %>% merge(exp.values.scaled_multi)) %>%
    do.call("rbind",.) %>%
    as.vector -> clusters.gaps.n
  
  annotation.colors <- list(cluster = c("1"=cols3[1], "2"=cols3[2], "3"=cols3[3], 
                                        "4"=cols3[4], "5"=cols3[5], "6"=cols3[6], 
                                        "7"=cols3[7], "8"=cols3[8], "9"=cols3[9], 
                                        "10"=cols3[10])[1:n_clusters])
  
  pdf(file = sprintf("%s/heatmap_%sgeneclusters_expr.pdf",
                     args$outdir,n_clusters), height = 6, width = 8)
  exp.values.scaled_multi %>%
    arrange(Gene) %>%
    full_join(clusters.n_ordered) %>%
    arrange(cluster) %>%
    select(-c("cluster","Gene")) %>%
    pheatmap(cluster_rows = F,
             cluster_cols = F,
             show_rownames = F,
             gaps_row = clusters.gaps.n,
             gaps_col = 5,
             annotation_colors = annotation.colors,
             color = colorRampPalette(cols3)(100),breaks=seq(from=-2,to=2,length.out=101))
  dev.off()
  
  exp.values.scaled_multi <- exp.values.scaled_multi %>% select(-Gene)
  exp.values.scaled_multi.melted <- melt(exp.values.scaled_multi, variable.name = "variable", value.name = "value")
  exp.values.scaled_multi.melted <- cbind(exp.values.scaled_multi.melted,
                                          rep(c("A","B","C","D","E","F","G", "I", "J")[clusters.n_ordered$cluster[sort(clusters.n_ordered$Gene,index.return=T)$ix]],6))
  colnames(exp.values.scaled_multi.melted)[3] <- "gene_cluster"
  dim(exp.values.scaled_multi.melted)
  head(exp.values.scaled_multi.melted)
  
  g <- ggplot(exp.values.scaled_multi.melted, aes(x=variable, y=value, fill=variable)) + 
    geom_violin(trim=FALSE) + 
    facet_wrap(~gene_cluster, scales="fixed") +
    ggplot_theme_NoAxes()
  
  pdf(file = sprintf("%s/violin_%sgeneclusters_expr.pdf",
                     args$outdir,n_clusters),
      height = 6, width = 8)
  print(g)
  dev.off()
}

##########################################
## Peak distribution in genomic context ##
##########################################

 if (grepl("PeakMatrix",args$matrix)){

  print("Peak distribution")

  peak_metadata.dt <- fread(args$atac_peak_metadata) %>%
    .[,peak:=sprintf("%s:%s-%s",chr,start,end)] %>%
    .[,c("peak","peakType","distToGeneStart","nearestGene")]
  
  colnames(peak_metadata.dt)[1] <- "Gene"

  print(dim(DEG_table_sorted))
  DEG_table_sorted <- merge(DEG_table_sorted,peak_metadata.dt,all.x=TRUE)
  print(dim(DEG_table_sorted))
  print(colnames(DEG_table_sorted))

  to.plot <- as.data.frame(table(DEG_table_sorted$gene_cluster, DEG_table_sorted$peakType))
  colnames(to.plot) <- c("gene_cluster", "peakType", "V1")

  print(head(to.plot))

  for (j in unique(to.plot$gene_cluster)){
    to.plot$V1[to.plot$gene_cluster==j] <- to.plot$V1[to.plot$gene_cluster==j]/sum(to.plot$V1[to.plot$gene_cluster==j])
  }

  to.plot$V1 <- to.plot$V1*100

  p <- ggbarplot(to.plot, x="gene_cluster", fill="peakType", y="V1", position=position_dodge(width = 0.75)) +
    labs(x="", y="% of peaks present in genomic context") +
    theme(
      axis.title = element_text(size=rel(0.80)),
      legend.title = element_blank(),
      axis.text.x = element_text(size=rel(0.80), color="black"),
      axis.text.y = element_text(size=rel(0.80), color="black")
    )

  pdf(file.path(args$outdir,"barplots_number_DAP_vs_genomic_context.pdf"), width=7, height=5)
  print(p)
  dev.off()
 }

##################
# Print warnings #
##################

warnings()