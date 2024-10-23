############################################
##                                        ##
##  Differential_expression_pseudobulk.R  ##
##                                        ##
############################################

source("/data/louisc/Project_Babraham/RNA/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',    type="character",    help='SingleCellExperiment file')
p$add_argument('--cluster_ids',            type="integer",     nargs='+',             help='Which cluster ids should be compared')
p$add_argument('--order_gene_clusters',            type="integer",     nargs='+',             help='Reorder gene clusters')
p$add_argument('--n_clusters',            type="integer",        help='How many gene clusters should be formed')
p$add_argument("--incl_samples",  type="character",    help='Which samples should be included')
p$add_argument("--rm_prev_res",  type="logical",default=T, help="Remove previous differential results?")
p$add_argument('--fdr_cutoff',        type="double",     default=0.01,     help='FDR cutoff to identify DE genes')
p$add_argument('--logfc_cutoff',        type="double",     default=1.00,     help='LogFC cutoff to identify DE genes')
p$add_argument('--min_expr',        type="integer",     default=4,     help='Minimal cpm value to retain expression values')
p$add_argument('--seed1',            type="integer",    default=42,                  help='Random seed')
p$add_argument('--outdir',   type="character",    help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

# START TEST
# args <- list()
# args$sce <- "/data/louisc/Project_Babraham/RNA/pseudobulk/cluster/SingleCellExperiment_pseudobulk_with_replicates.rds"
# args$sce_full <-  "/data/louisc/Project_Babraham/RNA/SingleCellExperiment.rds"
# args$cluster_ids <- 1:5
# args$rm_prev_res <- TRUE
# args$plot_expression <- FALSE
# args$fdr_cutoff <- 0.01
# args$logfc_cutoff <- 1
# args$min_expr <- 4 # 2**4 = 16, at least an average of 8 counts per milion for each group
# args$seed <- 42 
# END TEST

#####################
## Define settings ##
#####################

if(grepl("nodiff",args$incl_samples)){
  args$filter_differentiated <- TRUE
} else {
  args$filter_differentiated <- FALSE
}

############################
## Option mutliple groups ##
############################

args$groups <- combn(args$cluster_ids,2)
print(args$groups)

for (i in 1:ncol(args$groups)){
  args$groupA <- as.character(args$groups[1,i])
  args$groupB <- as.character(args$groups[2,i])
  args$outfile <-  sprintf("%s/res_cluster%svsclusters%s.txt", args$outdir, args$groupA, args$groupB)
  
  print(paste0("cluster ",args$groupA," vs cluster ",args$groupB))
  print(args$outfile)
  
  if(args$rm_prev_res & i==1){
    unlink(dirname(args$outfile), recursive = T, force = T)
  }
  
  dir.create(dirname(args$outfile), showWarnings = F)
  
  #####################
  ## Define settings ##
  #####################
  
  # Define groups
  opts <- list()
  opts$groups <- c(args$groupA,args$groupB)
  
  # stupid stuff but otherwise the snakemake pipeline doesn't work
  if (args$groupA==args$groupB) {
    out <- data.table(feature=NA, logFC=NA, padj_fdr=NA)
    fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
    warning("groupA and groupB are the same, saving an empty file...")
    quit(status=0)
  }
  
  #########################
  ## Load RNA expression ##
  #########################
  
  # Load SingleCellExperiment object
  sce <- readRDS(args$sce)
  
  if (i==1){
    print(sce)
  }
  
  # temporary
  if (!"celltype" %in% colnames(colData(sce))) {
    sce$celltype <- colnames(sce) %>% strsplit("_rep") %>% map_chr(1)
    sce$batch <- c("Jasmin","SEQ1","SEQ2")
  }
  
  sce <- sce[,sce$celltype %in% opts$groups]
  
  sce$celltype <- factor(sce$celltype, levels=opts$groups)
  print(table(sce$celltype))
  
  if (i == 1){
    DE_genes <- matrix(rep(F,nrow(sce)*ncol(args$groups)),ncol=ncol(args$groups))
    rownames(DE_genes) <- sort(rownames(sce))
    colnames(DE_genes) <- paste0(args$groups[1,],rep("vs",ncol(args$groups)),args$groups[2,])
  }

  #########################################
  ## Calculate average expression levels ##
  #########################################
  
  expr.dt <- data.table(
    gene = rownames(sce),
    mean_groupA = rowMeans(logcounts(sce[,sce$celltype==args$groupA])) %>% round(2),
    mean_groupB = rowMeans(logcounts(sce[,sce$celltype==args$groupB])) %>% round(2)
  )
  
  #######################
  ## Feature selection ##
  #######################
  
  genes.to.use <- expr.dt[mean_groupA>=args$min_expr | mean_groupB>=args$min_expr,gene]
  
  ################################################
  ## Differential expression testing with edgeR ##
  ################################################
  
  # Convert SCE to DGEList
  sce_edger <- scran::convertTo(sce[genes.to.use,], type="edgeR")
  
  # Define design matrix (with intercept)
  design <- model.matrix(~sce$celltype)
  
  # Estimate dispersions
  sce_edger  <- estimateDisp(sce_edger,design)
  
  # Fit GLM
  fit <- glmQLFit(sce_edger,design)
  
  # Likelihood ratio test
  lrt <- glmQLFTest(fit)
  
  # Construct output data.frame
  out <- topTags(lrt, n=nrow(lrt))$table %>% as.data.table(keep.rownames=T) %>%
    setnames(c("gene","logFC","logCPM","LR","p.value","padj_fdr")) %>%
    .[,c("logCPM","LR","p.value"):=NULL] %>%
    .[,c("padj_fdr","logFC"):=list(signif(padj_fdr,digits=3), round(logFC,3))] %>%
    merge(expr.dt, by="gene", all.y=TRUE) %>%
    setorder(padj_fdr, na.last=T)
  
  print(head(out))
  
  mask_out <- (out$padj_fdr<args$fdr_cutoff) & (abs(out$logFC)>args$logfc_cutoff)
  DE_genes[rownames(DE_genes)%in%out$gene[mask_out],i] <- T
  
  print(sum(mask_out))
  
  ##################
  ## Save results ##
  ##################
  
  fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
}

print(dim(DE_genes))
print(sum(rowSums(DE_genes)>0))
print(head(DE_genes))
saveRDS(DE_genes,paste0(args$outdir,"/DE_genes.rds"))

###################
# Gene clustering #
###################

sce <- readRDS(args$sce)
exp.values <- NULL
for (i in 1:5){
  exp.values <- cbind(exp.values,rowMeans(assays(sce)[["logcounts"]][rownames(sce)%in%rownames(DE_genes)[rowSums(DE_genes)>0],sce$celltype==i]))
}
colnames(exp.values) <- paste0(rep("cluster",5),1:5)

max_exp_row <- apply(assays(sce)[["logcounts"]],1,max)
pdf(file = paste0(args$outdir,"/hist_expr.pdf"), height = 6, width = 8)
hist(max_exp_row,breaks=75,xlim=c(0.1,12))
# hist(max_exp_row,breaks=75)
dev.off()

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

pdf(file =  paste0(args$outdir,"/elbow.plot.pdf"), height = 6, width = 8)
plot(1:50, wss, 
     type="b", 
     xlab="Number of Clusters",
     ylab="Within groups sum of squares")
dev.off()


# Then we perform clustering and select various numbers of clusters
set.seed(args$seed1)
Rclusterpp.hclust(exp.values.scaled) -> clustering.genes

saveRDS(clustering.genes, paste0(args$outdir,"/clustering.genes.rds"))

### This is to visulalise the tree
# plot(clustering.genes)

### Then we look at various options, how many clusters to extract

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
                                        "6"=cols3[6], "7"=cols3[7], "8"=cols3[8], "9"=cols3[9], "10"=cols3[10]
  ))
  annotation.colors <- annotation.colors[1:i]
  
  pdf(file = sprintf("%s/heatmap_%sgeneclusters.pdf", 
                     args$outdir, i),
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
pairw_comp <- paste0(args$groups[1,],rep("vs",ncol(args$groups)),args$groups[2,])

# Average expression
sce <- readRDS(args$sce)
for (i in args$cluster_ids){
  DEG_table <- cbind(DEG_table,rowMeans(assays(sce)[["logcounts"]][,sce$celltype==i]))
  c_names <- c(c_names,paste0("log_av_",i))
}

DEG_table_sorted <- DEG_table[sort(rownames(DEG_table),index.return=T)$ix,]

# DE results
for (i in 1:length(pairw_comp)){
  res_i <- read.table(sprintf("%s/res_cluster%svsclusters%s.txt",
                              args$outdir, substr(pairw_comp[i],1,1), substr(pairw_comp[i],4,4)),
                      header=T, sep="\t")
  res_i_sorted <- res_i[sort(res_i$gene,index.return=T)$ix,]
  DEG_table_sorted <- cbind(DEG_table_sorted,res_i_sorted[,2:3],
                            (abs(res_i_sorted$logFC)>1) & (res_i_sorted$padj_fdr<0.01))
  c_names <- c(c_names,paste0(c("logFC_","FDR","Sign"),rep(pairw_comp[i],3)))
}

# Gene clustering
cutree(clustering.genes, k=args$n_clusters) %>% 
  as_tibble(rownames = "Gene") %>%
  dplyr::rename(cluster = value) -> clusters.n

saveRDS(clusters.n,paste0(args$outdir,"/clusters.n.rds"))

clusters.n_sorted <- clusters.n[sort(clusters.n$Gene,index.return=T)$ix,]
clusters.n_sorted$cluster <- c("A","B","C","D","E","F","G")[args$order_gene_clusters][clusters.n_sorted$cluster]

DEG_table_sorted <- cbind(DEG_table_sorted,rep(NA,nrow(DEG_table_sorted)))
DEG_table_sorted[rownames(DEG_table_sorted)%in%clusters.n_sorted$Gene,ncol(DEG_table_sorted)] <- clusters.n_sorted$cluster 

c_names <- c(c_names,"gene_cluster")

# Cluster spec
for (i in args$cluster_ids){
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
clusters.n_ordered$cluster_new <- args$order_gene_clusters[as.numeric(clusters.n_ordered$cluster)]
clusters.n_ordered_sorted <- clusters.n_ordered[sort(clusters.n_ordered$Gene,index.return=T)$ix,]
table(clusters.n_ordered$cluster,clusters.n_ordered$cluster_new)
clusters.n_ordered$cluster <- clusters.n_ordered$cluster_new
clusters.n_ordered <- clusters.n_ordered[,-3]

lapply(c(1:args$n_clusters), count.gaps, clusters.n_ordered) %>%
  do.call("rbind",.) %>% 
  as.vector -> clusters.gaps.n

cols3 <- rev(brewer.pal(10, "Spectral"))

annotation.colors <- list(cluster = c("1"=cols3[1], "2"=cols3[2], "3"=cols3[3], "4"=cols3[4], "5"=cols3[5],
                                      "6"=cols3[6], "7"=cols3[7], "8"=cols3[8], "9"=cols3[9], "10"=cols3[10]
))
annotation.colors <- annotation.colors[1:7]

pdf(file = paste0(args$outdir,"/heatmap_",args$n_clusters,"geneclusters_final.pdf"),
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

clusters.n_ordered_sorted$cluster_letter <- clusters.n_sorted$cluster
table(clusters.n_ordered_sorted$cluster,clusters.n_ordered_sorted$cluster_letter)
table(clusters.n_ordered_sorted$cluster_new,clusters.n_ordered_sorted$cluster_letter)

#################################
# Violin plots per gene cluster #
#################################

exp.values.scaled.melted <- exp.values.scaled %>% melt
exp.values.scaled.melted <- cbind(exp.values.scaled.melted,
                                  rep(c("A","B","C","D","E","F","G")[clusters.n_ordered$cluster],5))
colnames(exp.values.scaled.melted)[3] <- "gene_cluster"

pdf(file = paste0(args$outdir,"/violin_",args$n_clusters,"geneclusters_final.pdf"),
    height = 6, width = 8)
ggplot(exp.values.scaled.melted, aes(x=variable, y=value, fill=variable)) + 
  geom_violin(trim=FALSE) + 
  facet_wrap(~gene_cluster, scales="fixed") +
  ggplot_theme_NoAxes()
dev.off()


