################
##            ##
##  NicheNet  ##
##            ##
################

source("/data/homes/louisc/Project_Babraham/RNA_ATAC/scripts/Settings.R")

args <- list()
args$grn_coef <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/gene_regulatory_network/metacells/trajectories/N2P/global_chip_GRN_coef_score0.06.txt.gz"
args$min_coef <- 0.25
args$max_pval <- 0.10
args$motif_annotation <- "CISBP"
args$motifmatcher <- sprintf("/data/homes/louisc/Project_Babraham/ATAC/archR/Annotations/%s-Scores.rds",args$motif_annotation)
args$motif2gene <- sprintf("/data/homes/louisc/Project_Babraham/ATAC/archR/Annotations/%s_motif2gene.txt.gz",args$motif_annotation)
args$outdir <- "/data/homes/louisc/Project_Babraham"

##############
## Load TFs ##
##############

# all_TFs
all_TFs <- read.table("/data/homes/louisc/Project_Babraham/TFs.txt",sep="\n",header=F,quote="")$V1
print("All TFs")
print(head(all_TFs))

##################################
## Load global GRN coefficients ##
##################################

GRN_coef.dt <- fread(args$grn_coef) %>% 
  .[,gene:=toupper(gene)] %>% # .[gene%in%unique(marker_TFs_filt.dt$gene) & tf%in%unique(marker_TFs_filt.dt$gene)] %>%
  .[pvalue<args$max_pval & abs(beta)>=args$min_coef]

head(GRN_coef.dt)


tf2gene_chip.dt <- fread(args$tf2gene_virtual_chip) %>%
  .[chip_score>=args$min_chip_score & dist<=args$max_distance] 
dim(tf2gene_chip.dt)
tf2gene_chip.dt <- tf2gene_chip.dt %>% merge(GRN_coef.dt,all=TRUE) 
dim(tf2gene_chip.dt)
head(tf2gene_chip.dt)
tf2gene_chip.dt <- tf2gene_chip.dt[tf2gene_chip.dt$pvalue<args$max_pval & abs(tf2gene_chip.dt$beta)>=args$min_coef,]
dim(tf2gene_chip.dt)
head(tf2gene_chip.dt)


tf2gene_chip.dt$motif_score <- NA
tf2gene_chip.dt$motif_counts <- NA
motifScores.dt <- assay(motifmatcher.se,"motifScores")
motifCounts.dt <- assay(motifmatcher.se,"motifCounts")
colnames(motifScores.dt) <- gsub("_.+$","",colnames(motifScores.dt))
colnames(motifCounts.dt) <- gsub("_.+$","",colnames(motifCounts.dt))
sum(rownames(motifScores.dt)==rownames(motifCounts.dt))
sum(colnames(motifScores.dt)==colnames(motifCounts.dt))
print(nrow(tf2gene_chip.dt))
for (i in 1:nrow(tf2gene_chip.dt)){
  if (i%%100){print(i)}
  mask_i_row <- rownames(motifScores.dt)==tf2gene_chip.dt$peak[i]
  mask_i_col <- colnames(motifScores.dt)==tf2gene_chip.dt$tf[i]
  tf2gene_chip.dt$motif_score[i] <- as.numeric(motifScores.dt[mask_i_row,mask_i_col])
  tf2gene_chip.dt$motif_counts[i] <- as.numeric(motifCounts.dt[mask_i_row,mask_i_col])
}

write.table(tf2gene_chip.dt,paste0(args$outdir,"/tf2gene_chip.dt.txt"),col.names=T,row.names=T,quote=F,sep="\t")

tf2gene_chip.dt[sort(tf2gene_chip.dt$peak,index.return=T)$ix,]

n_repr <- 0
n_act <- 0
n_biv <- 0
for (i in unique(tf2gene_chip.dt$peak)){


n_repr <- 0
n_act <- 0
n_biv <- 0
for (i in unique(tf2gene_chip.dt$tf)){
  print(i)
  tf2gene_chip.dt_tf_i <- tf2gene_chip.dt[tf2gene_chip.dt$tf==i,]
  if (sum(duplicated(tf2gene_chip.dt_tf_i$peak))>0){
    print(i)
  }
  mask_cor_pos_tf_i <- tf2gene_chip.dt_tf_i$beta>0
  if(sum(mask_cor_pos_tf_i)==0){
    print("Solely repressive links")
    n_repr <- n_repr + 1
    print(sprintf("Av # peaks per gene: %.2f",mean(table(tf2gene_chip.dt_tf_i$gene[!mask_cor_pos_tf_i]))))
  } else if (sum(mask_cor_pos_tf_i)==nrow(tf2gene_chip.dt_tf_i)) {
    print("Solely activating links")
    n_act <- n_act + 1
    print(sprintf("Av # peaks per gene: %.2f",mean(table(tf2gene_chip.dt_tf_i$gene[mask_cor_pos_tf_i]))))
  } else {
    n_biv <- n_biv + 1
    if (length(table(tf2gene_chip.dt_tf_i$gene[mask_cor_pos_tf_i]))==1){
      print(sprintf("Av # peaks per gene: %.2f vs %.2f (p=%.2f)",
                    mean(table(tf2gene_chip.dt_tf_i$gene[mask_cor_pos_tf_i])),
                    mean(table(tf2gene_chip.dt_tf_i$gene[!mask_cor_pos_tf_i])),
                    t.test(table(tf2gene_chip.dt_tf_i$gene[!mask_cor_pos_tf_i]),mu=as.vector(table(tf2gene_chip.dt_tf_i$gene[mask_cor_pos_tf_i])))$p.value))
      if (var(tf2gene_chip.dt_tf_i$motif_score)!=0){
        print(sprintf("Av score: %.2f vs %.2f (p=%.2f)",
                      mean(tf2gene_chip.dt_tf_i$motif_score[mask_cor_pos_tf_i]),
                      mean(tf2gene_chip.dt_tf_i$motif_score[!mask_cor_pos_tf_i]),
                      t.test(tf2gene_chip.dt_tf_i$motif_score[!mask_cor_pos_tf_i],mu=tf2gene_chip.dt_tf_i$motif_score[mask_cor_pos_tf_i])$p.value))
      } else {
        print(sprintf("Av score: %.2f",
                      mean(tf2gene_chip.dt_tf_i$motif_score)))
      }
    } else if (length(table(tf2gene_chip.dt_tf_i$gene[!mask_cor_pos_tf_i]))==1){
      print(sprintf("Av # peaks per gene: %.2f vs %.2f (p=%.2f)",
                    mean(table(tf2gene_chip.dt_tf_i$gene[!mask_cor_pos_tf_i])),
                    mean(table(tf2gene_chip.dt_tf_i$gene[!mask_cor_pos_tf_i])),
                    t.test(table(tf2gene_chip.dt_tf_i$gene[mask_cor_pos_tf_i]),mu=as.vector(table(tf2gene_chip.dt_tf_i$gene[!mask_cor_pos_tf_i])))$p.value))
      if (var(tf2gene_chip.dt_tf_i$motif_score)!=0){
        print(sprintf("Av score: %.2f vs %.2f (p=%.2f)",
                      mean(tf2gene_chip.dt_tf_i$motif_score[!mask_cor_pos_tf_i]),
                      mean(tf2gene_chip.dt_tf_i$motif_score[!mask_cor_pos_tf_i]),
                      t.test(tf2gene_chip.dt_tf_i$motif_score[mask_cor_pos_tf_i],mu=tf2gene_chip.dt_tf_i$motif_score[!mask_cor_pos_tf_i])$p.value))
      } else {
        print(sprintf("Av score: %.2f",
                      mean(tf2gene_chip.dt_tf_i$motif_score)))
      }
    } else {
      print(sprintf("Av # peaks per gene: %.2f vs %.2f (p=%.2f)",
                    mean(table(tf2gene_chip.dt_tf_i$gene[mask_cor_pos_tf_i])),
                    mean(table(tf2gene_chip.dt_tf_i$gene[!mask_cor_pos_tf_i])),
                    t.test(table(tf2gene_chip.dt_tf_i$gene[mask_cor_pos_tf_i]),table(tf2gene_chip.dt_tf_i$gene[!mask_cor_pos_tf_i]))$p.value))
      if (var(tf2gene_chip.dt_tf_i$motif_score)!=0){
        print(sprintf("Av score: %.2f vs %.2f (p=%.2f)",
                      mean(tf2gene_chip.dt_tf_i$motif_score[mask_cor_pos_tf_i]),
                      mean(tf2gene_chip.dt_tf_i$motif_score[!mask_cor_pos_tf_i]),
                      t.test(tf2gene_chip.dt_tf_i$motif_score[mask_cor_pos_tf_i],tf2gene_chip.dt_tf_i$motif_score[!mask_cor_pos_tf_i])$p.value))
      } else {
        print(sprintf("Av score: %.2f",
                      mean(tf2gene_chip.dt_tf_i$motif_score)))
      }
    }
  }
}


tf2gene_chip.dt_ESRRB <- tf2gene_chip.dt[tf2gene_chip.dt$tf=="ESRRB",]
sum(duplicated(tf2gene_chip.dt_ESRRB$peak))
motifmatcher.se <- readRDS(args$motifmatcher)]