source("/data/homes/louisc/Project_Babraham/RNA_ATAC/scripts/Settings.R")

#alt1 <- fread("/data/homes/louisc/Project_Babraham/RNA_ATAC/gene_regulatory_network/metacells/trajectories/N2P/sign_alt1/global_chip_GRN_sign_alt1_coef_score0.06.txt.gz")
alt2_m <- fread("/data/homes/louisc/Project_Babraham/RNA_ATAC/gene_regulatory_network/metacells/trajectories/N2P/TFandSign/alt2/merged/global_chip_GRN_sign_alt2_coef_score0.06.txt.gz")
alt2_um <- fread("/data/homes/louisc/Project_Babraham/RNA_ATAC/gene_regulatory_network/metacells/trajectories/N2P/TFandSign/alt2/unmerged/global_chip_GRN_sign_alt2_coef_score0.06.txt.gz")

alt4_m <- fread("/data/homes/louisc/Project_Babraham/RNA_ATAC/gene_regulatory_network/metacells/trajectories/N2P/TFandSign/alt4/merged/global_chip_GRN_sign_alt4_coef_score0.06.txt.gz")

alt5_m <- fread("/data/homes/louisc/Project_Babraham/RNA_ATAC/gene_regulatory_network/metacells/trajectories/N2P/TFandSign/alt5/merged/global_chip_GRN_sign_alt5_coef_score0.06.txt.gz")
alt5_um <- fread("/data/homes/louisc/Project_Babraham/RNA_ATAC/gene_regulatory_network/metacells/trajectories/N2P/TFandSign/alt5/unmerged/global_chip_GRN_sign_alt5_coef_score0.06.txt.gz")

genes_delay <- c("ZFP42", "ZNF281", "SOX4", "YBX1", "HMGA1", "ESRRB", "IRF3", "SOX11", "SOX2", "ZIC2")
genes_acceleration <- c("KLF5", "PBX4", "TFCP2L1", "PPARG", "NFKB1")

info_table <- matrix(rep(0,10),nrow=2)
colnames(info_table) <- 1:5
rownames(info_table) <- c("equal_coef","diff_coef")

alt1$sender_cluster[alt1$sender_cluster==5] <- 4
alt1$target_cluster[alt1$target_cluster==5] <- 4

alt21$sender_cluster[alt2$sender_cluster==5] <- 4
alt2$target_cluster[alt2$target_cluster==5] <- 4

mask_sender_eq_target <- (alt1$sender_cluster==alt1$target_cluster) & !is.na(alt1$sender_cluster_marker) & !is.na(alt1$target_cluster_marker)
info_table_all <- info_table
for (i in which(mask_sender_eq_target)){
  if (sum(paste(alt2$sender,alt2$target,sep="_")==paste0(alt1$sender[i],"_",alt1$target[i]))==0){
    next()
  }
  i_2 <- which(paste(alt2$sender,alt2$target,sep="_")==paste0(alt1$sender[i],"_",alt1$target[i]))
  if ((sign(alt1$beta[i])!=sign(alt2$beta[i_2])) & (abs(alt1$beta[i]-alt2$beta[i_2])>0.2) & ((alt1$pvalue[i]<0.1) & (alt2$pvalue[i_2]<0.1))){
    counter <- counter + 1 
    print("Opposite coef!")
    print(paste0("Alt1: From ",alt1$sender[i], " (",alt1$sender_cluster[i],") to ",
                 alt1$target[i]," (",alt1$target_cluster[i],
                 "), coef=",alt1$beta[i]," (pval=",alt1$pvalue[i],")."))
    print(paste0("Alt2: From ",alt2$sender[i_2], " (",alt2$sender_cluster[i_2],") to ",
                 alt2$target[i_2]," (",alt2$sender_cluster[i_2],
                 "), coef=",alt2$beta[i_2]," (pval=",alt2$pvalue[i_2],")."))
    info_table_all["diff_coef",alt1$sender_cluster[i]] <- info_table_all["diff_coef",alt1$sender_cluster[i]] + 1
  } else {
    info_table_all["equal_coef",alt1$sender_cluster[i]] <- info_table_all["equal_coef",alt1$sender_cluster[i]] + 1
  }
}   


mask_sender_eq_target_acceleration <- (alt1$sender_cluster==alt1$target_cluster)  & !is.na(alt1$sender_cluster_marker) & !is.na(alt1$target_cluster_marker) &
  ((alt1$sender%in%genes_acceleration) | (alt1$target%in%genes_acceleration))
info_table_acceleration <- info_table
for (i in which(mask_sender_eq_target_acceleration)){
  if (sum(paste(alt2$sender,alt2$target,sep="_")==paste0(alt1$sender[i],"_",alt1$target[i]))==0){
    next()
  }
  i_2 <- which(paste(alt2$sender,alt2$target,sep="_")==paste0(alt1$sender[i],"_",alt1$target[i]))
  if ((sign(alt1$beta[i])!=sign(alt2$beta[i_2])) & (abs(alt1$beta[i]-alt2$beta[i_2])>0.2) & ((alt1$pvalue[i]<0.1) & (alt2$pvalue[i_2]<0.1))){
    counter <- counter + 1 
    print("Opposite coef!")
    print(paste0("Alt1: From ",alt1$sender[i], " (",alt1$sender_cluster[i],") to ",
                 alt1$target[i]," (",alt1$target_cluster[i],
                 "), coef=",alt1$beta[i]," (pval=",alt1$pvalue[i],")."))
    print(paste0("Alt2: From ",alt2$sender[i_2], " (",alt2$sender_cluster[i_2],") to ",
                 alt2$target[i_2]," (",alt2$sender_cluster[i_2],
                 "), coef=",alt2$beta[i_2]," (pval=",alt2$pvalue[i_2],")."))
    info_table_acceleration["diff_coef",alt1$sender_cluster[i]] <- info_table_acceleration["diff_coef",alt1$sender_cluster[i]] + 1
  } else {
    info_table_acceleration["equal_coef",alt1$sender_cluster[i]] <- info_table_acceleration["equal_coef",alt1$sender_cluster[i]] + 1
  }
} 


mask_sender_eq_target_delay <- (alt1$sender_cluster==alt1$target_cluster)  & !is.na(alt1$sender_cluster_marker) & !is.na(alt1$target_cluster_marker) &
  ((alt1$sender%in%genes_delay) | (alt1$target%in%genes_delay))
info_table_delay <- info_table
for (i in which(mask_sender_eq_target_delay)){
  if (sum(paste(alt2$sender,alt2$target,sep="_")==paste0(alt1$sender[i],"_",alt1$target[i]))==0){
    next()
  }
  i_2 <- which(paste(alt2$sender,alt2$target,sep="_")==paste0(alt1$sender[i],"_",alt1$target[i]))
  if ((sign(alt1$beta[i])!=sign(alt2$beta[i_2])) & (abs(alt1$beta[i]-alt2$beta[i_2])>0.2) & ((alt1$pvalue[i]<0.1) & (alt2$pvalue[i_2]<0.1))){
    counter <- counter + 1 
    print("Opposite coef!")
    print(paste0("Alt1: From ",alt1$sender[i], " (",alt1$sender_cluster[i],") to ",
                 alt1$target[i]," (",alt1$target_cluster[i],
                 "), coef=",alt1$beta[i]," (pval=",alt1$pvalue[i],")."))
    print(paste0("Alt2: From ",alt2$sender[i_2], " (",alt2$sender_cluster[i_2],") to ",
                 alt2$target[i_2]," (",alt2$sender_cluster[i_2],
                 "), coef=",alt2$beta[i_2]," (pval=",alt2$pvalue[i_2],")."))
    info_table_delay["diff_coef",alt1$sender_cluster[i]] <- info_table_delay["diff_coef",alt1$sender_cluster[i]] + 1
  } else {
    info_table_delay["equal_coef",alt1$sender_cluster[i]] <- info_table_delay["equal_coef",alt1$sender_cluster[i]] + 1
  }
} 

info_table_all
info_table_delay
info_table_acceleration

alt1[(alt1$sender%in%genes_delay) & (!duplicated(alt1$sender)),c("sender","sender_cluster")]
alt1[(alt1$sender%in%genes_acceleration) & (!duplicated(alt1$sender)),c("sender","sender_cluster")]


########################
## Make hybrid object ##
######################

l_1 <- unique(paste(alt1$sender[alt1$sender_cluster==alt1$target_cluster],alt1$target[alt1$sender_cluster==alt1$target_cluster],sep="_"))
l_2 <- unique(paste(alt2$sender[alt2$sender_cluster==alt2$target_cluster],alt2$target[alt2$sender_cluster==alt2$target_cluster],sep="_"))

l_1[!(l_1%in%l_2)]
l_2[!(l_2%in%l_1)]

hybrid <- alt1 

########################
## Checks ##
######################

sum(paste(alt2$sender,alt2$target,sep="_")==paste(alt4$sender,alt4$target,sep="_"))
sum(alt2$beta==alt4$beta & !is.na(alt2$beta))

head(alt2[alt2$beta!=alt4$beta & !is.na(alt2$beta),])
head(alt4[alt2$beta!=alt4$beta & !is.na(alt2$beta),])

head(alt5_comp[alt5_comp$beta!=alt4$beta & !is.na(alt2$beta),])

# Alt 4_m vs 5_m
############

dim(alt4_m)
dim(alt5_m)

mask_alt4_both <- paste(alt4_m$sender,alt4_m$target,sep="_")%in%paste(alt5_m$sender,alt5_m$target,sep="_")
mask_alt5_both <- paste(alt5_m$sender,alt5_m$target,sep="_")%in%paste(alt4_m$sender,alt4_m$target,sep="_")

sum(paste(alt5_m$sender,alt5_m$target,sep="_")[mask_alt5_both]==paste(alt4_m$sender,alt4_m$target,sep="_")[mask_alt4_both])

alt4_m_comp <- alt4_m[mask_alt4_both,]
alt5_m_comp <- alt5_m[mask_alt5_both,]

head(alt4_m_comp[alt4_m_comp$beta!=alt5_m_comp$beta & !is.na(alt4_m_comp$beta),])
sum(alt4_m_comp$beta!=alt5_m_comp$beta & !is.na(alt4_m_comp$beta) & (alt4_m_comp$sender_cluster==alt4_m_comp$target_cluster))
sum(alt4_m_comp$beta!=alt5_m_comp$beta & !is.na(alt4_m_comp$beta) & (alt4_m_comp$sender_cluster<=alt4_m_comp$target_cluster))
sum(alt4_m_comp$beta!=alt5_m_comp$beta & !is.na(alt4_m_comp$beta) & (alt4_m_comp$target_cluster==1))

# Alt 5_m vs 5_um
############

dim(alt5_um)
dim(alt5_m)

sum(paste(alt5_m$sender,alt5_m$target,sep="_")==paste(alt5_um$sender,alt5_um$target,sep="_"))

head(alt5_m[alt5_m$beta!=alt5_um$beta & !is.na(alt5_m$beta),])
head(alt5_um[alt5_um$beta!=alt5_m$beta & !is.na(alt5_um$beta),])

sum(alt5_m$beta!=alt5_um$beta & !is.na(alt5_m$beta) & !(alt5_m$sender_cluster%in%c(4,5)) & !(alt5_m$target_cluster%in%c(4,5)))

sum(alt5_m$beta!=alt5_um$beta & !is.na(alt5_m$beta) & (alt5_m$sender_cluster==4))
sum(alt5_m$beta!=alt5_um$beta & !is.na(alt5_m$beta) & (alt5_m$sender_cluster==5))

sum(alt5_m$beta!=alt5_um$beta & !is.na(alt5_m$beta) & (alt5_m$target_cluster==4))
sum(alt5_m$beta!=alt5_um$beta & !is.na(alt5_m$beta) & (alt5_m$target_cluster==5))

# 3 to 4/5 connections alt2_um vs alt4_m
############

alt2_um_c3toc45 <- alt2_um[(alt2_um$sender_cluster_marker==3) & (alt2_um$target_cluster_marker%in%c(4:5)) & 
  !is.na(alt2_um$sender_cluster_marker) & !is.na(alt2_um$target_cluster_marker),]
dim(alt2_um_c3toc45)

alt4_m_c3toc45 <- alt4_m[(alt4_m$sender_cluster_marker==3) & (alt4_m$target_cluster_marker%in%c(4:5)) & 
  !is.na(alt4_m$sender_cluster_marker) & !is.na(alt4_m$target_cluster_marker),]
dim(alt4_m_c3toc45)

sum(paste(alt2_um_c3toc45$sender,alt2_um_c3toc45$target,sep="_")==paste(alt4_m_c3toc45$sender,alt4_m_c3toc45$target,sep="_"))
# data_diff <- data.frame(cbind(alt2_um_c3toc45$sender,alt2_um_c3toc45$target,alt2_um_c3toc45$sender_cluster,alt2_um_c3toc45$target_cluster,alt2_um_c3toc45$beta,
#                               alt4_m_c3toc45$beta)[alt2_um_c3toc45$beta!=alt4_m_c3toc45$beta,])
data_diff <- data.frame(cbind(alt2_um_c3toc45$sender,alt2_um_c3toc45$target,alt2_um_c3toc45$sender_cluster,alt2_um_c3toc45$target_cluster,
                              alt2_um_c3toc45$beta,alt2_um_c3toc45$pvalue,alt4_m_c3toc45$beta,alt4_m_c3toc45$pvalue)[abs(alt2_um_c3toc45$beta-alt4_m_c3toc45$beta)>0.2,])
colnames(data_diff) <- c("sender","target","cluster_sender","cluster_target","beta_alt2_um","pvalue_alt2_um","beta_alt4_m","pvalue_alt4_m")                              