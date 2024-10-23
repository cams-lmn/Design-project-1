##########
## conv ##
##########

myfun_conv <- function(i){ #conv
  if(which(genes==i)%%round(length(genes)/10)==0){print(paste0(round(which(genes==i)/length(genes)*100),"% done!"))}
  tfs <- tmp[tmp$gene==i,tf]
  dat_return <- NULL

  for (j in tfs){
    x <- rna_tf.mtx[j,]
    y <- rna_targets.mtx[i,]
    
    if (args$plot_correlations & ((i%in%markers_TF_undupl$TF) & (j%in%markers_TF_undupl$TF))){
      to.plot <- data.frame(tf=rep(j,length(rna_tf.mtx[j,])),gene=rep(i,length(rna_tf.mtx[j,])),
                            x=rna_tf.mtx[j,],y=rna_targets.mtx[i,],
                            sample=factor(gsub("#.+$","",names(rna_targets.mtx[i,])),
                                          levels=c("d0","d1","d3","d5","d7","d10","d14","d18")))
      
      p1 <- ggscatter(to.plot, x="x", y="y", fill="sample", size=5, shape=21, 
                      add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
        stat_cor(method = "pearson", label.x.npc = 0.5, label.y.npc = 0.05) +
        stat_regline_equation(label.x.npc = 0.5, label.y.npc = 0.1) +
        scale_fill_manual(values=opts$color_scheme[1:length(unique(to.plot$sample))]) +
        labs(x=sprintf("%s expression",j), y=sprintf("%s expression",i)) +
        theme_classic() +
        theme(
          axis.text = element_text(size=rel(0.7))
        )
      
      pdf(sprintf("%s/plots_cor_score%s/%s_vs_%s.pdf",args$outdir,args$min_chip_score,j,i), 
          width = 7, height = 5)
      suppressMessages(print(p1))
      dev.off()
    }
    
    lm.fit <- lm(y~x)
    if (is.na(coef(lm.fit)[2])){
      dat_return <- rbind(dat_return,c(j, i, NA, NA, paste(1:5,collapse="-")))
    } else {
      dat_return <- rbind(dat_return,c(j, i, 
                                      round(coef(lm.fit)[[2]],3), 
                                      format(summary(lm.fit)$coefficients[2,4], digits=3),
                                      paste(1:5,collapse="-")))
    }
  }
  if (is.null(dim(dat_return))){
    return(NULL)
  } else {
    colnames(dat_return) <- c("tf","gene","beta","pvalue","clusters_for_calc")
    dat_return <- data.frame(dat_return,stringsAsFactors = F)
    return(dat_return)
  }
}

##########
## alt1 ##
##########

myfun_alt1 <- function(i){ #alt1
  if(which(genes==i)%%round(length(genes)/10)==0){print(paste0(round(which(genes==i)/length(genes)*100),"% done!"))}
  tfs <- tmp[tmp$gene==i,tf]
  dat_return <- NULL

  clust_i <- gene_cluster_annot$cluster[gene_cluster_annot$gene==i]
  
  # print("Number of NULLs")
  # print(sum(is.null(tfs)))
  for (j in tfs){
    clust_j <- gene_cluster_annot$cluster[gene_cluster_annot$gene==j]

    cluster_range <- sort(c(clust_i,clust_j))
    cluster_range <- cluster_range[1]:cluster_range[2]

    if (sum(gene_cluster_annot$cluster==5)==0){
      cluster_sce_mc <- sce_mc$cluster
      cluster_sce_mc[cluster_sce_mc==5] <- 4
      mask_cluster <- cluster_sce_mc%in%cluster_range
    } else {
      mask_cluster <- sce_mc$cluster%in%cluster_range
    }
    
    x <- rna_tf.mtx[j,mask_cluster]
    y <- rna_targets.mtx[i,mask_cluster]
    
    if (args$plot_correlations & ((i%in%c(markers_TF_undupl$TF,sign_senders)) & (j%in%c(markers_TF_undupl$TF,sign_senders)))){
      # print(paste0(j," vs ",i))
      to.plot <- data.frame(tf=rep(j,length(rna_tf.mtx[j,])),gene=rep(i,length(rna_tf.mtx[j])),
                            x=rna_tf.mtx[j,],y=rna_targets.mtx[i,],
                            sample=factor(gsub("#.+$","",names(rna_targets.mtx[i,])),
                                          levels=c("d0","d1","d3","d5","d7","d10","d14","d18")))
      # print(head(to.plot))
      
      p1 <- ggscatter(to.plot, x="x", y="y", fill="sample", size=5, shape=21)  +
              scale_fill_manual(values=opts$color_scheme[1:length(unique(to.plot$sample))]) +
              labs(x=sprintf("%s expression",j), y=sprintf("%s expression",i)) +
              theme_classic() +
              theme(axis.text = element_text(size=rel(0.7))) +
              geom_smooth(data = to.plot[mask_cluster,], mapping=aes(x=x,y=y), method = "lm",  inherit.aes = FALSE) +
              stat_regline_equation(data = to.plot[mask_cluster,], mapping=aes(x=x,y=y), label.x.npc = 0.5, label.y.npc = 0.1,  inherit.aes = FALSE) +
              stat_cor(data = to.plot[mask_cluster,], mapping=aes(x=x, y=y), method = "pearson", label.x.npc = 0.5, label.y.npc = 0.05) 
        
      pdf(sprintf("%s/plots_cor_score%s/%s_vs_%s.pdf",args$outdir,args$min_chip_score,j,i), 
          width = 7, height = 5)
      suppressMessages(print(p1))
      dev.off()
    }
    
    lm.fit <- lm(y~x)
    if (is.na(coef(lm.fit)[2])){
      dat_return <- rbind(dat_return,c(j, i, NA, NA, paste(cluster_range,collapse="-")))
    } else {
      dat_return <- rbind(dat_return,c(j, i, 
                                      round(coef(lm.fit)[[2]],3), 
                                      format(summary(lm.fit)$coefficients[2,4], digits=3),
                                      paste(cluster_range,collapse="-")))
    }
  }
  if (is.null(dim(dat_return))){
    return(NULL)
  } else {
    colnames(dat_return) <- c("tf","gene","beta","pvalue","clusters_for_calc")
    dat_return <- data.frame(dat_return,stringsAsFactors = F)
    return(dat_return)
  }
}

##########
## alt2 ##
##########

myfun_alt2 <- function(i){ #alt2
  if(which(genes==i)%%round(length(genes)/10)==0){print(paste0(round(which(genes==i)/length(genes)*100),"% done!"))}
  tfs <- tmp[tmp$gene==i,tf]
  dat_return <- NULL

  clust_i <- gene_cluster_annot$cluster[gene_cluster_annot$gene==i]
  
  # print("Number of NULLs")
  # print(sum(is.null(tfs)))
  for (j in tfs){
    clust_j <- gene_cluster_annot$cluster[gene_cluster_annot$gene==j]

    cluster_range <- max(c(clust_i,clust_j))
    cluster_range <- 1:cluster_range

    if (sum(gene_cluster_annot$cluster==5)==0){
      cluster_sce_mc <- sce_mc$cluster
      cluster_sce_mc[cluster_sce_mc==5] <- 4
      mask_cluster <- cluster_sce_mc%in%cluster_range
    } else {
      mask_cluster <- sce_mc$cluster%in%cluster_range
    }
    
    x <- rna_tf.mtx[j,mask_cluster]
    y <- rna_targets.mtx[i,mask_cluster]
    
    if (args$plot_correlations & ((i%in%c(markers_TF_undupl$TF,sign_senders)) & (j%in%c(markers_TF_undupl$TF,sign_senders)))){
      # print(paste0(j," vs ",i))
      to.plot <- data.frame(tf=rep(j,length(rna_tf.mtx[j,])),gene=rep(i,length(rna_tf.mtx[j,])),
                            x=rna_tf.mtx[j,],y=rna_targets.mtx[i,],
                            sample=factor(gsub("#.+$","",names(rna_targets.mtx[i,])),
                                          levels=c("d0","d1","d3","d5","d7","d10","d14","d18")))
      #print(head(to.plot))
      
      p1 <- ggscatter(to.plot, x="x", y="y", fill="sample", size=5, shape=21)  +
              scale_fill_manual(values=opts$color_scheme[1:length(unique(to.plot$sample))]) +
              labs(x=sprintf("%s expression",j), y=sprintf("%s expression",i)) +
              theme_classic() +
              theme(axis.text = element_text(size=rel(0.7))) +
              geom_smooth(data = to.plot[mask_cluster,], mapping=aes(x=x,y=y), method = "lm",  inherit.aes = FALSE) +
              stat_regline_equation(data = to.plot[mask_cluster,], mapping=aes(x=x,y=y), label.x.npc = 0.5, label.y.npc = 0.1,  inherit.aes = FALSE) +
              stat_cor(data = to.plot[mask_cluster,], mapping=aes(x=x, y=y), method = "pearson", label.x.npc = 0.5, label.y.npc = 0.05) 
        
      pdf(sprintf("%s/plots_cor_score%s/%s_vs_%s.pdf",args$outdir,args$min_chip_score,j,i), 
          width = 7, height = 5)
      suppressMessages(print(p1))
      dev.off()
    }
    
    lm.fit <- lm(y~x)
    if (is.na(coef(lm.fit)[2])){
      dat_return <- rbind(dat_return,c(j, i, NA, NA, paste(cluster_range,collapse="-")))
    } else {
      dat_return <- rbind(dat_return,c(j, i, 
                                      round(coef(lm.fit)[[2]],3), 
                                      format(summary(lm.fit)$coefficients[2,4], digits=3),
                                      paste(cluster_range,collapse="-")))
    }
  }
  if (is.null(dim(dat_return))){
    return(NULL)
  } else {
    colnames(dat_return) <- c("tf","gene","beta","pvalue","clusters_for_calc")
    # print(dat_return[1,])
    dat_return <- data.frame(dat_return,stringsAsFactors = F)
    return(dat_return)
  }
}

##########
## alt3 ##
##########

myfun_alt3 <- function(i){ #alt3
  if(which(genes==i)%%round(length(genes)/10)==0){print(paste0(round(which(genes==i)/length(genes)*100),"% done!"))}
  tfs <- tmp[tmp$gene==i,tf]
  dat_return <- NULL

  clust_i <- gene_cluster_annot$cluster[gene_cluster_annot$gene==i]
  
  # print("Number of NULLs")
  # print(sum(is.null(tfs)))
  for (j in tfs){
    clust_j <- gene_cluster_annot$cluster[gene_cluster_annot$gene==j]

    cluster_range <- sort(c(clust_i,clust_j))
    cluster_range <- cluster_range[1]:cluster_range[2]
    
    if (sum(gene_cluster_annot$cluster==5)==0){
      cluster_sce_mc <- sce_mc$cluster
      cluster_sce_mc[cluster_sce_mc==5] <- 4
      mask_cluster <- cluster_sce_mc%in%cluster_range
    } else {
      mask_cluster <- sce_mc$cluster%in%cluster_range
    }
    
    x <- rna_tf.mtx[j,mask_cluster]
    y <- rna_targets.mtx[i,mask_cluster]

    clusters_x <- sce_mc$cluster[mask_cluster]

    if (length(unique(clusters_x))>1){
      lm.fit <- lm(y~x+as.factor(clusters_x))
    } else {
      lm.fit <- lm(y~x)
    }
    
    if (args$plot_correlations & ((i%in%c(markers_TF_undupl$TF,sign_senders)) & (j%in%c(markers_TF_undupl$TF,sign_senders)))){
      # print(paste0(j," vs ",i))
      to.plot <- data.frame(tf=rep(j,length(rna_tf.mtx[j,mask_cluster])),gene=rep(i,length(rna_tf.mtx[j,mask_cluster])),
                            x=rna_tf.mtx[j,mask_cluster],y=rna_targets.mtx[i,mask_cluster],
                            sample=factor(gsub("#.+$","",names(rna_targets.mtx[i,mask_cluster])),
                                          levels=c("d0","d1","d3","d5","d7","d10","d14","d18")),
                            cluster=as.factor(clusters_x))
      # print(head(to.plot))
      
      p1 <- ggscatter(to.plot, x="x", y="y", fill="cluster", size=5, shape=21, 
                      add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
        stat_cor(method = "pearson", label.x.npc = 0.5, label.y.npc = 0.05) +
        stat_regline_equation(label.x.npc = 0.5, label.y.npc = 0.1) +
        scale_fill_manual(values=opts$color_scheme[1:length(unique(to.plot$cluster))]) +
        labs(x=sprintf("%s expression",j), y=sprintf("%s expression",i)) +
        theme_classic() +
        theme(
          axis.text = element_text(size=rel(0.7))
        ) + geom_abline(intercept=lm.fit$coef[1],slope=lm.fit$coef[2],color=opts$color_scheme[1]) +
        geom_abline(intercept=lm.fit$coef[1]+lm.fit$coef[3],slope=lm.fit$coef[2],color=opts$color_scheme[2]) +
        geom_abline(intercept=lm.fit$coef[1]+lm.fit$coef[4],slope=lm.fit$coef[2],color=opts$color_scheme[3]) +
        geom_abline(intercept=lm.fit$coef[1]+lm.fit$coef[5],slope=lm.fit$coef[2],color=opts$color_scheme[4]) +
        geom_abline(intercept=lm.fit$coef[1]+lm.fit$coef[6],slope=lm.fit$coef[2],color=opts$color_scheme[5])
      
      pdf(sprintf("%s/plots_cor_score%s/%s_vs_%s.pdf",args$outdir,args$min_chip_score,j,i), 
          width = 7, height = 5)
      suppressMessages(print(p1))
      dev.off()
    }

    if (is.na(coef(lm.fit)[2])){
      dat_return <- rbind(dat_return,c(j, i, NA, NA, paste(cluster_range,collapse="-")))
    } else {
      dat_return <- rbind(dat_return,c(j, i, 
                                      round(coef(lm.fit)[[2]],3), 
                                      format(summary(lm.fit)$coefficients[2,4], digits=3),
                                      paste(cluster_range,collapse="-")))
    }
  }
  if (is.null(dim(dat_return))){
    return(NULL)
  } else {
    colnames(dat_return) <- c("tf","gene","beta","pvalue","clusters_for_calc")
    # print(dat_return[1,])
    dat_return <- data.frame(dat_return,stringsAsFactors = F)
    return(dat_return)
  }
}


##########
## alt4 ##
##########

myfun_alt4 <- function(i){ #alt4
  if(which(genes==i)%%round(length(genes)/10)==0){print(paste0(round(which(genes==i)/length(genes)*100),"% done!"))}
  tfs <- tmp[tmp$gene==i,tf]
  dat_return <- NULL

  clust_i <- gene_cluster_annot$cluster[gene_cluster_annot$gene==i]
  
  # print("Number of NULLs")
  # print(sum(is.null(tfs)))
  for (j in tfs){
    clust_j <- gene_cluster_annot$cluster[gene_cluster_annot$gene==j]

    if (clust_i==clust_j){
      cluster_range <- sort(c(clust_i,clust_j))
      cluster_range <- cluster_range[1]:cluster_range[2]
    } else {
      cluster_range <- max(c(clust_i,clust_j))
      cluster_range <- 1:cluster_range
    }

    if (sum(gene_cluster_annot$cluster==5)==0){
      cluster_sce_mc <- sce_mc$cluster
      cluster_sce_mc[cluster_sce_mc==5] <- 4
      mask_cluster <- cluster_sce_mc%in%cluster_range
    } else {
      mask_cluster <- sce_mc$cluster%in%cluster_range
    }

    if(length(mask_cluster)!=ncol(rna_tf.mtx)){print("Alarm!")}
    
    x <- rna_tf.mtx[j,mask_cluster]
    y <- rna_targets.mtx[i,mask_cluster]
    
    if (args$plot_correlations & ((i%in%c(markers_TF_undupl$TF,sign_senders)) & (j%in%c(markers_TF_undupl$TF,sign_senders)))){
        # print(paste0(j," vs ",i))
      to.plot <- data.frame(tf=rep(j,length(rna_tf.mtx[j,])),gene=rep(i,length(rna_tf.mtx[j])),
                            x=rna_tf.mtx[j,],y=rna_targets.mtx[i,],
                            sample=factor(sce_mc$cluster))
      # print(head(to.plot))
      
      p1 <- ggscatter(to.plot, x="x", y="y", fill="sample", size=5, shape=21)  +
              scale_fill_manual(values=opts$celltype.colors[1:length(unique(to.plot$sample))]) +
              labs(x=sprintf("%s expression",j), y=sprintf("%s expression",i)) +
              theme_classic() +
              theme(axis.text = element_text(size=rel(0.7))) +
              geom_smooth(data = to.plot[mask_cluster,], mapping=aes(x=x,y=y), method = "lm",  inherit.aes = FALSE) +
              stat_regline_equation(data = to.plot[mask_cluster,], mapping=aes(x=x,y=y), label.x.npc = 0.5, label.y.npc = 0.1,  inherit.aes = FALSE) +
              stat_cor(data = to.plot[mask_cluster,], mapping=aes(x=x, y=y), method = "pearson", label.x.npc = 0.5, label.y.npc = 0.05) 
        
      pdf(sprintf("%s/plots_cor_score%s/%s_vs_%s.pdf",args$outdir,args$min_chip_score,j,i), 
          width = 7, height = 5)
      suppressMessages(print(p1))
      dev.off()
    }
    
    lm.fit <- lm(y~x)
    if (is.na(coef(lm.fit)[2])){
      dat_return <- rbind(dat_return,c(j, i, NA, NA, paste(cluster_range,collapse="-")))
    } else {
      dat_return <- rbind(dat_return,c(j, i, 
                                      round(coef(lm.fit)[[2]],3), 
                                      format(summary(lm.fit)$coefficients[2,4], digits=3),
                                      paste(cluster_range,collapse="-")))
    }
  }
  if (is.null(dim(dat_return))){
    return(NULL)
  } else {
    colnames(dat_return) <- c("tf","gene","beta","pvalue","clusters_for_calc")
    dat_return <- data.frame(dat_return,stringsAsFactors = F)
    return(dat_return)
  }
}


##########
## alt5 ##
##########

myfun_alt5 <- function(i){ #alt5
  if(which(genes==i)%%round(length(genes)/10)==0){print(paste0(round(which(genes==i)/length(genes)*100),"% done!"))}
  tfs <- tmp[tmp$gene==i,tf]
  dat_return <- NULL

  clust_i <- gene_cluster_annot$cluster[gene_cluster_annot$gene==i]
  
  # print("Number of NULLs")
  # print(sum(is.null(tfs)))
  for (j in tfs){
    clust_j <- gene_cluster_annot$cluster[gene_cluster_annot$gene==j]

    if (clust_j>=clust_i){
      cluster_range <- sort(c(clust_i,clust_j))
      cluster_range <- cluster_range[1]:cluster_range[2]
    } else {
      cluster_range <- max(c(clust_i,clust_j))
      cluster_range <- 1:cluster_range
    }

    if (sum(gene_cluster_annot$cluster==5)==0){
      cluster_sce_mc <- sce_mc$cluster
      cluster_sce_mc[cluster_sce_mc==5] <- 4
      mask_cluster <- cluster_sce_mc%in%cluster_range
    } else {
      mask_cluster <- sce_mc$cluster%in%cluster_range
    }

    if(length(mask_cluster)!=ncol(rna_tf.mtx)){print("Alarm!")}
    
    x <- rna_tf.mtx[j,mask_cluster]
    y <- rna_targets.mtx[i,mask_cluster]
    
    if (args$plot_correlations & ((i%in%c(markers_TF_undupl$TF,sign_senders)) & (j%in%c(markers_TF_undupl$TF,sign_senders)))){
      # print(paste0(j," vs ",i))
      to.plot <- data.frame(tf=rep(j,length(rna_tf.mtx[j,])),gene=rep(i,length(rna_tf.mtx[j])),
                            x=rna_tf.mtx[j,],y=rna_targets.mtx[i,],
                            sample=factor(sce_mc$cluster))
      # print(head(to.plot))
      
      p1 <- ggscatter(to.plot, x="x", y="y", fill="sample", size=5, shape=21)  +
              scale_fill_manual(values=opts$celltype.colors[1:length(unique(to.plot$sample))]) +
              labs(x=sprintf("%s expression",j), y=sprintf("%s expression",i)) +
              theme_classic() +
              theme(axis.text = element_text(size=rel(0.7))) +
              geom_smooth(data = to.plot[mask_cluster,], mapping=aes(x=x,y=y), method = "lm",  inherit.aes = FALSE) +
              stat_regline_equation(data = to.plot[mask_cluster,], mapping=aes(x=x,y=y), label.x.npc = 0.5, label.y.npc = 0.1,  inherit.aes = FALSE) +
              stat_cor(data = to.plot[mask_cluster,], mapping=aes(x=x, y=y), method = "pearson", label.x.npc = 0.5, label.y.npc = 0.05) 
        
      pdf(sprintf("%s/plots_cor_score%s/%s_vs_%s.pdf",args$outdir,args$min_chip_score,j,i), 
          width = 7, height = 5)
      suppressMessages(print(p1))
      dev.off()
    }
    
    lm.fit <- lm(y~x)
    if (is.na(coef(lm.fit)[2])){
      dat_return <- rbind(dat_return,c(j, i, NA, NA, paste(cluster_range,collapse="-")))
    } else {
      dat_return <- rbind(dat_return,c(j, i, 
                                      round(coef(lm.fit)[[2]],3), 
                                      format(summary(lm.fit)$coefficients[2,4], digits=3),
                                      paste(cluster_range,collapse="-")))
    }
  }
  if (is.null(dim(dat_return))){
    return(NULL)
  } else {
    colnames(dat_return) <- c("tf","gene","beta","pvalue","clusters_for_calc")
    dat_return <- data.frame(dat_return,stringsAsFactors = F)
    return(dat_return)
  }
}



##########
## alt7 ##
##########

myfun_alt7 <- function(i){ #alt7
  if(which(genes==i)%%round(length(genes)/10)==0){print(paste0(round(which(genes==i)/length(genes)*100),"% done!"))}
  tfs <- tmp[tmp$gene==i,tf]
  dat_return <- NULL

  clust_i <- gene_cluster_annot$cluster[gene_cluster_annot$gene==i]
  
  # print("Number of NULLs")
  # print(sum(is.null(tfs)))
  for (j in tfs){
    clust_j <- gene_cluster_annot$cluster[gene_cluster_annot$gene==j]

    if (j %in% sign_senders){
      cluster_range <- max(c(clust_i,clust_j))
      cluster_range <- 1:cluster_range
    } else {
      if (clust_j>=clust_i){
        cluster_range <- sort(c(clust_i,clust_j))
        cluster_range <- cluster_range[1]:cluster_range[2]
      } else {
        cluster_range <- max(c(clust_i,clust_j))
        cluster_range <- 1:cluster_range
      }
    } 

    if (sum(gene_cluster_annot$cluster==5)==0){
      cluster_sce_mc <- sce_mc$cluster
      cluster_sce_mc[cluster_sce_mc==5] <- 4
      mask_cluster <- cluster_sce_mc%in%cluster_range
    } else {
      mask_cluster <- sce_mc$cluster%in%cluster_range
    }

    if(length(mask_cluster)!=ncol(rna_tf.mtx)){print("Alarm!")}
    
    x <- rna_tf.mtx[j,mask_cluster]
    y <- rna_targets.mtx[i,mask_cluster]
    
    if (args$plot_correlations & ((i%in%c(markers_TF_undupl$TF,sign_senders)) & (j%in%c(markers_TF_undupl$TF,sign_senders)))){
      # print(paste0(j," vs ",i))
      to.plot <- data.frame(tf=rep(j,length(rna_tf.mtx[j,])),gene=rep(i,length(rna_tf.mtx[j])),
                            x=rna_tf.mtx[j,],y=rna_targets.mtx[i,],
                            sample=factor(sce_mc$cluster))
      # print(head(to.plot))
      
      p1 <- ggscatter(to.plot, x="x", y="y", fill="sample", size=5, shape=21)  +
              scale_fill_manual(values=opts$celltype.colors[1:length(unique(to.plot$sample))]) +
              labs(x=sprintf("%s expression",j), y=sprintf("%s expression",i)) +
              theme_classic() +
              theme(axis.text = element_text(size=rel(0.7))) +
              geom_smooth(data = to.plot[mask_cluster,], mapping=aes(x=x,y=y), method = "lm",  inherit.aes = FALSE) +
              stat_regline_equation(data = to.plot[mask_cluster,], mapping=aes(x=x,y=y), label.x.npc = 0.5, label.y.npc = 0.1,  inherit.aes = FALSE) +
              stat_cor(data = to.plot[mask_cluster,], mapping=aes(x=x, y=y), method = "pearson", label.x.npc = 0.5, label.y.npc = 0.05) 
        
      pdf(sprintf("%s/plots_cor_score%s/%s_vs_%s.pdf",args$outdir,args$min_chip_score,j,i), 
          width = 7, height = 5)
      suppressMessages(print(p1))
      dev.off()
    }
    
    lm.fit <- lm(y~x)
    if (is.na(coef(lm.fit)[2])){
      dat_return <- rbind(dat_return,c(j, i, NA, NA, paste(cluster_range,collapse="-")))
    } else {
      dat_return <- rbind(dat_return,c(j, i, 
                                      round(coef(lm.fit)[[2]],3), 
                                      format(summary(lm.fit)$coefficients[2,4], digits=3),
                                      paste(cluster_range,collapse="-")))
    }
  }
  if (is.null(dim(dat_return))){
    return(NULL)
  } else {
    colnames(dat_return) <- c("tf","gene","beta","pvalue","clusters_for_calc")
    dat_return <- data.frame(dat_return,stringsAsFactors = F)
    return(dat_return)
  }
}