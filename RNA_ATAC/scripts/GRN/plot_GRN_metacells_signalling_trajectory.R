#####################################################
##                                                 ##
##  Plot_GRN_metacells_signalling_trajectories.R   ##
##                                                 ##
#####################################################

source("/data/homes/louisc/Project_Babraham/RNA_ATAC/scripts/Settings.R")

suppressMessages(library(GGally))
suppressMessages(library(network))
suppressMessages(library(sna))
suppressMessages(library(ggraph))
suppressMessages(library(igraph))
suppressMessages(library(tidygraph))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--trajectory_name',  type="character",              help='Trajectory name')
p$add_argument('--rna_metacells.sce',  type="character",              help='SingleCellExperiment (metacells)')
p$add_argument('--rna_pseudobulk.sce',       type="character",                help='SingleCellExperiment (pseudobulk)')
p$add_argument('--grn_coef',       type="character",                help='GRN coefficients')
p$add_argument('--min_coef',  type="double",            default=0.25,      help='Minimal coefficient for linear regression')
p$add_argument('--max_pval',  type="double",            default=0.10,      help='Maximal pvalue for linear regression')
p$add_argument('--min_chip_score',  type="double",            default=0.25,      help='Minimal chip score')
p$add_argument('--TFs_to_plot',       type="character",                help='File containing TF of interest to plot')
p$add_argument('--force',       type="logical",  default=F         ,      help='Should replotting be forced?')
p$add_argument('--sign_connections',       type="character",                help='Metadata metacells')
p$add_argument('--outdir',       type="character",                help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

# I/O
# args <- list()
# args$rna_metacells.sce <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/metacells/SingleCellExperiment_metacells_nodiff.rds"
# args$rna_pseudobulk.sce < "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/SingleCellExperiment_pseudobulk.rds"
# args$trajectory_name <- "N2P"
# args$TFs_to_plot <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/chromvar_chip/GRN_TFs_to_plot.txt"
# args$grn_coef <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/gene_regulatory_network/metacells/trajectories/N2P/sign/global_chip_GRN_coef_score0.06.txt.gz"
# args$outdir <-  "/data/homes/louisc/Project_Babraham/RNA_ATAC/gene_regulatory_network/metacells/trajectories/N2P/sign/marker_TF/fig"
# args$min_coef <- 0.25
# args$max_pval <- 0.10
# args$min_chip_score <- 0.06
# args$force <- TRUE
# args$sign_connections <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/NicheNet/GRN_coef_signalling_to_add.txt"

print(args)

#####################
## Define settings ##
#####################

args$trajectory <- sprintf("/data/homes/louisc/Project_Babraham/RNA_ATAC/trajectories/%s/%s_trajectory.txt.gz",
                           args$trajectory_name,args$trajectory_name)

if(!dir.exists(file.path(args$outdir,sprintf("Score%s_Coef%s",args$min_chip_score,args$min_coef)))){
  dir.create(args$outdir, showWarnings = F)
  dir.create(file.path(args$outdir,sprintf("Score%s_Coef%s",args$min_chip_score,args$min_coef)))
}

# Options
opts$celltype_trajectory_dic <- list(
  "N2P" = c(1,2,3,4,5)
)

stopifnot(args$trajectory_name%in%names(opts$celltype_trajectory_dic))
celltypes.to.plot <- opts$celltype_trajectory_dic[[args$trajectory_name]]

##############
## Load TFs ##
##############

# all_TFs
all_TFs <- read.table("/data/homes/louisc/Project_Babraham/TFs.txt",sep="\n",header=F,quote="")$V1
print("All TFs")
print(head(all_TFs))

# marker TFs
markers_TF <- NULL
for (i in celltypes.to.plot){
  marker_tf_i <- read.table(sprintf("/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/RNA/marker_TFs/Marker_TF_cluster%s.txt",i),header=T,sep="\t")
  markers_TF <- rbind(markers_TF,cbind(marker_tf_i$Gene,rep(i,length(marker_tf_i$Gene))))
}
markers_TF <- data.frame(markers_TF)
colnames(markers_TF) <- c("TF","celltype")
markers_TF_undupl <- markers_TF[!duplicated(markers_TF$TF),]
print("Marker TFs")
print(head(markers_TF_undupl))

###########################
## Load signalling & add ##
###########################

sign2TF <- fread(args$sign_connections) 
lines_to_add <- sign2TF[!duplicated(sign2TF$from),c("from","fr_max_expr")]
colnames(lines_to_add) <- c("TF","celltype")
markers_TF_undupl <- rbind(markers_TF_undupl,lines_to_add)

##############################
## Load RNA expression data ##
##############################

# SingleCellExperiment at metacell resolution
rna_metacells.sce <- readRDS(args$rna_metacells.sce)
print("rna_metacells.sce")
rna_metacells.sce

# SingleCellExperiment at pseudobulk resolution
rna_pseudobulk.sce <- readRDS(args$rna_pseudobulk.sce)#[,celltypes.to.plot]
print("rna_pseudobulk.sce")
rna_pseudobulk.sce

##################################
## Load global GRN coefficients ##
##################################

GRN_coef.dt <- fread(args$grn_coef) 
print(head(GRN_coef.dt))
print(dim(GRN_coef.dt))

##################
## Create table ##
##################

if (!file.exists(file.path(args$outdir,sprintf("Score%s_Coef%s/table_interactions.txt",args$min_chip_score,args$min_coef))) | args$force){
  table_interactions <- GRN_coef.dt[,c("sender","target","beta","pvalue")] %>% setnames(c("from","to","beta","pval")) %>% .[!from==to]

  table_interactions$from_cluster <- NA
  table_interactions$to_cluster <- NA
  for (i in 1:nrow(table_interactions)){
    if (sum(markers_TF_undupl$TF==table_interactions$from[i])>0){
      table_interactions$from_cluster[i] <- markers_TF_undupl$celltype[markers_TF_undupl$TF==table_interactions$from[i]]
    }
    if (sum(markers_TF_undupl$TF==table_interactions$to[i])>0){
      table_interactions$to_cluster[i] <- markers_TF_undupl$celltype[markers_TF_undupl$TF==table_interactions$to[i]]
    }
  }
  
  table_interactions$from_supercluster <- c("A","B","C")[(table_interactions$from_cluster>1)+(table_interactions$from_cluster>3)+1]
  table_interactions$to_supercluster <- c("A","B","C")[(table_interactions$to_cluster>1)+(table_interactions$to_cluster>3)+1]
  
  table_interactions$from_cluster_all <- NA
  table_interactions$to_cluster_all <- NA
  for (i in 1:nrow(table_interactions)){
    if (sum(markers_TF$TF==table_interactions$from[i])>0){
      table_interactions$from_cluster_all[i] <- paste(markers_TF$celltype[markers_TF$TF==table_interactions$from[i]],collapse=",")
    }
    if (sum(markers_TF$TF==table_interactions$to[i])>0){
      table_interactions$to_cluster_all[i] <- paste(markers_TF$celltype[markers_TF$TF==table_interactions$to[i]],collapse=",")
    }
  }
  
  print("Table output")
  print(head(table_interactions))
  print("Marker TF interactions")
  print(dim(table_interactions[(!is.na(table_interactions$from_cluster)) & (!is.na(table_interactions$to_cluster)),]))
  
  write.table(table_interactions,file=file.path(args$outdir,sprintf("Score%s_Coef%s/table_interactions.txt",args$min_chip_score,args$min_coef)),
              col.names = T, row.names = F, sep="\t", quote=F)
} else {
  table_interactions <- read.table(file=file.path(args$outdir,sprintf("Score%s_Coef%s/table_interactions.txt",args$min_chip_score,args$min_coef)),
                                   header=T,sep="\t")
  print("Skipped creating Table")
}

# ##
# ## THE REMAINDER OF THIS SCRIPT IS NOT YET ADAPTED TO SIGNALLING PW!!!!! 
# ##

# ##############################
# ## Plot network per cluster ##
# ##############################

# for (j in 1:5){
#   print(paste0("Plot network for cluster ",j))
  
#   # Create node and edge data.frames
#   TFs_marker_j <- TFs_marker[TFs_marker%in%markers_TF_undupl$TF[markers_TF_undupl$celltype==j]]
#   node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
#   edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$gene%in%TFs_marker)][,c("tf","gene","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]
  
#   # Create igraph object
#   igraph.net_marker <- graph_from_data_frame(d = edge_list.dt_marker)
  
#   # Create tbl_graph object for ggraph
#   igraph.tbl_marker <- as_tbl_graph(igraph.net_marker) %>%
#     activate(nodes) %>%
#     mutate(tf=names(V(igraph.net_marker))) %>%
#     mutate(degree=igraph::degree(igraph.net_marker)) %>%
#     mutate(eigen_centrality=eigen_centrality(igraph.net_marker)$vector) %>%
#     activate(edges) %>%
#     mutate(sign=ifelse(E(igraph.net_marker)$weight>0,"Positive","Negative"))
  
#   celltypes_markers <- rep(NA,length(names(V(igraph.net_marker))))
#   for (i in 1:length(celltypes_markers)){
#     if (sum(markers_TF_undupl$TF==names(V(igraph.net_marker))[i])>0){
#       celltypes_markers[i] <- markers_TF_undupl$celltype[markers_TF_undupl$TF==names(V(igraph.net_marker))[i]]
#     }
#   }

#   igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(celltype=celltypes_markers)
  
#   celltypes_present <- sort(as.numeric(unique(celltypes_markers)))
#   celltypes_present <- celltypes_present[!is.na(celltypes_present)]
#   print(celltypes_present)
  
#   if(length(celltypes_present)==0){
#     print(paste0("No marker TF left for cluster ",j," with coefficient cutoff at ",args$min_coef))
#     next()
#   }
  
#   set.seed(42)
#   p <- ggraph(igraph.tbl_marker, 'stress') +
#     geom_edge_link(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
#                    arrow=arrow(length=unit(1.5,'mm')), end_cap=circle(1,'mm')) +
#     geom_node_point(aes(fill=celltype), size=3, stroke=0, shape=21, alpha=0.75) +
#     geom_node_text(aes(label=name), size=1) +
#     scale_edge_colour_manual(values=c("Positive"="darkred", "Negative"="darkblue")) +
#     scale_fill_manual(values=opts$celltype.colors[celltypes_present],na.value="white") +
#     theme_graph(base_family = 'Helvetica')
  
#   pdf(file.path(args$outdir,paste0("Score",args$min_chip_score,"_Coef",args$min_coef,"/cluster",j,"_celltype_markerTFs.pdf")),
#       width = 5.5, height = 5.75)
#   print(p)
#   dev.off()
  
#   # Create node and edge data.frames
#   node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
#   edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$gene%in%all_TFs)][,c("tf","gene","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]
  
#   # Create igraph object
#   igraph.net_marker <- graph_from_data_frame(d = edge_list.dt_marker)
  
#   # Create tbl_graph object for ggraph
#   igraph.tbl_marker <- as_tbl_graph(igraph.net_marker) %>%
#     activate(nodes) %>%
#     mutate(tf=names(V(igraph.net_marker))) %>%
#     mutate(degree=igraph::degree(igraph.net_marker)) %>%
#     mutate(eigen_centrality=eigen_centrality(igraph.net_marker)$vector) %>%
#     activate(edges) %>%
#     mutate(sign=ifelse(E(igraph.net_marker)$weight>0,"Positive","Negative"))
#   # igraph.tbl
  
#   celltypes_markers <- rep(NA,length(names(V(igraph.net_marker))))
#   for (i in 1:length(celltypes_markers)){
#     if (sum(markers_TF_undupl$TF==names(V(igraph.net_marker))[i])>0){
#       celltypes_markers[i] <- markers_TF_undupl$celltype[markers_TF_undupl$TF==names(V(igraph.net_marker))[i]]
#     }
#   }

#   igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(celltype=celltypes_markers)
  
#   celltypes_present <- sort(as.numeric(unique(celltypes_markers)))
#   celltypes_present <- celltypes_present[!is.na(celltypes_present)]
#   print(celltypes_present)
  
#   set.seed(42)
#   p <- ggraph(igraph.tbl_marker, 'stress') +
#     geom_edge_link(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
#                    arrow=arrow(length=unit(1.5,'mm')), end_cap=circle(1,'mm')) +
#     geom_node_point(aes(fill=celltype), size=3, stroke=0, shape=21, alpha=0.75) +
#     geom_node_text(aes(label=name), size=1) +
#     scale_edge_colour_manual(values=c("Positive"="darkred", "Negative"="darkblue")) +
#     scale_fill_manual(values=opts$celltype.colors[celltypes_present],na.value="white") +
#     theme_graph(base_family = 'Helvetica')
  
#   pdf(file.path(args$outdir,paste0("Score",args$min_chip_score,"_Coef",args$min_coef,"/cluster",j,"_celltype_allTFs.pdf")),
#       width = 5.5, height = 5.75)
#   print(p)
#   dev.off()
# }

# ##############################################
# ## Plot network per cluster on general plot ##
# ##############################################

# for (j in 1:5){
#   for (k in 1:2){
#     if (k == 1){
#       interactions <- "activatory"
#     } else {
#       interactions <- "repressive"
#     }
    
#     print(paste0("Plot global network for cluster ",j))
    
#     # Create node and edge data.frames
#     TFs_marker_j <- TFs_marker[TFs_marker%in%markers_TF_undupl$TF[markers_TF_undupl$celltype==j]]
#     node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
#     edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$tf%in%TFs_marker) & (GRN_coef.dt$gene%in%TFs_marker)][,c("tf","gene","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]
    
#     # Create igraph object
#     igraph.net_marker <- graph_from_data_frame(d = edge_list.dt_marker)
    
#     # Create tbl_graph object for ggraph
#     igraph.tbl_marker <- as_tbl_graph(igraph.net_marker) %>%
#       activate(nodes) %>%
#       mutate(tf=names(V(igraph.net_marker))) %>%
#       mutate(degree=igraph::degree(igraph.net_marker)) %>%
#       mutate(eigen_centrality=eigen_centrality(igraph.net_marker)$vector) %>%
#       activate(edges) %>%
#       mutate(sign=ifelse(E(igraph.net_marker)$weight>0,"Positive","Negative"))
#     # print(length(E(igraph.net_marker)$weight))
#     mask_not <- (!(edge_list.dt_marker$from%in%TFs_marker_j))
#     mask_pos <- (E(igraph.tbl_marker)$weight>0) 
#     mask_neg <- (E(igraph.tbl_marker)$weight<0) 
#     # print(sum(mask_pos & mask_not))
#     # print(sum(mask_neg & mask_not))
#     E(igraph.tbl_marker)$sign[mask_not] <- paste0("Not in cluster (n=",sum(mask_not),")")
#     E(igraph.tbl_marker)$sign[mask_pos & !mask_not] <- paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")")
#     E(igraph.tbl_marker)$sign[mask_neg & !mask_not] <- paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")")
#     E(igraph.tbl_marker)$sign <- factor(E(igraph.tbl_marker)$sign,levels=c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
#                                                                            paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
#                                                                            paste0("Not in cluster (n=",sum(mask_not),")")))
    
#     labs_plot <- c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
#                    paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
#                    paste0("Not in cluster (n=",sum(mask_not),")"))
    
#     print(labs_plot)
    
#     celltypes_markers <- rep(NA,length(V(igraph.tbl_marker)$tf))
#     for (i in 1:length(celltypes_markers)){
#       if (sum(markers_TF_undupl$TF==V(igraph.tbl_marker)$tf[i])>0){
#         celltypes_markers[i] <- markers_TF_undupl$celltype[markers_TF_undupl$TF==V(igraph.tbl_marker)$tf[i]]
#       }
#     }

#     if (k == 1){
#       celltypes_markers[!((celltypes_markers==j) | (V(igraph.tbl_marker)$tf%in%GRN_coef.dt$gene[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$beta>0)]))] <- NA
#     } else {
#       celltypes_markers[!((celltypes_markers==j) | (V(igraph.tbl_marker)$tf%in%GRN_coef.dt$gene[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$beta<0)]))] <- NA
#     }
    
#     igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(celltype=celltypes_markers)

#     if (k == 1){
#       col_arrows <- c("darkred","gray","gray")
#     } else{
#       col_arrows <- c("gray","darkblue","gray")
#     }
    
#     celltypes_present <- sort(as.numeric(unique(celltypes_markers)))
#     celltypes_present <- celltypes_present[!is.na(celltypes_present)]
#     print(celltypes_present)
    
#     if(length(celltypes_present)==0){
#       print(paste0("No marker TF left for cluster ",j," with coefficient cutoff at ",args$min_coef))
#       next()
#     }
     
#     set.seed(42)
#     p <- ggraph(igraph.tbl_marker, 'stress') +
#       geom_edge_link(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
#                      arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(2,'mm'), start_cap=circle(2,'mm')) +
#       geom_node_point(aes(fill=celltype), size=5, stroke=0, shape=21, alpha=0.75) +
#       geom_node_text(aes(label=name), size=1.5,) +
#       scale_edge_colour_manual(values=col_arrows,
#                                breaks=waiver(),
#                                labels=labs_plot) +
#       scale_fill_manual(values=opts$celltype.colors[celltypes_present],na.value="gray") +
#       theme_graph(base_family = 'Helvetica')
    
#     if (k == 1){
#       pdf(file.path(args$outdir,sprintf("Score%s_Coef%s/global_network_markerTFs_cluster%s_%s.pdf",args$min_chip_score,args$min_coef,j,interactions)), 
#           width = 5.5, height = 5.75)
#     } else {
#       pdf(file.path(args$outdir,sprintf("Score%s_Coef%s/global_network_markerTFs_cluster%s_%s.pdf",args$min_chip_score,args$min_coef,j,interactions)), 
#           width = 5.5, height = 5.75)
#     }
#     print(p)
#     dev.off()
#   }
# }


# ##################################################################
# ## Plot network per cluster on general plot (no double markers) ##
# ##################################################################

# for (j in 1:5){
#   for (k in 1:2){
#     if (k == 1){
#       interactions <- "activatory"
#     } else {
#       interactions <- "repressive"
#     }
    
#     print(paste0("Plot global network for cluster ",j))
    
#     # Create node and edge data.frames
#     TFs_marker_j <- TFs_marker[TFs_marker%in%table_interactions$from[!is.na(table_interactions$from_cluster_all) & (table_interactions$from_cluster_all==j) & !grepl(",",table_interactions$from_cluster_all)]]
#     print(TFs_marker_j)
#     node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
#     edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$tf%in%TFs_marker) & (GRN_coef.dt$gene%in%TFs_marker)][,c("tf","gene","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]
    
#     # Create igraph object
#     igraph.net_marker <- graph_from_data_frame(d = edge_list.dt_marker)
    
#     # Create tbl_graph object for ggraph
#     igraph.tbl_marker <- as_tbl_graph(igraph.net_marker) %>%
#       activate(nodes) %>%
#       mutate(tf=names(V(igraph.net_marker))) %>%
#       mutate(degree=igraph::degree(igraph.net_marker)) %>%
#       mutate(eigen_centrality=eigen_centrality(igraph.net_marker)$vector) %>%
#       activate(edges) %>%
#       mutate(sign=ifelse(E(igraph.net_marker)$weight>0,"Positive","Negative"))
#     # print(length(E(igraph.net_marker)$weight))
#     mask_pos <- (E(igraph.tbl_marker)$weight>0) 
#     mask_not <- (!(edge_list.dt_marker$from%in%TFs_marker_j))
#     mask_neg <- (E(igraph.tbl_marker)$weight<0)
#     # print(sum(mask_pos & mask_not))
#     # print(sum(mask_neg & mask_not))
#     E(igraph.tbl_marker)$sign[mask_not] <- paste0("Not in cluster (n=",sum(mask_not),")")
#     E(igraph.tbl_marker)$sign[mask_pos & !mask_not] <- paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")")
#     E(igraph.tbl_marker)$sign[mask_neg & !mask_not] <- paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")")
#     E(igraph.tbl_marker)$sign <- factor(E(igraph.tbl_marker)$sign,levels=c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
#                                                                            paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
#                                                                            paste0("Not in cluster (n=",sum(mask_not),")")))
    
#     labs_plot <- c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
#                    paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
#                    paste0("Not in cluster (n=",sum(mask_not),")"))
    
#     print(labs_plot)
    
#     celltypes_markers <- rep(NA,length(V(igraph.tbl_marker)$tf))
#     for (i in 1:length(celltypes_markers)){
#       if (sum(markers_TF_undupl$TF==V(igraph.tbl_marker)$tf[i])>0){
#         celltypes_markers[i] <- markers_TF_undupl$celltype[markers_TF_undupl$TF==V(igraph.tbl_marker)$tf[i]]
#       }
#     }
    
#     if (k == 1){
#       celltypes_markers[!((V(igraph.tbl_marker)$tf%in%table_interactions$from[!is.na(table_interactions$from_cluster_all) & (table_interactions$from_cluster_all==j) & !grepl(",",table_interactions$from_cluster_all) & (table_interactions$beta>0)]))] <- NA
#     } else {
#       celltypes_markers[!((V(igraph.tbl_marker)$tf%in%table_interactions$from[!is.na(table_interactions$from_cluster_all) & (table_interactions$from_cluster_all==j) & !grepl(",",table_interactions$from_cluster_all) & (table_interactions$beta<0)]))] <- NA
#     }
    
#     igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(celltype=celltypes_markers)
    
#     if (k == 1){
#       col_arrows <- c("darkred","gray","gray")
#     } else{
#       col_arrows <- c("gray","darkblue","gray")
#     }
    
#     celltypes_present <- sort(as.numeric(unique(celltypes_markers)))
#     celltypes_present <- celltypes_present[!is.na(celltypes_present)]
#     print(celltypes_present)
    
#     if(length(celltypes_present)==0){
#       print(paste0("No marker TF left for cluster ",j," with coefficient cutoff at ",args$min_coef))
#       next()
#     }
    
#     set.seed(42)
#     p <- ggraph(igraph.tbl_marker, 'stress') +
#       geom_edge_link(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
#                      arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(2,'mm'), start_cap=circle(2,'mm')) +
#       geom_node_point(aes(fill=celltype), size=5, stroke=0, shape=21, alpha=0.75) +
#       geom_node_text(aes(label=name), size=1.5,) +
#       scale_edge_colour_manual(values=col_arrows,
#                                breaks=waiver(),
#                                labels=labs_plot) +
#       scale_fill_manual(values=opts$celltype.colors[celltypes_present],na.value="gray") +
#       theme_graph(base_family = 'Helvetica')
    
#     if (k == 1){
#       pdf(file.path(args$outdir,sprintf("Score%s_Coef%s/global_network_markerTFs_nodupl_cluster%s_%s.pdf",args$min_chip_score,args$min_coef,j,interactions)), 
#           width = 5.5, height = 5.75)
#     } else {
#       pdf(file.path(args$outdir,sprintf("Score%s_Coef%s/global_network_markerTFs_nodupl_cluster%s_%s.pdf",args$min_chip_score,args$min_coef,j,interactions)), 
#           width = 5.5, height = 5.75)
#     }
#     print(p)
#     dev.off()
#   }
# }

# #################################
# ## Plot network per transition ##
# #################################

# for (j in 1:4){
#   print(paste0("Plot network for transition from ",j, " to ",j+1))
  
#   # Create node and edge data.frames
#   TFs_marker_j <- TFs_marker[TFs_marker%in%markers_TF_undupl$TF[markers_TF_undupl$celltype==j]]
#   TFs_marker_j1 <- TFs_marker[TFs_marker%in%markers_TF_undupl$TF[markers_TF_undupl$celltype==(j+1)]]
#   node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
#   edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$gene%in%TFs_marker_j1)][,c("tf","gene","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]
  
#   # Create igraph object
#   igraph.net_marker <- graph_from_data_frame(d = edge_list.dt_marker)
  
#   # Create tbl_graph object for ggraph
#   igraph.tbl_marker <- as_tbl_graph(igraph.net_marker) %>%
#     activate(nodes) %>%
#     mutate(tf=names(V(igraph.net_marker))) %>%
#     mutate(degree=igraph::degree(igraph.net_marker)) %>%
#     mutate(eigen_centrality=eigen_centrality(igraph.net_marker)$vector) %>%
#     activate(edges) %>%
#     mutate(sign=ifelse(E(igraph.net_marker)$weight>0,"Positive","Negative"))
#   # igraph.tbl
  
#   celltypes_markers <- rep(NA,length(names(V(igraph.net_marker))))
#   for (i in 1:length(celltypes_markers)){
#     if (sum(markers_TF_undupl$TF==names(V(igraph.net_marker))[i])>0){
#       celltypes_markers[i] <- markers_TF_undupl$celltype[markers_TF_undupl$TF==names(V(igraph.net_marker))[i]]
#     }
#   }

#   igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(celltype=celltypes_markers)
  
#   celltypes_present <- sort(as.numeric(unique(celltypes_markers)))
#   celltypes_present <- celltypes_present[!is.na(celltypes_present)]
#   print(celltypes_present)
  
#   if(length(celltypes_present)==0){
#     print(paste0("No marker TF left for cluster ",j," with coefficient cutoff at ",args$min_coef))
#     next()
#   }
  
#   set.seed(42)
#   p <- ggraph(igraph.tbl_marker, 'stress') +
#     geom_edge_link(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
#                    arrow=arrow(length=unit(1.5,'mm')), end_cap=circle(1,'mm')) +
#     geom_node_point(aes(fill=celltype), size=3, stroke=0, shape=21, alpha=0.75) +
#     geom_node_text(aes(label=name), size=1) +
#     scale_edge_colour_manual(values=c("Positive"="darkred", "Negative"="darkblue")) +
#     scale_fill_manual(values=opts$celltype.colors[celltypes_present],na.value="white") +
#     theme_graph(base_family = 'Helvetica')
  
#   pdf(file.path(args$outdir,paste0("Score",args$min_chip_score,"_Coef",args$min_coef,"/transition",j,"to",j+1,"_celltype_markerTFs.pdf")),
#       width = 5.5, height = 5.75)
#   print(p)
#   dev.off()
# }

# ###################################
# ## Plot network per supercluster ##
# ###################################

# for (j in c("A","B","C")){
#   print(paste0("Plot network for supercluster ",j))
  
#   # Create node and edge data.frames
#   if (j=="A"){
#     TFs_marker_j <- TFs_marker[TFs_marker%in%markers_TF_undupl$TF[markers_TF_undupl$celltype==1]]
#   } else if (j=="B"){
#     TFs_marker_j <- TFs_marker[TFs_marker%in%markers_TF_undupl$TF[markers_TF_undupl$celltype%in%c(2,3)]]
#   } else if (j=="C"){
#     TFs_marker_j <- TFs_marker[TFs_marker%in%markers_TF_undupl$TF[markers_TF_undupl$celltype%in%c(4,5)]]
#   }
#   node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
#   edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$gene%in%TFs_marker)][,c("tf","gene","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]
  
#   # Create igraph object
#   igraph.net_marker <- graph_from_data_frame(d = edge_list.dt_marker)
  
#   # Create tbl_graph object for ggraph
#   igraph.tbl_marker <- as_tbl_graph(igraph.net_marker) %>%
#     activate(nodes) %>%
#     mutate(tf=names(V(igraph.net_marker))) %>%
#     mutate(degree=igraph::degree(igraph.net_marker)) %>%
#     mutate(eigen_centrality=eigen_centrality(igraph.net_marker)$vector) %>%
#     activate(edges) %>%
#     mutate(sign=ifelse(E(igraph.net_marker)$weight>0,"Positive","Negative"))
#   # igraph.tbl
  
#   celltypes_markers <- rep(NA,length(names(V(igraph.net_marker))))
#   for (i in 1:length(celltypes_markers)){
#     if (sum(markers_TF_undupl$TF==names(V(igraph.net_marker))[i])>0){
#       celltypes_markers[i] <- markers_TF_undupl$celltype[markers_TF_undupl$TF==names(V(igraph.net_marker))[i]]
#     }
#   }
  
#   celltypes_present <- sort(as.numeric(unique(celltypes_markers)))
#   celltypes_present <- celltypes_present[!is.na(celltypes_present)]
#   celltypes_present <- unique(c(1,3,3,5,5)[celltypes_present])
#   print(celltypes_present)
  
#   if(length(celltypes_present)==0){
#     print(paste0("No marker TF left for cluster ",j," with coefficient cutoff at ",args$min_coef))
#     next()
#   }
  
#   celltypes_markers[celltypes_markers==1] <- "A"
#   celltypes_markers[celltypes_markers%in%c(2,3)] <- "B"
#   celltypes_markers[celltypes_markers%in%c(4,5)] <- "C" 
  
#   igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(celltype=celltypes_markers)
  
#   set.seed(42)
#   p <- ggraph(igraph.tbl_marker, 'stress') +
#     geom_edge_link(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
#                    arrow=arrow(length=unit(1.5,'mm')), end_cap=circle(1,'mm')) +
#     geom_node_point(aes(fill=celltype), size=3, stroke=0, shape=21, alpha=0.75) +
#     geom_node_text(aes(label=name), size=1) +
#     scale_edge_colour_manual(values=c("Positive"="darkred", "Negative"="darkblue")) +
#     scale_fill_manual(values=opts$celltype.colors[celltypes_present],na.value="white") +
#     theme_graph(base_family = 'Helvetica')
  
#   pdf(file.path(args$outdir,paste0("Score",args$min_chip_score,"_Coef",args$min_coef,"/supercluster",j,"_celltype_markerTFs.pdf")),
#       width = 5.5, height = 5.75)
#   print(p)
#   dev.off()
# }

# ###################################################
# ## Plot network per supercluster on general plot ##
# ###################################################

# for (j in c("A","B","C")){
#   for (k in 1:2){
#     if (k == 1){
#       interactions <- "activatory"
#     } else {
#       interactions <- "repressive"
#     }
    
#     print(paste0("Plot global network for supercluster ",j))
    
#     # Create node and edge data.frames
#     if (j=="A"){
#       TFs_marker_j <- TFs_marker[TFs_marker%in%markers_TF_undupl$TF[markers_TF_undupl$celltype==1]]
#     } else if (j=="B"){
#       TFs_marker_j <- TFs_marker[TFs_marker%in%markers_TF_undupl$TF[markers_TF_undupl$celltype%in%c(2,3)]]
#     } else if (j=="C"){
#       TFs_marker_j <- TFs_marker[TFs_marker%in%markers_TF_undupl$TF[markers_TF_undupl$celltype%in%c(4,5)]]
#     }
#     node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
#     edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$tf%in%TFs_marker) & (GRN_coef.dt$gene%in%TFs_marker)][,c("tf","gene","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]
    
#     # Create igraph object
#     igraph.net_marker <- graph_from_data_frame(d = edge_list.dt_marker)
    
#     # Create tbl_graph object for ggraph
#     igraph.tbl_marker <- as_tbl_graph(igraph.net_marker) %>%
#       activate(nodes) %>%
#       mutate(tf=names(V(igraph.net_marker))) %>%
#       mutate(degree=igraph::degree(igraph.net_marker)) %>%
#       mutate(eigen_centrality=eigen_centrality(igraph.net_marker)$vector) %>%
#       activate(edges) %>%
#       mutate(sign=ifelse(E(igraph.net_marker)$weight>0,"Positive","Negative"))
#     # print(length(E(igraph.net_marker)$weight))
#     mask_pos <- (E(igraph.tbl_marker)$weight>0) 
#     mask_not <- (!(edge_list.dt_marker$from%in%TFs_marker_j))
#     mask_neg <- (E(igraph.tbl_marker)$weight<0) 
#     # print(sum(mask_pos & mask_not))
#     # print(sum(mask_neg & mask_not))
#     E(igraph.tbl_marker)$sign[mask_not] <- paste0("Not in cluster (n=",sum(mask_not),")")
#     E(igraph.tbl_marker)$sign[mask_pos & !mask_not] <- paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")")
#     E(igraph.tbl_marker)$sign[mask_neg & !mask_not] <- paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")")
#     E(igraph.tbl_marker)$sign <- factor(E(igraph.tbl_marker)$sign,levels=c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
#                                                                            paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
#                                                                            paste0("Not in cluster (n=",sum(mask_not),")")))
    
#     labs_plot <- c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
#                    paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
#                    paste0("Not in cluster (n=",sum(mask_not),")"))
    
#     print(labs_plot)
#     celltypes_markers <- rep(NA,length(V(igraph.tbl_marker)$tf))
#     for (i in 1:length(celltypes_markers)){
#       if (sum(markers_TF_undupl$TF==V(igraph.tbl_marker)$tf[i])>0){
#         celltypes_markers[i] <- markers_TF_undupl$celltype[markers_TF_undupl$TF==V(igraph.tbl_marker)$tf[i]]
#       }
#     }
    
#     print(celltypes_markers)
#     print(length(celltypes_markers))
    
#     if (k == 1){
#       if (j == "A"){
#         celltypes_markers[!((celltypes_markers==1) | (V(igraph.tbl_marker)$tf%in%GRN_coef.dt$gene[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$beta>0)]))] <- NA
#       } else if (j == "B"){
#         celltypes_markers[!((celltypes_markers%in%c(2,3)) | (V(igraph.tbl_marker)$tf%in%GRN_coef.dt$gene[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$beta>0)]))] <- NA
#       } else if (j == "C"){
#         celltypes_markers[!((celltypes_markers%in%c(4,5)) | (V(igraph.tbl_marker)$tf%in%GRN_coef.dt$gene[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$beta>0)]))] <- NA
#       }
#     } else {
#       if (j == "A"){
#         celltypes_markers[!((celltypes_markers==1) | (V(igraph.tbl_marker)$tf%in%GRN_coef.dt$gene[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$beta<0)]))] <- NA
#       } else if (j == "B"){
#         celltypes_markers[!((celltypes_markers%in%c(2,3)) | (V(igraph.tbl_marker)$tf%in%GRN_coef.dt$gene[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$beta<0)]))] <- NA
#       } else if (j == "C"){
#         celltypes_markers[!((celltypes_markers%in%c(4,5)) | (V(igraph.tbl_marker)$tf%in%GRN_coef.dt$gene[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$beta<0)]))] <- NA
#       }
#     }

#     celltypes_present <- sort(as.numeric(unique(celltypes_markers)))
#     celltypes_present <- celltypes_present[!is.na(celltypes_present)]
#     celltypes_present <- unique(c(1,3,3,5,5)[celltypes_present])
#     print(celltypes_present)
    
#     if(length(celltypes_present)==0){
#       print(paste0("No marker TF left for cluster ",j," with coefficient cutoff at ",args$min_coef))
#       next()
#     }
    
#     celltypes_markers[celltypes_markers==1] <- "A"
#     celltypes_markers[celltypes_markers%in%c(2,3)] <- "B"
#     celltypes_markers[celltypes_markers%in%c(4,5)] <- "C" 
    
#     print(celltypes_markers)
#     print(length(celltypes_markers))
    
#     igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(celltype=celltypes_markers)
    
#     if (k == 1){
#       col_arrows <- c("darkred","gray","gray")
#     } else{
#       col_arrows <- c("gray","darkblue","gray")
#     }
    
#     set.seed(42)
#     p <- ggraph(igraph.tbl_marker, 'stress') +
#       geom_edge_link(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
#                      arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(2,'mm'), start_cap=circle(2,'mm')) +
#       geom_node_point(aes(fill=celltype), size=5, stroke=0, shape=21, alpha=0.75) +
#       geom_node_text(aes(label=name), size=1.5,) +
#       scale_edge_colour_manual(values=col_arrows,
#                                breaks=waiver(),
#                                labels=labs_plot) +
#       scale_fill_manual(values=opts$celltype.colors[celltypes_present],na.value="gray") +
#       theme_graph(base_family = 'Helvetica')
    
#     if (k == 1){
#       pdf(file.path(args$outdir,sprintf("Score%s_Coef%s/global_network_markerTFs_supercluster%s_%s.pdf",args$min_chip_score,args$min_coef,j,interactions)), 
#           width = 5.5, height = 5.75)
#     } else {
#       pdf(file.path(args$outdir,sprintf("Score%s_Coef%s/global_network_markerTFs_supercluster%s_%s.pdf",args$min_chip_score,args$min_coef,j,interactions)), 
#           width = 5.5, height = 5.75)
#     }
#     print(p)
#     dev.off()
#   }
# }

# ######################################
# ## Plot network per supertransition ##
# ######################################

# for (j in 1:2){
#   print(paste0("Plot network for transition from ",c("A","B","C")[j], " to ",c("A","B","C")[j+1]))
  
#   # Create node and edge data.frames
#   if (j==1){
#     TFs_marker_j <- TFs_marker[TFs_marker%in%markers_TF_undupl$TF[markers_TF_undupl$celltype==1]]
#     TFs_marker_j1 <- TFs_marker[TFs_marker%in%markers_TF_undupl$TF[markers_TF_undupl$celltype%in%c(2,3)]]
#   } else if (j==2){
#     TFs_marker_j <- TFs_marker[TFs_marker%in%markers_TF_undupl$TF[markers_TF_undupl$celltype%in%c(2,3)]]
#     TFs_marker_j1 <- TFs_marker[TFs_marker%in%markers_TF_undupl$TF[markers_TF_undupl$celltype%in%c(4,5)]]
#   }
#   node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
#   edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$gene%in%TFs_marker_j1)][,c("tf","gene","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]
  
#   # Create igraph object
#   igraph.net_marker <- graph_from_data_frame(d = edge_list.dt_marker)
  
#   # Create tbl_graph object for ggraph
#   igraph.tbl_marker <- as_tbl_graph(igraph.net_marker) %>%
#     activate(nodes) %>%
#     mutate(tf=names(V(igraph.net_marker))) %>%
#     mutate(degree=igraph::degree(igraph.net_marker)) %>%
#     mutate(eigen_centrality=eigen_centrality(igraph.net_marker)$vector) %>%
#     activate(edges) %>%
#     mutate(sign=ifelse(E(igraph.net_marker)$weight>0,"Positive","Negative"))
#   # igraph.tbl
  
#   celltypes_markers <- rep(NA,length(names(V(igraph.net_marker))))
#   for (i in 1:length(celltypes_markers)){
#     if (sum(markers_TF_undupl$TF==names(V(igraph.net_marker))[i])>0){
#       celltypes_markers[i] <- markers_TF_undupl$celltype[markers_TF_undupl$TF==names(V(igraph.net_marker))[i]]
#     }
#   }

#   celltypes_present <- sort(as.numeric(unique(celltypes_markers)))
#   celltypes_present <- celltypes_present[!is.na(celltypes_present)]
#   celltypes_present <- unique(c(1,3,3,5,5)[celltypes_present])
#   print(celltypes_present)
  
#   if(length(celltypes_present)==0){
#     print(paste0("No marker TF left for cluster ",j," with coefficient cutoff at ",args$min_coef))
#     next()
#   }
  
#   celltypes_markers[celltypes_markers==1] <- "A"
#   celltypes_markers[celltypes_markers%in%c(2,3)] <- "B"
#   celltypes_markers[celltypes_markers%in%c(4,5)] <- "C" 
  
#   igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(celltype=celltypes_markers)
  
#   set.seed(42)
#   p <- ggraph(igraph.tbl_marker, 'stress') +
#     geom_edge_link(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
#                    arrow=arrow(length=unit(1.5,'mm')), end_cap=circle(1,'mm')) +
#     geom_node_point(aes(fill=celltype), size=3, stroke=0, shape=21, alpha=0.75) +
#     geom_node_text(aes(label=name), size=1) +
#     scale_edge_colour_manual(values=c("Positive"="darkred", "Negative"="darkblue")) +
#     scale_fill_manual(values=opts$celltype.colors[celltypes_present],na.value="white") +
#     theme_graph(base_family = 'Helvetica')
  
#   pdf(file.path(args$outdir,paste0("Score",args$min_chip_score,"_Coef",args$min_coef,"/transition",c("A","B","C")[j],"to",c("A","B","C")[j+1],"_celltype_markerTFs.pdf")),
#       width = 5.5, height = 5.75)
#   print(p)
#   dev.off()
# }

# #########################
# ## Plot global network ##
# #########################

# print("Plot global network: celltype")

# # Create node and edge data.frames
# TFs <- unique(c(GRN_coef.dt$tf,GRN_coef.dt$gene))
# node_list.dt <- data.table(node_id=1:length(TFs), node_name=TFs)
# edge_list.dt <- GRN_coef.dt[,c("tf","gene","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]

# print("Summary edges:")
# n_pos_edges <- sum(edge_list.dt$weight>0)
# n_neg_edges <- sum(edge_list.dt$weight<0)
# print(paste0("positive: ",n_pos_edges))
# print(paste0("negative: ",n_neg_edges))

# # Create igraph object
# igraph.net <- graph_from_data_frame(d = edge_list.dt)

# # Create tbl_graph object for ggraph
# igraph.tbl <- as_tbl_graph(igraph.net) %>%
#   activate(nodes) %>%
#   mutate(tf=names(V(igraph.net))) %>%
#   mutate(degree=igraph::degree(igraph.net)) %>%
#   mutate(eigen_centrality=eigen_centrality(igraph.net)$vector) %>%
#   activate(edges) %>%
#   mutate(sign=ifelse(E(igraph.net)$weight>0,paste0("Positive (n=",n_pos_edges,")"),paste0("Negative (n=",n_pos_edges,")")))

# # # Define node color
# # tmp <- marker_TFs_all.dt[gene%in%names(V(igraph.net))] %>% .[,.SD[which.max(score)][,"celltype"],by="gene"] %>% setkey(gene)
# # igraph.tbl <- igraph.tbl %>% activate(nodes) %>% mutate(celltype=tmp[names(V(igraph.net)),celltype])
# # print(names(V(igraph.net)))

# celltypes_markers <- rep(NA,length(names(V(igraph.net))))
# for (i in 1:length(celltypes_markers)){
#   if (sum(markers_TF_undupl$TF==names(V(igraph.net))[i])>0){
#     celltypes_markers[i] <- markers_TF_undupl$celltype[markers_TF_undupl$TF==names(V(igraph.net))[i]]
#   }
# }
# celltypes_markers <- factor(as.character(celltypes_markers),levels=as.character(1:5))

# igraph.tbl <- igraph.tbl %>% activate(nodes) %>% mutate(celltype=celltypes_markers)

# if(!file.exists(file.path(args$outdir,sprintf("Score%s_Coef%s/global_network_celltype.pdf",args$min_chip_score,args$min_coef))) | args$force){
#   # Plot
#   set.seed(42)
#   p <- ggraph(igraph.tbl, 'stress') +
#     geom_edge_link(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
#                    arrow=arrow(length=unit(1.5,'mm')), end_cap=circle(4,'mm')) +
#     geom_node_point(aes(fill=celltype), size=6, stroke=0, shape=21, alpha=0.75) +
#     geom_node_text(aes(label=name), size=2) +
#     scale_edge_colour_manual(values=c("darkblue","darkred")) +
#     scale_fill_manual(values=opts$celltype.colors,na.value="white") +
#     theme_graph(base_family = 'Helvetica')
  
#   pdf(file.path(args$outdir,sprintf("Score%s_Coef%s/global_network_celltype.pdf",args$min_chip_score,args$min_coef)), width = 20, height = 20)
#   print(p)
#   dev.off()
# } else {
#   print("Skipped actual plotting")
# }

# print("Plot global network: celltype (all TFs)")

# # Create node and edge data.frames
# node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
# edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$tf%in%all_TFs) & (GRN_coef.dt$gene%in%all_TFs)][,c("tf","gene","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]

# print("Summary edges:")
# print(paste0("positive: ",sum(edge_list.dt_marker$weight>0)))
# print(paste0("negative: ",sum(edge_list.dt_marker$weight<0)))

# # Create igraph object
# igraph.net_marker <- graph_from_data_frame(d = edge_list.dt_marker)

# # Create tbl_graph object for ggraph
# igraph.tbl_marker <- as_tbl_graph(igraph.net_marker) %>%
#   activate(nodes) %>%
#   mutate(tf=names(V(igraph.net_marker))) %>%
#   mutate(degree=igraph::degree(igraph.net_marker)) %>%
#   mutate(eigen_centrality=eigen_centrality(igraph.net_marker)$vector) %>%
#   activate(edges) %>%
#   mutate(sign=ifelse(E(igraph.net_marker)$weight>0,"Positive","Negative"))

# celltypes_markers <- rep(NA,length(names(V(igraph.net_marker))))
# for (i in 1:length(celltypes_markers)){
#   if (sum(markers_TF_undupl$TF==names(V(igraph.net_marker))[i])>0){
#     celltypes_markers[i] <- markers_TF_undupl$celltype[markers_TF_undupl$TF==names(V(igraph.net_marker))[i]]
#   }
# }
# celltypes_markers <- factor(as.character(celltypes_markers),levels=as.character(1:5))

# igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(celltype=celltypes_markers)

# if(!file.exists(file.path(args$outdir,sprintf("Score%s_Coef%s/global_network_allTFs_celltype.pdf",args$min_chip_score,args$min_coef))) | args$force){
#   # Plot
#   set.seed(42)
#   p <- ggraph(igraph.tbl_marker, 'stress') +
#     geom_edge_link(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
#                    arrow=arrow(length=unit(1.5,'mm')), end_cap=circle(1,'mm')) +
#     geom_node_point(aes(fill=celltype), size=3, stroke=0, shape=21, alpha=0.75) +
#     geom_node_text(aes(label=name), size=1) +
#     scale_edge_colour_manual(values=c("Positive"="darkred", "Negative"="darkblue")) +
#     scale_fill_manual(values=opts$celltype.colors,na.value="white") +
#     theme_graph(base_family = 'Helvetica')
  
#   pdf(file.path(args$outdir,sprintf("Score%s_Coef%s/global_network_allTFs_celltype.pdf",args$min_chip_score,args$min_coef)), width = 5.5, height = 5.75)
#   print(p)
#   dev.off()
# } else {
#   print("Skipped actual plotting")
# }

# print("Plot global network: celltype (marker TFs)")

# # Create node and edge data.frames
# node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
# edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$tf%in%TFs_marker) & (GRN_coef.dt$gene%in%TFs_marker)][,c("tf","gene","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]

# print("Summary edges:")
# n_pos_edges <- sum(edge_list.dt_marker$weight>0)
# n_neg_edges <- sum(edge_list.dt_marker$weight<0)
# print(paste0("positive: ",n_pos_edges))
# print(paste0("negative: ",n_neg_edges))

# # Create igraph object
# igraph.net_marker <- graph_from_data_frame(d = edge_list.dt_marker)

# # Create tbl_graph object for ggraph
# igraph.tbl_marker <- as_tbl_graph(igraph.net_marker) %>%
#   activate(nodes) %>%
#   mutate(tf=names(V(igraph.net_marker))) %>%
#   mutate(degree=igraph::degree(igraph.net_marker)) %>%
#   mutate(eigen_centrality=eigen_centrality(igraph.net_marker)$vector) %>%
#   activate(edges) %>%
#   mutate(sign=ifelse(E(igraph.net_marker)$weight>0,
#                      paste0("Positive (n=",n_pos_edges,")"),paste0("Negative (n=",n_neg_edges,")")))
# # igraph.tbl

# celltypes_markers <- rep(NA,length(names(V(igraph.net_marker))))
# for (i in 1:length(celltypes_markers)){
#   if (sum(markers_TF_undupl$TF==names(V(igraph.net_marker))[i])>0){
#     celltypes_markers[i] <- markers_TF_undupl$celltype[markers_TF_undupl$TF==names(V(igraph.net_marker))[i]]
#   }
# }
# celltypes_markers <- factor(as.character(celltypes_markers),levels=as.character(1:5))

# igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(celltype=celltypes_markers)

# set.seed(42)
# p <- ggraph(igraph.tbl_marker, 'stress') +
#   geom_edge_link(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
#                  arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(2,'mm'), start_cap=circle(2,'mm')) +
#   geom_node_point(aes(fill=celltype), size=5, stroke=0, shape=21, alpha=0.75) +
#   geom_node_text(aes(label=name), size=1.5) +
#   scale_edge_colour_manual(values=c("darkblue","darkred")) +
#   scale_fill_manual(values=opts$celltype.colors,na.value="white") +
#   theme_graph(base_family = 'Helvetica')

# pdf(file.path(args$outdir,sprintf("Score%s_Coef%s/global_network_markerTFs_celltype.pdf",args$min_chip_score,args$min_coef)), width = 5.5, height = 5.75)
# print(p)
# dev.off()

# #################################################
# ## Plot global network, coloured by centrality ##
# #################################################

# print("Plot global network: centrality")

# if(!file.exists(file.path(args$outdir,"Score%s_Coef%s/global_network_eigen_centrality.pdf")) | args$force){
#   # Plot
#   set.seed(42)
#   p <- ggraph(igraph.tbl, 'stress') +
#     geom_edge_link(edge_colour = "grey66", edge_alpha=0.50, edge_width=0.10, 
#                    arrow=arrow(length=unit(1.5,'mm')), end_cap=circle(4,'mm')) +
#     # geom_node_point(aes(fill=degree), size=8, shape=21) +
#     geom_node_point(aes(fill=eigen_centrality), size=4, shape=21, alpha=0.80) +
#     geom_node_text(aes(label=name), size=1) +
#     scale_fill_gradient(low = "white", high = "orange") +
#     theme_graph(base_family = 'Helvetica')
  
#   pdf(file.path(args$outdir,sprintf("Score%s_Coef%s/global_network_eigen_centrality_bis.pdf",args$min_chip_score,args$min_coef)), width = 4.5, height = 5)
#   print(p)
#   dev.off()
  
#   to.plot <- data.table(tf=names(V(igraph.net)), centrality = eigen_centrality(igraph.net)$vector) %>% 
#     .[sample(x=.N, size=round(.N/1.5))] %>%
#     setorder(-centrality) %>%
#     .[,tf:=factor(tf, levels=tf)] 
  
#   p <- ggplot(to.plot, aes_string(x="centrality", y="tf")) +
#     geom_point(aes(size=centrality)) +
#     geom_segment(aes(yend=tf), xend=0, size=0.20) +
#     scale_size_continuous(range = c(1.25,3.75)) +
#     labs(y="", x="Eigenvalue centrality") +
#     theme_classic() +
#     guides(size="none") +
#     theme(
#       axis.text.x = element_text(size=rel(1.0), color="black"),
#       axis.text.y = element_text(size=rel(0.75), color="black")
#     )
  
#   pdf(file.path(args$outdir,sprintf("Score%s_Coef%s/global_network_eigen_centrality.pdf",args$min_chip_score,args$min_coef)), width = 4, height = 4.5)
#   print(p)
#   dev.off()
# } else {
#   print("Skipped actual plotting")
# }


# ##################################
# ## Plot repressive interactions ##
# ##################################

# print("Plot repressive interactions")

# if(!file.exists(file.path(args$outdir,sprintf("Score%s_Coef%s/global_network_repressive_interactions.pdf",args$min_chip_score,args$min_coef))) | args$force){
#   # Plot
#   set.seed(42)
#   p <- ggraph(igraph.tbl, 'stress') +
#     geom_edge_link(aes(edge_colour = sign, edge_alpha = sign), edge_width=0.40, 
#                    arrow=arrow(length=unit(1.5,'mm')), end_cap=circle(4,'mm')) +
#     #geom_node_point(aes(fill=celltype), size=8, shape=21) +
#     geom_node_text(aes(label=name), size=1) +
#     scale_edge_alpha_manual(values=c("Positive"=0, "Negative"=0.80)) +
#     scale_edge_colour_manual(values=c("Positive"="darkred", "Negative"="darkblue")) +
#     scale_fill_manual(values=opts$celltype.colors) +
#     theme_graph(base_family = 'Helvetica')
  
#   pdf(file.path(args$outdir,sprintf("Score%s_Coef%s/global_network_repressive_interactions.pdf",args$min_chip_score,args$min_coef)), width = 5, height = 5.5)
#   print(p)
#   dev.off()
# } else {
#   print("Skipped actual plotting")
# }

# ##################################
# ## Plot activatory interactions ##
# ##################################

# print("Plot activatory interactions")

# if(!file.exists(file.path(args$outdir,sprintf("Score%s_Coef%s/global_network_activatory_interactions.pdf",args$min_chip_score,args$min_coef))) | args$force){
#   # Plot
#   set.seed(42)
#   p <- ggraph(igraph.tbl, 'stress') +
#     geom_edge_link(aes(edge_colour = sign, edge_alpha = sign), edge_width=0.20, 
#                    arrow=arrow(length=unit(1.5,'mm')), end_cap=circle(4,'mm')) +
#     #geom_node_point(aes(fill=celltype), size=8, shape=21) +
#     geom_node_text(aes(label=name), size=1) +
#     scale_edge_alpha_manual(values=c("Positive"=0.80, "Negative"=0)) +
#     scale_edge_colour_manual(values=c("Positive"="darkred", "Negative"="darkblue")) +
#     scale_fill_manual(values=opts$celltype.colors) +
#     theme_graph(base_family = 'Helvetica')
  
#   pdf(file.path(args$outdir,sprintf("Score%s_Coef%s/global_network_activatory_interactions.pdf",args$min_chip_score,args$min_coef)), width = 5.5, height = 5.5)
#   print(p)
#   dev.off()
# } else {
#   print("Skipped actual plotting")
# }

# ######################################################################################
# ## Highlight a TF with pleiotropy effects (both positive and negative interactions) ##
# ######################################################################################

# print("TF of interest")

# TFs_to_plot <- read.table(args$TFs_to_plot,header=F,sep="\n")[[1]]
# for (i in 1:length(TFs_to_plot)){
#   if (i == 1){
#     if(!dir.exists(file.path(args$outdir,sprintf("Score%s_Coef%s/TF_of_interest",args$min_chip_score,args$min_coef)))){
#       dir.create(file.path(args$outdir,sprintf("Score%s_Coef%s/TF_of_interest",args$min_chip_score,args$min_coef)))
#     } 
#   }
#   TF_of_interest <- TFs_to_plot[i]
#   print(TF_of_interest)
  
#   idx <- which(names(V(igraph.net))==TF_of_interest)
#   if(length(idx)==0){
#     print("not found!")
#     next()
#   }
#   test <- igraph.tbl %>% activate(edges) %>% filter(from==idx)
  
#   tmp <- names(V(igraph.net))[test %>% activate(edges) %>% as_tibble() %>% .$to]
#   test <- test %>% activate(nodes) %>% filter(name%in%c(TF_of_interest,tmp))
  
#   set.seed(42)
#   p <- ggraph(test, 'stress') +
#     geom_edge_link(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.35, 
#                    arrow=arrow(length=unit(1.5,'mm')), end_cap=circle(1,'mm'), 
#                    start_cap=circle(3,"mm")) +
#     # geom_node_point(aes(fill=celltype), size=11, shape=21, stroke=0, alpha=0.75) +
#     geom_node_text(aes(label=name), size=1.5) +
#     scale_edge_colour_manual(values=c("Positive"="darkred", "Negative"="darkblue")) +
#     # scale_fill_manual(values=opts$celltype.colors) +
#     theme_graph(base_family = 'Helvetica') + 
#     theme(legend.title = element_text(size=4), legend.text=element_text(size=4))
  
#   pdf(file.path(args$outdir,paste0("Score",args$min_chip_score,"_Coef",args$min_coef,"/TF_of_interest/",TF_of_interest,"_network.pdf")), width = 4, height = 4)
#   print(p)
#   dev.off()
# }

# ######################################################################################
# ## Highlight a TF with pleiotropy effects (both positive and negative interactions) ##
# ######################################################################################

# donators <- list(c("HIVEP3","ARID5B","ETV1","ESRRB"))
# for (j in 1:length(donators)){
#   for (k in 1:2){
#     if (k == 1){
#       interactions <- "activatory"
#     } else {
#       interactions <- "repressive"
#     }
    
#     print(paste0("Plot global network for cluster ",j))
    
#     # Create node and edge data.frames
#     TFs_marker_j <- TFs_marker[TFs_marker%in%donators[[j]]]
#     node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
#     edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$tf%in%TFs_marker) & (GRN_coef.dt$gene%in%TFs_marker)][,c("tf","gene","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]
#     print(TFs_marker_j)
    
#     # Create igraph object
#     igraph.net_marker <- graph_from_data_frame(d = edge_list.dt_marker)
    
#     # Create tbl_graph object for ggraph
#     igraph.tbl_marker <- as_tbl_graph(igraph.net_marker) %>%
#       activate(nodes) %>%
#       mutate(tf=names(V(igraph.net_marker))) %>%
#       mutate(degree=igraph::degree(igraph.net_marker)) %>%
#       mutate(eigen_centrality=eigen_centrality(igraph.net_marker)$vector) %>%
#       activate(edges) %>%
#       mutate(sign=ifelse(E(igraph.net_marker)$weight>0,"Positive","Negative"))
#     # print(length(E(igraph.net_marker)$weight))
#     mask_pos <- (E(igraph.tbl_marker)$weight>0) & !grepl("^Not",E(igraph.tbl_marker)$sign)
#     mask_not <- (!(edge_list.dt_marker$from%in%TFs_marker_j))
#     mask_neg <- (E(igraph.tbl_marker)$weight<0) & !grepl("^Not",E(igraph.tbl_marker)$sign)
#     # print(sum(mask_pos & mask_not))
#     # print(sum(mask_neg & mask_not))
#     E(igraph.tbl_marker)$sign[mask_not] <- paste0("Not in cluster (n=",sum(mask_not),")")
#     E(igraph.tbl_marker)$sign[mask_pos & !mask_not] <- paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")")
#     E(igraph.tbl_marker)$sign[mask_neg & !mask_not] <- paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")")
#     E(igraph.tbl_marker)$sign <- factor(E(igraph.tbl_marker)$sign,levels=c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
#                                                                            paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
#                                                                            paste0("Not in cluster (n=",sum(mask_not),")")))
    
#     labs_plot <- c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
#                    paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
#                    paste0("Not in cluster (n=",sum(mask_not),")"))
    
#     print(labs_plot)
    
#     celltypes_markers <- rep(NA,length(V(igraph.tbl_marker)$tf))
#     for (i in 1:length(celltypes_markers)){
#       if (sum(markers_TF_undupl$TF==V(igraph.tbl_marker)$tf[i])>0){
#         celltypes_markers[i] <- markers_TF_undupl$celltype[markers_TF_undupl$TF==V(igraph.tbl_marker)$tf[i]]
#       }
#     }
    
#     if (k == 1){
#       celltypes_markers[!((V(igraph.tbl_marker)$tf%in%TFs_marker_j) | (V(igraph.tbl_marker)$tf%in%GRN_coef.dt$gene[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$beta>0)]))] <- NA
#     } else {
#       celltypes_markers[!((V(igraph.tbl_marker)$tf%in%TFs_marker_j) | (V(igraph.tbl_marker)$tf%in%GRN_coef.dt$gene[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$beta<0)]))] <- NA
#     }
    
#     igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(celltype=celltypes_markers)
    
#     if (k == 1){
#       col_arrows <- c("darkred","gray","gray")
#     } else{
#       col_arrows <- c("gray","darkblue","gray")
#     }
    
#     celltypes_present <- sort(as.numeric(unique(celltypes_markers)))
#     celltypes_present <- celltypes_present[!is.na(celltypes_present)]
#     print(celltypes_present)
    
#     set.seed(42)
#     p <- ggraph(igraph.tbl_marker, 'stress') +
#       geom_edge_link(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
#                      arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(2,'mm'), start_cap=circle(2,'mm')) +
#       geom_node_point(aes(fill=celltype), size=5, stroke=0, shape=21, alpha=0.75) +
#       geom_node_text(aes(label=name), size=1.5,) +
#       scale_edge_colour_manual(values=col_arrows,
#                                breaks=waiver(),
#                                labels=labs_plot) +
#       scale_fill_manual(values=opts$celltype.colors[celltypes_present],na.value="gray") +
#       theme_graph(base_family = 'Helvetica')
    
#     if (k == 1){
#       pdf(file.path(args$outdir,sprintf("Score%s_Coef%s/global_network_markerTFs_tailored%s_%s.pdf",args$min_chip_score,args$min_coef,j,interactions)), 
#           width = 5.5, height = 5.75)
#     } else {
#       pdf(file.path(args$outdir,sprintf("Score%s_Coef%s/global_network_markerTFs_tailored%s_%s.pdf",args$min_chip_score,args$min_coef,j,interactions)), 
#           width = 5.5, height = 5.75)
#     }
#     print(p)
#     dev.off()
#   }
# }

# # ################################
# # ## Plot network per cell type ##
# # ################################
# # 
# # celltypes.to.plot <- c("Spinal_cord","NMP","Somitic_mesoderm")
# # 
# # for (i in celltypes.to.plot) {
# #   
# #   # Define node color based on TF expr
# #   expr.values <- rna_tf_pseudobulk_scaled.mtx[,i]
# #   expr.values[expr.values<=0.1] <- 0.1; expr.values[expr.values>=0.9] <- 0.9
# #   igraph.tbl <- igraph.tbl %>% activate(nodes) %>% mutate(expr=expr.values[names(V(igraph.net))])
# #   
# #   # (TO-DO) Gray out edges of genes that are not expressed
# #   
# #   set.seed(42)
# #   p <- ggraph(igraph.tbl, 'stress') +
# #     geom_edge_link(edge_colour="gray70", edge_alpha=0.50, edge_width=0.10, arrow=arrow(length=unit(1.5,'mm')), end_cap=circle(4,'mm')) +
# #     geom_node_point(aes(color=expr), size=8, stroke=0, alpha=0.75) +
# #     geom_node_text(aes(label=name), size=3.75) +
# #     # scale_color_gradient(low = "gray80", high = "#008B00") +
# #     viridis::scale_color_viridis(option="viridis") +
# #     theme_graph(base_family = 'Helvetica')+ 
# #     theme(
# #       legend.position = "right"
# #     )
# #   
# #   pdf(file.path(args$outdir,sprintf("network_coloured_by_%s_expr.pdf",i)), width = 5, height = 5.5)
# #   print(p)
# #   dev.off()
# #   
# # }

######################
## Completion token ##
######################

file.create(file.path(args$outdir,sprintf("Score%s_Coef%s/completed.txt",args$min_chip_score,args$min_coef)))

##################
# Print warnings #
##################

warnings()