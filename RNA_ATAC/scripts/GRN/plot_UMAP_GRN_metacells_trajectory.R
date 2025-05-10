##############################################
##                                          ##
##  Plot_UMAP_GRN_metacells_trajectories.R  ##
##                                          ##
##############################################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA_ATAC/scripts/Settings.R")

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
p$add_argument('--rna_pseudobulk_cluster.sce',       type="character",                help='SingleCellExperiment (pseudobulk) based on cluster')
p$add_argument('--rna_pseudobulk_sample.sce',       type="character",                help='SingleCellExperiment (pseudobulk) based on sample')
p$add_argument('--grn_coef',       type="character",                help='GRN coefficients')
p$add_argument('--min_coef',  type="double",            default=0.25,      help='Minimal coefficient for linear regression')
p$add_argument('--max_distance',  type="integer",            default=1e5,      help='Maximum distance for a linkage between a peak and a gene')
p$add_argument('--max_pval',  type="double",            default=0.25,      help='Maximal pvalue for linear regression')
p$add_argument('--min_chip_score',  type="double",            default=0.25,      help='Minimal chip score')
p$add_argument('--TFs_to_plot',       type="character",                help='File containing TF of interest to plot')
p$add_argument('--force',       type="logical",  default=F         ,      help='Should replotting be forced?')
p$add_argument('--seed',       type="integer",  default=42         ,      help='Random seed')
p$add_argument('--tf2gene_virtual_chip',       type="character",                help='Links between tfs and genes based on virtual chipseq')
p$add_argument('--outdir',       type="character",                help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

args$trajectory <- sprintf("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA_ATAC/trajectories/%s/%s_trajectory.txt.gz",
                           args$trajectory_name,args$trajectory_name)

if(!dir.exists(file.path(args$outdir,sprintf("Score%s",args$min_chip_score)))){
  dir.create(args$outdir, showWarnings = F)
  dir.create(file.path(args$outdir,sprintf("Score%s",args$min_chip_score)))
}

# Options
opts$celltype_trajectory_dic <- list(
  "N2P" = c(1,2,3,4,5,6),
  "superclusterA" = c(1,2,3),
  "superclusterB" = c(4,5),
  "superclusterC" = 6
)

stopifnot(args$trajectory_name%in%names(opts$celltype_trajectory_dic))
celltypes.to.plot <- opts$celltype_trajectory_dic[[args$trajectory_name]]

##############
## Load TFs ##
##############

# all_TFs
all_TFs <- read.table("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/MOUSE_TF_aertslab.txt",sep="\n",header=F,quote="")$V1
print("All TFs")
print(head(all_TFs))

# marker TFs
markers_TF <- NULL
for (i in 1:6){
  marker_tf_i <- read.table(sprintf("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA_ATAC/pseudobulk/cluster/RNA/marker_TFs/Marker_TF_cluster%s.txt",i),header=T,sep="\t")
  markers_TF <- rbind(markers_TF,cbind(marker_tf_i$Gene,rep(i,length(marker_tf_i$Gene))))
}
markers_TF <- data.frame(markers_TF)
colnames(markers_TF) <- c("TF","celltype")
markers_TF_undupl <- markers_TF[!duplicated(markers_TF$TF),]
print("Marker TFs")
print(nrow(markers_TF_undupl))
print(head(markers_TF_undupl))

##############################
## Load RNA expression data ##
##############################

# SingleCellExperiment at metacell resolution
rna_metacells.sce <- readRDS(args$rna_metacells.sce)

if(grepl("supercluster",args$trajectory_name)){
  metadata <- fread('/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA_ATAC/clustering/PeakMatrix/sample_metadata_all_after_clustering.txt.gz')
  metadata <- metadata[metadata$sample %in% colnames(rna_metacells.sce),]
  metadata <- metadata[match(colnames(rna_metacells.sce), metadata$sample),]

  if (grepl("A",args$trajectory_name)){
    clusters <- c(1, 2, 3)
  } else if (grepl("B",args$trajectory_name)) {
    clusters <- c(4, 5)
  } else if (grepl("C",args$trajectory_name)) {
    clusters <- 6
  }

  selected_cells <- metadata$sample[metadata$cluster %in% clusters]
  rna_metacells.sce <- rna_metacells.sce[, colnames(rna_metacells.sce) %in% selected_cells]
}
print("rna_metacells.sce")
rna_metacells.sce

# SingleCellExperiment at pseudobulk resolution
rna_pseudobulk.sce <- readRDS(args$rna_pseudobulk_cluster.sce)#[,celltypes.to.plot]
print("rna_pseudobulk.sce")
rna_pseudobulk.sce

# SingleCellExperiment at pseudobulk resolution
rna_pseudobulk_sample.sce <- readRDS(args$rna_pseudobulk_sample.sce)#[,celltypes.to.plot]
print("rna_pseudobulk_sample.sce")
rna_pseudobulk_sample.sce

##################################
## Load global GRN coefficients ##
##################################

GRN_coef.dt <- fread(args$grn_coef) %>% 
  # .[,gene:=toupper(gene)] %>% # .[gene%in%unique(marker_TFs_filt.dt$gene) & tf%in%unique(marker_TFs_filt.dt$gene)] %>%
  .[pvalue<args$max_pval & abs(beta)>=args$min_coef]
print(head(GRN_coef.dt))
print(dim(GRN_coef.dt))

##########################
## Filter TFs and genes ##
##########################

TFs_marker <- markers_TF_undupl$TF

length(TFs_marker)

TFs <- TFs_marker
print("TFs")
head(TFs)
print("GRN_coef.dt")
head(GRN_coef.dt)
GRN_coef.dt <- GRN_coef.dt[tf%in%TFs & gene%in%TFs]
print(dim(GRN_coef.dt))

# more TFs in TFs object than in sce data (Check how this is possible!)
print("Length TFs and amounts in sce data")
length(TFs)
sum(TFs%in%rownames(logcounts(rna_pseudobulk.sce)))
sum(TFs%in%rownames(logcounts(rna_metacells.sce)))
TFs <- TFs[(TFs%in%rownames(logcounts(rna_pseudobulk.sce)))]

# Fetch RNA expression matrices
# head(logcounts(rna_pseudobulk.sce))
rna_tf_pseudobulk.mtx <- logcounts(rna_pseudobulk.sce)[TFs,]
# head(logcounts(rna_metacells.sce))
rna_tf_metacells.mtx <- logcounts(rna_metacells.sce)[TFs,]


################################
## Load in silico Chip scores ##
################################

tf2gene_chip.dt <- fread(args$tf2gene_virtual_chip) %>%
  .[chip_score>=args$min_chip_score & dist<=args$max_distance] %>% 
  .[,c("tf","gene","chip_score")] %>% unique

print("In silico ChipSeq")
print(head(tf2gene_chip.dt))
print(dim(tf2gene_chip.dt))
# print(sum(tf2gene_chip.dt$tf%in%TFs))
# print(sum(tf2gene_chip.dt$gene%in%TFs))

################################
## Determine UMAP coordinates ##
################################

print("Create Beta matrix")

args$beta_mtx_cte <- 0.2

if (!file.exists(file.path(args$outdir,sprintf("Score%s/betas_mtx.txt",args$min_chip_score))) | args$force){
  betas_mtx <- matrix(rep(0,length(TFs)*length(TFs)),ncol=length(TFs))
  colnames(betas_mtx) <- sort(TFs)
  rownames(betas_mtx) <- sort(TFs)
  print(ncol(betas_mtx))
  for (i in 1:ncol(betas_mtx)){
    if((i%%1000)==0){print(i)}
    if ((sum(GRN_coef.dt$tf==colnames(betas_mtx)[i])>0) & (sum(tf2gene_chip.dt$tf==colnames(betas_mtx)[i])>0)){
      GRN_coef.dt_TF_i <- GRN_coef.dt[GRN_coef.dt$tf==colnames(betas_mtx)[i],]
      GRN_coef.dt_TF_i <- GRN_coef.dt_TF_i[GRN_coef.dt_TF_i$gene%in%TFs,]
      GRN_coef.dt_TF_i <- GRN_coef.dt_TF_i[sort(GRN_coef.dt_TF_i$gene,index.return=T)$ix,]
      betas_mtx[rownames(betas_mtx)%in%GRN_coef.dt_TF_i$gene,i] <- GRN_coef.dt_TF_i$beta
      
      tf2gene_chip.dt_TF_i <- tf2gene_chip.dt[tf2gene_chip.dt$tf==colnames(betas_mtx)[i],]
      tf2gene_chip.dt_TF_i <- tf2gene_chip.dt_TF_i[tf2gene_chip.dt_TF_i$gene%in%TFs,]
      tf2gene_chip.dt_TF_i <- tf2gene_chip.dt_TF_i[sort(tf2gene_chip.dt_TF_i$chip_score,index.return=T,decreasing = T)$ix,]
      tf2gene_chip.dt_TF_i <- tf2gene_chip.dt_TF_i[!duplicated(tf2gene_chip.dt_TF_i$gene),]
      tf2gene_chip.dt_TF_i <- tf2gene_chip.dt_TF_i[sort(tf2gene_chip.dt_TF_i$gene,index.return=T)$ix,]
      # betas_mtx[rownames(betas_mtx)%in%tf2gene_chip.dt_TF_i$gene,i] <-
      # betas_mtx[rownames(betas_mtx)%in%tf2gene_chip.dt_TF_i$gene,i] * sqrt(abs(tf2gene_chip.dt_TF_i$chip_score))
      betas_mtx[rownames(betas_mtx)%in%tf2gene_chip.dt_TF_i$gene,i] <-
       (betas_mtx[rownames(betas_mtx)%in%tf2gene_chip.dt_TF_i$gene,i] * sqrt(abs(tf2gene_chip.dt_TF_i$chip_score)) + args$beta_mtx_cte)
    }
  }
  print("Betas table output")
  print(dim(betas_mtx))
  print(head(betas_mtx[,1:length(celltypes.to.plot)]))

  write.table(betas_mtx,file=file.path(args$outdir,sprintf("Score%s/betas_mtx.txt",args$min_chip_score)),
              col.names = T, row.names = T, sep="\t", quote=F)
} else {
  print("Skipped creating Betas, loaded from memory. Use 'Force=T' to force a rerun.")
  betas_mtx <- read.table(file=file.path(args$outdir,sprintf("Score%s/betas_mtx.txt",args$min_chip_score)),
                          header=T,sep="\t",row.names=1)
  print(dim(betas_mtx))
  print(head(betas_mtx[,1:celltypes.to.plot]))
}

# Run
set.seed(args$seed)
umap_embedding <- uwot::umap(betas_mtx, spread = 1, min_dist = 0.5, n_neighbors=10, pca = 15)
# umap_embedding <- uwot::umap(betas_mtx, pca = 20)

print(dim(umap_embedding))
print(head(umap_embedding))

##############################################
## Plot network per cluster on general plot ##
##############################################

for (j in celltypes.to.plot){
  for (k in 1:2){
    if (k == 1){
      interactions <- "activatory"
    } else {
      interactions <- "repressive"
    }
    
    print(paste0("Plot global network for cluster ",j))
    
    # Create node and edge data.frames
    TFs_marker_j <- TFs_marker[TFs_marker%in%markers_TF_undupl$TF[markers_TF_undupl$celltype==j]]
    node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
    edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$tf%in%TFs_marker) & (GRN_coef.dt$gene%in%TFs_marker)][,c("tf","gene","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]
    print("Node and edge statistics")
    print(dim(node_list.dt_marker))
    print(head(node_list.dt_marker))
    print(dim(edge_list.dt_marker))
    print(head(edge_list.dt_marker))
    
    # Create igraph object
    igraph.net_marker <- graph_from_data_frame(d = edge_list.dt_marker)
    
    # Create tbl_graph object for ggraph
    igraph.tbl_marker <- as_tbl_graph(igraph.net_marker) %>%
      activate(nodes) %>%
      mutate(tf=names(V(igraph.net_marker))) %>%
      mutate(degree=igraph::degree(igraph.net_marker)) %>%
      mutate(eigen_centrality=eigen_centrality(igraph.net_marker)$vector) %>%
      activate(edges) %>%
      mutate(sign=ifelse(E(igraph.net_marker)$weight>0,"Positive","Negative"))
    # print(length(E(igraph.net_marker)$weight))
    mask_pos <- (E(igraph.tbl_marker)$weight>0) & !grepl("^Not",E(igraph.tbl_marker)$sign)
    mask_not <- (!(edge_list.dt_marker$from%in%TFs_marker_j))
    mask_neg <- (E(igraph.tbl_marker)$weight<0) & !grepl("^Not",E(igraph.tbl_marker)$sign)
    # print(sum(mask_pos & mask_not))
    # print(sum(mask_neg & mask_not))
    E(igraph.tbl_marker)$sign[mask_not] <- paste0("Not in cluster (n=",sum(mask_not),")")
    E(igraph.tbl_marker)$sign[mask_pos & !mask_not] <- paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")")
    E(igraph.tbl_marker)$sign[mask_neg & !mask_not] <- paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")")
    E(igraph.tbl_marker)$sign <- factor(E(igraph.tbl_marker)$sign,levels=c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
                                                                           paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
                                                                           paste0("Not in cluster (n=",sum(mask_not),")")))
    
    labs_plot <- c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
                   paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
                   paste0("Not in cluster (n=",sum(mask_not),")"))
    
    print(labs_plot)
    
    celltypes_markers <- rep(NA,length(V(igraph.tbl_marker)$tf))
    for (i in 1:length(celltypes_markers)){
      if (sum(markers_TF_undupl$TF==V(igraph.tbl_marker)$tf[i])>0){
        celltypes_markers[i] <- markers_TF_undupl$celltype[markers_TF_undupl$TF==V(igraph.tbl_marker)$tf[i]]
      }
    }
    
    if (k == 1){
      celltypes_markers[!((celltypes_markers==j) | (V(igraph.tbl_marker)$tf%in%GRN_coef.dt$gene[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$beta>0)]))] <- NA
    } else {
      celltypes_markers[!((celltypes_markers==j) | (V(igraph.tbl_marker)$tf%in%GRN_coef.dt$gene[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$beta<0)]))] <- NA
    }
    
    igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(celltype=celltypes_markers)
    
    if (k == 1){
      col_arrows <- c("darkred","gray","gray")
    } else{
      col_arrows <- c("gray","darkblue","gray")
    }
    
    celltypes_present <- sort(as.numeric(unique(celltypes_markers)))
    celltypes_present <- celltypes_present[!is.na(celltypes_present)]
    print(celltypes_present)
    
    set.seed(42)
    p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
                y=umap_embedding[names(V(igraph.net_marker)),1]) +
      geom_edge_diagonal(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
                     arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
      geom_node_point(aes(fill=celltype), size=7.5, stroke=0, shape=21, alpha=0.75) +
      geom_node_text(aes(label=name), size=2.25,) +
      scale_edge_colour_manual(values=col_arrows,
                               breaks=waiver(),
                               labels=labs_plot) +
      scale_fill_manual(values=opts$celltype.colors[celltypes_present],na.value="gray") +
      theme_graph(base_family = 'Helvetica') + theme(legend.position = "none")
    
    if (k == 1){
      pdf(file.path(args$outdir,sprintf("Score%s/global_network_markerTFs_cluster%s_%s.pdf",args$min_chip_score,j,interactions)), 
          width = 6, height = 5.5)
    } else {
      pdf(file.path(args$outdir,sprintf("Score%s/global_network_markerTFs_cluster%s_%s.pdf",args$min_chip_score,j,interactions)), 
          width = 6, height = 5.5)
    }
    print(p)
    dev.off()
  }
}

if(grepl("N2P",args$trajectory_name)){
  ###################################################
  ## Plot network per supercluster on general plot ##
  ###################################################

  # define super‐clusters mapping
  superclusters <- list(
    A = 1:3,
    B = 4:5,
    C = 6
  )

  for (sc in names(superclusters)) {
    clusts <- superclusters[[sc]]
    
    for (k in 1:2) {
      if (k == 1) {
        interactions <- "activatory"
      } else {
        interactions <- "repressive"
      }
      
      print(paste0("Plot global network for super‐cluster ",sc))
      
      # Create node and edge data.frames
      TFs_marker_j <- TFs_marker[TFs_marker %in% markers_TF_undupl$TF[markers_TF_undupl$celltype %in% clusts]]
      node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
      edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$tf%in%TFs_marker) & (GRN_coef.dt$gene%in%TFs_marker)][,c("tf","gene","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]
      print("Node and edge statistics")
      print(dim(node_list.dt_marker))
      print(head(node_list.dt_marker))
      print(dim(edge_list.dt_marker))
      print(head(edge_list.dt_marker))
      
      # Create igraph object
      igraph.net_marker <- graph_from_data_frame(d = edge_list.dt_marker)
      
      # Create tbl_graph object for ggraph
      igraph.tbl_marker <- as_tbl_graph(igraph.net_marker) %>%
        activate(nodes) %>%
        mutate(tf=names(V(igraph.net_marker))) %>%
        mutate(degree=igraph::degree(igraph.net_marker)) %>%
        mutate(eigen_centrality=eigen_centrality(igraph.net_marker)$vector) %>%
        activate(edges) %>%
        mutate(sign=ifelse(E(igraph.net_marker)$weight>0,"Positive","Negative"))
      # print(length(E(igraph.net_marker)$weight))
      mask_pos <- (E(igraph.tbl_marker)$weight>0) & !grepl("^Not",E(igraph.tbl_marker)$sign)
      mask_not <- (!(edge_list.dt_marker$from%in%TFs_marker_j))
      mask_neg <- (E(igraph.tbl_marker)$weight<0) & !grepl("^Not",E(igraph.tbl_marker)$sign)
      # print(sum(mask_pos & mask_not))
      # print(sum(mask_neg & mask_not))
      E(igraph.tbl_marker)$sign[mask_not] <- paste0("Not in cluster (n=",sum(mask_not),")")
      E(igraph.tbl_marker)$sign[mask_pos & !mask_not] <- paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")")
      E(igraph.tbl_marker)$sign[mask_neg & !mask_not] <- paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")")
      E(igraph.tbl_marker)$sign <- factor(E(igraph.tbl_marker)$sign,levels=c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
                                                                            paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
                                                                            paste0("Not in cluster (n=",sum(mask_not),")")))
      
      labs_plot <- c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
                    paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
                    paste0("Not in cluster (n=",sum(mask_not),")"))
      
      print(labs_plot)
      
      celltypes_markers <- rep(NA,length(V(igraph.tbl_marker)$tf))
      for (i in 1:length(celltypes_markers)){
        if (sum(markers_TF_undupl$TF==V(igraph.tbl_marker)$tf[i])>0){
          celltypes_markers[i] <- markers_TF_undupl$celltype[markers_TF_undupl$TF==V(igraph.tbl_marker)$tf[i]]
        }
      }
      
      if (k == 1){
        celltypes_markers[!((celltypes_markers==j) | (V(igraph.tbl_marker)$tf%in%GRN_coef.dt$gene[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$beta>0)]))] <- NA
      } else {
        celltypes_markers[!((celltypes_markers==j) | (V(igraph.tbl_marker)$tf%in%GRN_coef.dt$gene[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$beta<0)]))] <- NA
      }
      
      igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(celltype=celltypes_markers)
      
      if (k == 1){
        col_arrows <- c("darkred","gray","gray")
      } else{
        col_arrows <- c("gray","darkblue","gray")
      }
      
      celltypes_present <- sort(as.numeric(unique(celltypes_markers)))
      celltypes_present <- celltypes_present[!is.na(celltypes_present)]
      print(celltypes_present)

      set.seed(42)
      p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
                  y=umap_embedding[names(V(igraph.net_marker)),1]) +
        geom_edge_diagonal(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
                      arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
        geom_node_point(aes(fill=celltype), size=7.5, stroke=0, shape=21, alpha=0.75) +
        geom_node_text(aes(label=name), size=2.25,) +
        scale_edge_colour_manual(values=col_arrows,
                                breaks=waiver(),
                                labels=labs_plot) +
        scale_fill_manual(values=opts$celltype.colors[celltypes_present],na.value="gray") +
        theme_graph(base_family = 'Helvetica') + theme(legend.position = "none")

      # save PDF
      pdf(file.path(args$outdir, sprintf("Score%s/global_network_markerTFs_supercluster_%s_%s.pdf",args$min_chip_score, sc, interactions)),
          width  = 6, height = 5.5)
      print(p)
      dev.off()
    }
  }
}

###############################################
## Plot outgoing per cluster on general plot ##
###############################################

for (j in celltypes.to.plot){
  for (k in 1:2) {
    if (k == 1){
      interactions <- "activatory"
    } else {
      interactions <- "repressive"
    }
    
    print(paste0("Plot global network for cluster ",j))
    
    # Create node and edge data.frames
    TFs_marker_j <- TFs_marker[TFs_marker%in%markers_TF_undupl$TF[markers_TF_undupl$celltype==j]]
    TFs_marker_jand1 <- TFs_marker[TFs_marker%in%markers_TF_undupl$TF[markers_TF_undupl$celltype==j+1]]
    
    node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
    edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$tf%in%TFs_marker) & (GRN_coef.dt$gene%in%TFs_marker)][,c("tf","gene","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]
    print(dim(edge_list.dt_marker))
    print(head(edge_list.dt_marker))
    
    # Create igraph object
    igraph.net_marker <- graph_from_data_frame(d = edge_list.dt_marker)
    
    for (m in 1:2){
      if (m == 1){
        # Create tbl_graph object for ggraph
        igraph.tbl_marker <- as_tbl_graph(igraph.net_marker) %>%
          activate(nodes) %>%
          mutate(tf=names(V(igraph.net_marker))) %>%
          mutate(degree=igraph::degree(igraph.net_marker)) %>%
          mutate(eigen_centrality=eigen_centrality(igraph.net_marker)$vector) %>%
          activate(edges) %>%
          mutate(sign=ifelse(E(igraph.net_marker)$weight>0,"Positive","Negative"))
        # print(length(E(igraph.net_marker)$weight))
        mask_pos <- (E(igraph.tbl_marker)$weight>0) & !grepl("^Not",E(igraph.tbl_marker)$sign)
        mask_not <- (!((edge_list.dt_marker$from%in%TFs_marker_j) & (edge_list.dt_marker$to%in%TFs_marker_j)))
        mask_neg <- (E(igraph.tbl_marker)$weight<0) & !grepl("^Not",E(igraph.tbl_marker)$sign)
        # print(sum(mask_pos & mask_not))
        # print(sum(mask_neg & mask_not))
        E(igraph.tbl_marker)$sign[mask_not] <- paste0("Not in cluster (n=",sum(mask_not),")")
        E(igraph.tbl_marker)$sign[mask_pos & !mask_not] <- paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")")
        E(igraph.tbl_marker)$sign[mask_neg & !mask_not] <- paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")")
        E(igraph.tbl_marker)$sign <- factor(E(igraph.tbl_marker)$sign,levels=c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
                                                                              paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
                                                                              paste0("Not in cluster (n=",sum(mask_not),")")))
        
        labs_plot <- c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
                      paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
                      paste0("Not in cluster (n=",sum(mask_not),")"))
        
        print(labs_plot)
        
        celltypes_markers <- rep(NA,length(V(igraph.tbl_marker)$tf))
        for (i in 1:length(celltypes_markers)){
          if (sum(markers_TF_undupl$TF==V(igraph.tbl_marker)$tf[i])>0){
            celltypes_markers[i] <- markers_TF_undupl$celltype[markers_TF_undupl$TF==V(igraph.tbl_marker)$tf[i]]
          }
        }
        
        if (k == 1){
          celltypes_markers[!((V(igraph.tbl_marker)$tf%in%GRN_coef.dt$tf[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$gene%in%TFs_marker_j) & (GRN_coef.dt$beta>0)]) |
                                (V(igraph.tbl_marker)$tf%in%GRN_coef.dt$gene[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$gene%in%TFs_marker_j) & (GRN_coef.dt$beta>0)]))] <- NA
        } else {
          celltypes_markers[!((V(igraph.tbl_marker)$tf%in%GRN_coef.dt$tf[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$gene%in%TFs_marker_j) & (GRN_coef.dt$beta<0)]) |
                                (V(igraph.tbl_marker)$tf%in%GRN_coef.dt$gene[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$gene%in%TFs_marker_j) & (GRN_coef.dt$beta<0)]))] <- NA
        }
        
        igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(celltype=celltypes_markers)
        
        if (k == 1){
          col_arrows <- c("darkred","gray","gray")
        } else{
          col_arrows <- c("gray","darkblue","gray")
        }
        
        celltypes_present <- sort(as.numeric(unique(celltypes_markers)))
        celltypes_present <- celltypes_present[!is.na(celltypes_present)]
        print(celltypes_present)
        
        set.seed(42)
        p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
                    y=umap_embedding[names(V(igraph.net_marker)),1]) +
          geom_edge_diagonal(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
                            arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
          geom_node_point(aes(fill=celltype), size=7.5, stroke=0, shape=21, alpha=0.75) +
          geom_node_text(aes(label=name), size=2.25,) +
          scale_edge_colour_manual(values=col_arrows,
                                  breaks=waiver(),
                                  labels=labs_plot) +
          scale_fill_manual(values=opts$celltype.colors[celltypes_present],na.value="gray") +
          theme_graph(base_family = 'Helvetica') + theme(legend.position = "none")
        
        if (k == 1){
          pdf(file.path(args$outdir,sprintf("Score%s/global_network_markerTFs_cluster%sto%s_%s.pdf",args$min_chip_score,j,j,interactions)), 
              width = 6, height = 5.5)
        } else {
          pdf(file.path(args$outdir,sprintf("Score%s/global_network_markerTFs_cluster%sto%s_%s.pdf",args$min_chip_score,j,j,interactions)), 
              width = 6, height = 5.5)
        }
        print(p)
        dev.off()
      } else {
        if (j == celltypes.to.plot[length(celltypes.to.plot)]){
          next()
        }
        # Create tbl_graph object for ggraph
        igraph.tbl_marker <- as_tbl_graph(igraph.net_marker) %>%
          activate(nodes) %>%
          mutate(tf=names(V(igraph.net_marker))) %>%
          mutate(degree=igraph::degree(igraph.net_marker)) %>%
          mutate(eigen_centrality=eigen_centrality(igraph.net_marker)$vector) %>%
          activate(edges) %>%
          mutate(sign=ifelse(E(igraph.net_marker)$weight>0,"Positive","Negative"))
        # print(length(E(igraph.net_marker)$weight))
        mask_pos <- (E(igraph.tbl_marker)$weight>0) & !grepl("^Not",E(igraph.tbl_marker)$sign)
        mask_not <- (!((edge_list.dt_marker$from%in%TFs_marker_j) & (edge_list.dt_marker$to%in%TFs_marker_jand1)))
        mask_neg <- (E(igraph.tbl_marker)$weight<0) & !grepl("^Not",E(igraph.tbl_marker)$sign)
        # print(sum(mask_pos & mask_not))
        # print(sum(mask_neg & mask_not))
        E(igraph.tbl_marker)$sign[mask_not] <- paste0("Not in cluster (n=",sum(mask_not),")")
        E(igraph.tbl_marker)$sign[mask_pos & !mask_not] <- paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")")
        E(igraph.tbl_marker)$sign[mask_neg & !mask_not] <- paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")")
        E(igraph.tbl_marker)$sign <- factor(E(igraph.tbl_marker)$sign,levels=c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
                                                                              paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
                                                                              paste0("Not in cluster (n=",sum(mask_not),")")))
        
        labs_plot <- c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
                      paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
                      paste0("Not in cluster (n=",sum(mask_not),")"))
        
        print(labs_plot)
        
        celltypes_markers <- rep(NA,length(V(igraph.tbl_marker)$tf))
        for (i in 1:length(celltypes_markers)){
          if (sum(markers_TF_undupl$TF==V(igraph.tbl_marker)$tf[i])>0){
            celltypes_markers[i] <- markers_TF_undupl$celltype[markers_TF_undupl$TF==V(igraph.tbl_marker)$tf[i]]
          }
        }
        
        if (k == 1){
          celltypes_markers[!((V(igraph.tbl_marker)$tf%in%GRN_coef.dt$tf[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$gene%in%TFs_marker_jand1) & (GRN_coef.dt$beta>0)]) |
                                (V(igraph.tbl_marker)$tf%in%GRN_coef.dt$gene[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$gene%in%TFs_marker_jand1) & (GRN_coef.dt$beta>0)]))] <- NA
        } else {
          celltypes_markers[!((V(igraph.tbl_marker)$tf%in%GRN_coef.dt$tf[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$gene%in%TFs_marker_jand1) & (GRN_coef.dt$beta<0)]) |
                                (V(igraph.tbl_marker)$tf%in%GRN_coef.dt$gene[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$gene%in%TFs_marker_jand1) & (GRN_coef.dt$beta<0)]))] <- NA
        }
        
        igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(celltype=celltypes_markers)
        
        if (k == 1){
          col_arrows <- c("darkred","gray","gray")
        } else{
          col_arrows <- c("gray","darkblue","gray")
        }
        
        celltypes_present <- sort(as.numeric(unique(celltypes_markers)))
        celltypes_present <- celltypes_present[!is.na(celltypes_present)]
        print(celltypes_present)
        
        set.seed(42)
        p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
                    y=umap_embedding[names(V(igraph.net_marker)),1]) +
          geom_edge_diagonal(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
                            arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
          geom_node_point(aes(fill=celltype), size=7.5, stroke=0, shape=21, alpha=0.75) +
          geom_node_text(aes(label=name), size=2.25,) +
          scale_edge_colour_manual(values=col_arrows,
                                  breaks=waiver(),
                                  labels=labs_plot) +
          scale_fill_manual(values=opts$celltype.colors[celltypes_present],na.value="gray") +
          theme_graph(base_family = 'Helvetica') + theme(legend.position = "none")
        
        if (k == 1){
          pdf(file.path(args$outdir,sprintf("Score%s/global_network_markerTFs_cluster%sto%s_%s.pdf",args$min_chip_score,j,j+1,interactions)), 
              width = 6, height = 5.5)
        } else {
          pdf(file.path(args$outdir,sprintf("Score%s/global_network_markerTFs_cluster%sto%s_%s.pdf",args$min_chip_score,j,j+1,interactions)), 
              width = 6, height = 5.5)
        }
        print(p)
        dev.off()
      }
    }
  }
}

if(grepl("N2P",args$trajectory_name)){
  ####################################################
  ## Plot outgoing per supercluster on general plot ##
  ####################################################

  for (sc in names(superclusters)) {
    clusts <- superclusters[[sc]]
    
    for (k in 1:2) {
      if (k == 1) {
        interactions <- "activatory"
      } else {
        interactions <- "repressive"
      }
    
      print(paste0("Plot global network for super‐cluster ",sc))
  
      # Create node and edge data.frames
      TFs_marker_j <- TFs_marker[TFs_marker%in%markers_TF_undupl$TF[markers_TF_undupl$celltype %in% clusts]]
      TFs_marker_jand1 <- TFs_marker[TFs_marker%in%markers_TF_undupl$TF[markers_TF_undupl$celltype==6]]
    
      node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
      edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$tf%in%TFs_marker) & (GRN_coef.dt$gene%in%TFs_marker)][,c("tf","gene","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]
      print(dim(edge_list.dt_marker))
      print(head(edge_list.dt_marker))
    
      # Create igraph object
      igraph.net_marker <- graph_from_data_frame(d = edge_list.dt_marker)
      
      for (m in 1:2){
        if (m == 1){
          # Create tbl_graph object for ggraph
          igraph.tbl_marker <- as_tbl_graph(igraph.net_marker) %>%
            activate(nodes) %>%
            mutate(tf=names(V(igraph.net_marker))) %>%
            mutate(degree=igraph::degree(igraph.net_marker)) %>%
            mutate(eigen_centrality=eigen_centrality(igraph.net_marker)$vector) %>%
            activate(edges) %>%
            mutate(sign=ifelse(E(igraph.net_marker)$weight>0,"Positive","Negative"))
          # print(length(E(igraph.net_marker)$weight))
          mask_pos <- (E(igraph.tbl_marker)$weight>0) & !grepl("^Not",E(igraph.tbl_marker)$sign)
          mask_not <- (!((edge_list.dt_marker$from%in%TFs_marker_j) & (edge_list.dt_marker$to%in%TFs_marker_j)))
          mask_neg <- (E(igraph.tbl_marker)$weight<0) & !grepl("^Not",E(igraph.tbl_marker)$sign)
          # print(sum(mask_pos & mask_not))
          # print(sum(mask_neg & mask_not))
          E(igraph.tbl_marker)$sign[mask_not] <- paste0("Not in cluster (n=",sum(mask_not),")")
          E(igraph.tbl_marker)$sign[mask_pos & !mask_not] <- paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")")
          E(igraph.tbl_marker)$sign[mask_neg & !mask_not] <- paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")")
          E(igraph.tbl_marker)$sign <- factor(E(igraph.tbl_marker)$sign,levels=c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
                                                                                paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
                                                                                paste0("Not in cluster (n=",sum(mask_not),")")))
          
          labs_plot <- c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
                        paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
                        paste0("Not in cluster (n=",sum(mask_not),")"))
          
          print(labs_plot)
          
          celltypes_markers <- rep(NA,length(V(igraph.tbl_marker)$tf))
          for (i in 1:length(celltypes_markers)){
            if (sum(markers_TF_undupl$TF==V(igraph.tbl_marker)$tf[i])>0){
              celltypes_markers[i] <- markers_TF_undupl$celltype[markers_TF_undupl$TF==V(igraph.tbl_marker)$tf[i]]
            }
          }
          
          if (k == 1){
            celltypes_markers[!((V(igraph.tbl_marker)$tf%in%GRN_coef.dt$tf[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$gene%in%TFs_marker_j) & (GRN_coef.dt$beta>0)]) |
                                  (V(igraph.tbl_marker)$tf%in%GRN_coef.dt$gene[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$gene%in%TFs_marker_j) & (GRN_coef.dt$beta>0)]))] <- NA
          } else {
            celltypes_markers[!((V(igraph.tbl_marker)$tf%in%GRN_coef.dt$tf[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$gene%in%TFs_marker_j) & (GRN_coef.dt$beta<0)]) |
                                  (V(igraph.tbl_marker)$tf%in%GRN_coef.dt$gene[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$gene%in%TFs_marker_j) & (GRN_coef.dt$beta<0)]))] <- NA
          }
          
          igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(celltype=celltypes_markers)
          
          if (k == 1){
            col_arrows <- c("darkred","gray","gray")
          } else{
            col_arrows <- c("gray","darkblue","gray")
          }
          
          celltypes_present <- sort(as.numeric(unique(celltypes_markers)))
          celltypes_present <- celltypes_present[!is.na(celltypes_present)]
          print(celltypes_present)
          
          set.seed(42)
          p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
                      y=umap_embedding[names(V(igraph.net_marker)),1]) +
            geom_edge_diagonal(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
                              arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
            geom_node_point(aes(fill=celltype), size=7.5, stroke=0, shape=21, alpha=0.75) +
            geom_node_text(aes(label=name), size=2.25,) +
            scale_edge_colour_manual(values=col_arrows,
                                    breaks=waiver(),
                                    labels=labs_plot) +
            scale_fill_manual(values=opts$celltype.colors[celltypes_present],na.value="gray") +
            theme_graph(base_family = 'Helvetica') + theme(legend.position = "none")
          
          # save PDF
          pdf(file.path(args$outdir, sprintf("Score%s/global_network_markerTFs_supercluster%sto%s_%s.pdf",args$min_chip_score, sc, sc, interactions)),
            width  = 6, height = 5.5)

          print(p)
          dev.off()
        } else {
          if (sc == "C"){
            next()
          }
          # Create tbl_graph object for ggraph
          igraph.tbl_marker <- as_tbl_graph(igraph.net_marker) %>%
            activate(nodes) %>%
            mutate(tf=names(V(igraph.net_marker))) %>%
            mutate(degree=igraph::degree(igraph.net_marker)) %>%
            mutate(eigen_centrality=eigen_centrality(igraph.net_marker)$vector) %>%
            activate(edges) %>%
            mutate(sign=ifelse(E(igraph.net_marker)$weight>0,"Positive","Negative"))
          # print(length(E(igraph.net_marker)$weight))
          mask_pos <- (E(igraph.tbl_marker)$weight>0) & !grepl("^Not",E(igraph.tbl_marker)$sign)
          mask_not <- (!((edge_list.dt_marker$from%in%TFs_marker_j) & (edge_list.dt_marker$to%in%TFs_marker_jand1)))
          mask_neg <- (E(igraph.tbl_marker)$weight<0) & !grepl("^Not",E(igraph.tbl_marker)$sign)
          # print(sum(mask_pos & mask_not))
          # print(sum(mask_neg & mask_not))
          E(igraph.tbl_marker)$sign[mask_not] <- paste0("Not in cluster (n=",sum(mask_not),")")
          E(igraph.tbl_marker)$sign[mask_pos & !mask_not] <- paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")")
          E(igraph.tbl_marker)$sign[mask_neg & !mask_not] <- paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")")
          E(igraph.tbl_marker)$sign <- factor(E(igraph.tbl_marker)$sign,levels=c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
                                                                                paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
                                                                                paste0("Not in cluster (n=",sum(mask_not),")")))
          
          labs_plot <- c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
                        paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
                        paste0("Not in cluster (n=",sum(mask_not),")"))
          
          print(labs_plot)
          
          celltypes_markers <- rep(NA,length(V(igraph.tbl_marker)$tf))
          for (i in 1:length(celltypes_markers)){
            if (sum(markers_TF_undupl$TF==V(igraph.tbl_marker)$tf[i])>0){
              celltypes_markers[i] <- markers_TF_undupl$celltype[markers_TF_undupl$TF==V(igraph.tbl_marker)$tf[i]]
            }
          }
          
          if (k == 1){
            celltypes_markers[!((V(igraph.tbl_marker)$tf%in%GRN_coef.dt$tf[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$gene%in%TFs_marker_jand1) & (GRN_coef.dt$beta>0)]) |
                                  (V(igraph.tbl_marker)$tf%in%GRN_coef.dt$gene[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$gene%in%TFs_marker_jand1) & (GRN_coef.dt$beta>0)]))] <- NA
          } else {
            celltypes_markers[!((V(igraph.tbl_marker)$tf%in%GRN_coef.dt$tf[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$gene%in%TFs_marker_jand1) & (GRN_coef.dt$beta<0)]) |
                                  (V(igraph.tbl_marker)$tf%in%GRN_coef.dt$gene[(GRN_coef.dt$tf%in%TFs_marker_j) & (GRN_coef.dt$gene%in%TFs_marker_jand1) & (GRN_coef.dt$beta<0)]))] <- NA
          }
          
          igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(celltype=celltypes_markers)
          
          if (k == 1){
            col_arrows <- c("darkred","gray","gray")
          } else{
            col_arrows <- c("gray","darkblue","gray")
          }
          
          celltypes_present <- sort(as.numeric(unique(celltypes_markers)))
          celltypes_present <- celltypes_present[!is.na(celltypes_present)]
          print(celltypes_present)
          
          set.seed(42)
          p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
                      y=umap_embedding[names(V(igraph.net_marker)),1]) +
            geom_edge_diagonal(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
                              arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
            geom_node_point(aes(fill=celltype), size=7.5, stroke=0, shape=21, alpha=0.75) +
            geom_node_text(aes(label=name), size=2.25,) +
            scale_edge_colour_manual(values=col_arrows,
                                    breaks=waiver(),
                                    labels=labs_plot) +
            scale_fill_manual(values=opts$celltype.colors[celltypes_present],na.value="gray") +
            theme_graph(base_family = 'Helvetica') + theme(legend.position = "none")
          
          # save PDF
          pdf(file.path(args$outdir, sprintf("Score%s/global_network_markerTFs_supercluster%stoC_%s.pdf",args$min_chip_score, sc, interactions)),
            width  = 6, height = 5.5)
            
          print(p)
          dev.off()
        }
      }
    }
  }
}

#########################
## Plot global network ##
#########################

print("Plot global network: celltype (marker TFs)")

# Create node and edge data.frames
node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$tf%in%TFs_marker) & (GRN_coef.dt$gene%in%TFs_marker)][,c("tf","gene","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]

print("Summary edges:")
n_pos_edges <- sum(edge_list.dt_marker$weight>0)
n_neg_edges <- sum(edge_list.dt_marker$weight<0)
print(paste0("positive: ",n_pos_edges))
print(paste0("negative: ",n_neg_edges))

# Create igraph object
igraph.net_marker <- graph_from_data_frame(d = edge_list.dt_marker)

# Create tbl_graph object for ggraph
igraph.tbl_marker <- as_tbl_graph(igraph.net_marker) %>%
  activate(nodes) %>%
  mutate(tf=names(V(igraph.net_marker))) %>%
  mutate(degree=igraph::degree(igraph.net_marker)) %>%
  mutate(eigen_centrality=eigen_centrality(igraph.net_marker)$vector) %>%
  activate(edges) %>%
  mutate(sign=ifelse(E(igraph.net_marker)$weight>0,
                     paste0("Positive (n=",n_pos_edges,")"),paste0("Negative (n=",n_neg_edges,")")))
# igraph.tbl

celltypes_markers <- rep(NA,length(names(V(igraph.net_marker))))
for (i in 1:length(celltypes_markers)){
  if (sum(markers_TF_undupl$TF==names(V(igraph.net_marker))[i])>0){
    celltypes_markers[i] <- markers_TF_undupl$celltype[markers_TF_undupl$TF==names(V(igraph.net_marker))[i]]
  }
}
celltypes_markers <- factor(as.character(celltypes_markers), levels=as.character(1:6))

igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(celltype=celltypes_markers)

set.seed(42)

p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
            y=umap_embedding[names(V(igraph.net_marker)),1]) +
  geom_edge_diagonal(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
                     arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
  geom_node_point(aes(fill=celltype), size=7.5, stroke=0, shape=21, alpha=0.75) +
  geom_node_text(aes(label=name), size=2.25) +
  scale_edge_colour_manual(values=c("darkblue","darkred")) +
  scale_fill_manual(values=opts$celltype.colors,na.value="white") +
  theme_graph(base_family = 'Helvetica') + theme(legend.position = "none")

pdf(file.path(args$outdir,sprintf("Score%s/global_network_markerTFs_celltype.pdf",args$min_chip_score)), 
    width = 6, height = 5.5)
print(p)
dev.off()

p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
            y=umap_embedding[names(V(igraph.net_marker)),1]) +
  geom_edge_diagonal(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
                     arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
  geom_node_point(aes(fill=celltype), size=7.5, stroke=0, shape=21, alpha=0.75) +
  scale_edge_colour_manual(values=c("darkblue","darkred")) +
  scale_fill_manual(values=opts$celltype.colors,na.value="white") +
  theme_graph(base_family = 'Helvetica') + theme(legend.position = "none")

pdf(file.path(args$outdir,sprintf("Score%s/global_network_markerTFs_celltype_blanco.pdf",args$min_chip_score)), 
    width = 6, height = 5.5)
print(p)
dev.off()

tf_celltype_df <- igraph.tbl_marker %>%
  activate(nodes) %>%
  as.data.frame() %>%
  data.table() %>%
  .[, .(tf = name, celltype)]

outfile <- file.path(args$outdir, sprintf("global_network_TFs_celltypes_%s.txt", args$trajectory_name))
fwrite(tf_celltype_df, file = outfile, sep = "\t", quote = FALSE)

##############################################
## Plot network per cluster on general plot ##
##############################################

cluster_names <- c("Myeloid","B cells","T cells","Fibroblasts","Endothelial","Tumor")
for (j in celltypes.to.plot){
  print(paste0("Plot global network for cluster ",j))
  
  # Create node and edge data.frames
  node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
  edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$tf%in%TFs_marker) & (GRN_coef.dt$gene%in%TFs_marker)][,c("tf","gene","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]
  print(dim(edge_list.dt_marker))
  print(head(edge_list.dt_marker))
  
  # Create igraph object
  igraph.net_marker <- graph_from_data_frame(d = edge_list.dt_marker)
  
  # Create tbl_graph object for ggraph
  igraph.tbl_marker <- as_tbl_graph(igraph.net_marker) %>%
    activate(nodes) %>%
    mutate(tf=names(V(igraph.net_marker))) %>%
    mutate(degree=igraph::degree(igraph.net_marker)) %>%
    mutate(eigen_centrality=eigen_centrality(igraph.net_marker)$vector) %>%
    activate(edges)
  
  expression_values <- rep(NA,length(V(igraph.tbl_marker)$tf))
  for (i in 1:length(V(igraph.tbl_marker))){
    expression_values[i] <- logcounts(rna_pseudobulk.sce[V(igraph.tbl_marker)$tf[i],as.character(j)])
  }
  
  z_score_values <- z_score(expression_values)
  
  expression_perc <- rep(NA,length(V(igraph.tbl_marker)$tf))
  for (i in 1:length(V(igraph.tbl_marker))){
    expression_perc[i] <- logcounts(rna_pseudobulk.sce[V(igraph.tbl_marker)$tf[i],as.character(j)])/
      max(logcounts(rna_pseudobulk.sce[V(igraph.tbl_marker)$tf[i],1:5]))
  }
  
  # reduction <- 0.10
  # red_cutoff <- quantile(z_score_values,c(reduction,1-reduction))
  # print(red_cutoff)
  red_cutoff <- c(-1,1.5)
  # print(red_cutoff)
  z_score_red <- z_score_values
  z_score_red[z_score_red<red_cutoff[1]] <- red_cutoff[1]
  z_score_red[z_score_red>red_cutoff[2]] <- red_cutoff[2]
  
  igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(expression=expression_values)
  igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(z_score=z_score_values)
  igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(z_score_red=z_score_red)
  igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(expression_perc=expression_perc)
  
  print(length(expression_values))
  print(length(expression_perc))
  print(length(z_score_red))
  print(length(z_score_values))
  print(str(igraph.tbl_marker))

  col_arrows <- "grey"
  
  set.seed(42)
  p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
              y=umap_embedding[names(V(igraph.net_marker)),1]) +
    geom_edge_diagonal(edge_colour="grey", edge_alpha=0.85, edge_width=0.15,
                       arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
    geom_node_point(aes(fill=expression), size=7.5, stroke=0, shape=21, alpha=0.75) +
    geom_node_text(aes(label=name), size=2.25,) +
    scale_colour_gradientn(colors=colorRampPalette(brewer.pal(9, "YlOrRd"))(100),aesthetics = "fill") +
    theme_graph(base_family = 'Helvetica') + theme(legend.position = "none") #+ 
  #annotate("text",  x=Inf, y = Inf, label = cluster_names[j], vjust=1, hjust=1)
  

  pdf(file.path(args$outdir,sprintf("Score%s/global_network_gradient_expr_cluster%s.pdf",args$min_chip_score,j)), 
      width = 6, height = 5.5)

  print(p)
  dev.off()
  
  set.seed(42)
  p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
              y=umap_embedding[names(V(igraph.net_marker)),1]) +
    geom_edge_diagonal(edge_colour="grey", edge_alpha=0.85, edge_width=0.15,
                       arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
    geom_node_point(aes(fill=z_score), size=7.5, stroke=0, shape=21, alpha=0.75) +
    geom_node_text(aes(label=name), size=2.25,) +
    scale_colour_gradientn(colors=colorRampPalette(brewer.pal(9, "YlOrRd"))(100),aesthetics = "fill") +
    theme_graph(base_family = 'Helvetica') + theme(legend.position = "none") #+ 
  #annotate("text",  x=Inf, y = Inf, label = cluster_names[j], vjust=1, hjust=1)
  

  pdf(file.path(args$outdir,sprintf("Score%s/global_network_gradient_zscore_cluster%s.pdf",args$min_chip_score,j)), 
      width = 6, height = 5.5)
  print(p)
  dev.off()
  
  set.seed(42)
  p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
              y=umap_embedding[names(V(igraph.net_marker)),1]) +
    geom_edge_diagonal(edge_colour="grey", edge_alpha=0.85, edge_width=0.15,
                       arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
    geom_node_point(aes(fill=z_score_red), size=7.5, stroke=0, shape=21, alpha=0.75) +
    geom_node_text(aes(label=name), size=2.25,) +
    scale_colour_gradientn(colors=colorRampPalette(brewer.pal(9, "YlOrRd"))(100),aesthetics = "fill") +
    theme_graph(base_family = 'Helvetica') + theme(legend.position = "none") #+ 
  #annotate("text",  x=Inf, y = Inf, label = cluster_names[j], vjust=1, hjust=1)
  
  
  pdf(file.path(args$outdir,sprintf("Score%s/global_network_gradient_zscore_red_cluster%s.pdf",args$min_chip_score,j)), 
      width = 6, height = 5.5)
  print(p)
  dev.off()
  
  set.seed(42)
  p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
              y=umap_embedding[names(V(igraph.net_marker)),1]) +
    geom_edge_diagonal(edge_colour="grey", edge_alpha=0.85, edge_width=0.15,
                       arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
    geom_node_point(aes(fill=expression_perc), size=7.5, stroke=0, shape=21, alpha=0.75) +
    geom_node_text(aes(label=name), size=2.25,) +
    scale_colour_gradientn(colors=colorRampPalette(brewer.pal(9, "YlOrRd"))(100),aesthetics = "fill") +
    theme_graph(base_family = 'Helvetica') + theme(legend.position = "none") #+ 
    #annotate("text",  x=Inf, y = Inf, label = cluster_names[j], vjust=1, hjust=1)
  
  
  pdf(file.path(args$outdir,sprintf("Score%s/global_network_gradient_expr_perc_cluster%s.pdf",args$min_chip_score,j)), 
      width = 6, height = 5.5)
  
  print(p)
  dev.off()
}


#############################################
## Plot network per sample on general plot ##
#############################################

sample_names <- colnames(rna_pseudobulk_sample.sce)[sort(as.numeric(gsub("^d","",colnames(rna_pseudobulk_sample.sce))),index.return=T)$ix]
print(sample_names)
for (j in sample_names){
  print(paste0("Plot global network for sample ",j))
  
  # Create node and edge data.frames
  node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
  edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$tf%in%TFs_marker) & (GRN_coef.dt$gene%in%TFs_marker)][,c("tf","gene","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]
  print(dim(edge_list.dt_marker))
  print(head(edge_list.dt_marker))
  
  # Create igraph object
  igraph.net_marker <- graph_from_data_frame(d = edge_list.dt_marker)
  
  # Create tbl_graph object for ggraph
  igraph.tbl_marker <- as_tbl_graph(igraph.net_marker) %>%
    activate(nodes) %>%
    mutate(tf=names(V(igraph.net_marker))) %>%
    mutate(degree=igraph::degree(igraph.net_marker)) %>%
    mutate(eigen_centrality=eigen_centrality(igraph.net_marker)$vector) %>%
    activate(edges)
  
  expression_values <- rep(NA,length(V(igraph.tbl_marker)$tf))
  for (i in 1:length(V(igraph.tbl_marker))){
    expression_values[i] <- logcounts(rna_pseudobulk_sample.sce[V(igraph.tbl_marker)$tf[i],j])
  }
  
  z_score_values <- z_score(expression_values)
  
  expression_perc <- rep(NA,length(V(igraph.tbl_marker)$tf))
  for (i in 1:length(V(igraph.tbl_marker))){
    expression_perc[i] <- logcounts(rna_pseudobulk_sample.sce[V(igraph.tbl_marker)$tf[i],j])/
      max(logcounts(rna_pseudobulk_sample.sce[V(igraph.tbl_marker)$tf[i],]))
  }
  
  # reduction <- 0.10
  # red_cutoff <- quantile(z_score_values,c(reduction,1-reduction))
  red_cutoff <- c(-1,1.5)
  z_score_red <- z_score_values
  z_score_red[z_score_red<red_cutoff[1]] <- red_cutoff[1]
  z_score_red[z_score_red>red_cutoff[2]] <- red_cutoff[2]
  
  igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(expression=expression_values)
  igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(z_score=z_score_values)
  igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(z_score_red=z_score_red)
  igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(expression_perc=expression_perc)
  
  col_arrows <- "grey"
  
  set.seed(42)
  p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
              y=umap_embedding[names(V(igraph.net_marker)),1]) +
    geom_edge_diagonal(edge_colour="grey", edge_alpha=0.85, edge_width=0.15,
                       arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
    geom_node_point(aes(fill=expression), size=7.5, stroke=0, shape=21, alpha=0.75) +
    geom_node_text(aes(label=name), size=2.25,) +
    scale_colour_gradientn(colors=colorRampPalette(brewer.pal(9, "YlOrRd"))(100),aesthetics = "fill") +
    theme_graph(base_family = 'Helvetica') + theme(legend.position = "none")#+ 
  # annotate("text",  x=Inf, y = Inf, label = j, vjust=1, hjust=1)
  
  
  pdf(file.path(args$outdir,sprintf("Score%s/global_network_gradient_expr_sample%s.pdf",args$min_chip_score,j)), 
      width = 6, height = 5.5)
  
  print(p)
  dev.off()
  
  set.seed(42)
  p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
              y=umap_embedding[names(V(igraph.net_marker)),1]) +
    geom_edge_diagonal(edge_colour="grey", edge_alpha=0.85, edge_width=0.15,
                       arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
    geom_node_point(aes(fill=z_score), size=7.5, stroke=0, shape=21, alpha=0.75) +
    geom_node_text(aes(label=name), size=2.25,) +
    scale_colour_gradientn(colors=colorRampPalette(brewer.pal(9, "YlOrRd"))(100),aesthetics = "fill") +
    theme_graph(base_family = 'Helvetica') + theme(legend.position = "none")#+ 
  # annotate("text",  x=Inf, y = Inf, label = j, vjust=1, hjust=1)
  
  
  pdf(file.path(args$outdir,sprintf("Score%s/global_network_gradient_zscore_sample%s.pdf",args$min_chip_score,j)), 
      width = 6, height = 5.5)
  print(p)
  dev.off()
  
  set.seed(42)
  p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
              y=umap_embedding[names(V(igraph.net_marker)),1]) +
    geom_edge_diagonal(edge_colour="grey", edge_alpha=0.85, edge_width=0.15,
                       arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
    geom_node_point(aes(fill=z_score_red), size=7.5, stroke=0, shape=21, alpha=0.75) +
    geom_node_text(aes(label=name), size=2.25,) +
    scale_colour_gradientn(colors=colorRampPalette(brewer.pal(9, "YlOrRd"))(100),aesthetics = "fill") +
    theme_graph(base_family = 'Helvetica') + theme(legend.position = "none") #+ 
  # annotate("text",  x=Inf, y = Inf, label = j, vjust=1, hjust=1)
  
  
  pdf(file.path(args$outdir,sprintf("Score%s/global_network_gradient_zscore_red_sample%s.pdf",args$min_chip_score,j)), 
      width = 6, height = 5.5)
  print(p)
  dev.off()
  
  set.seed(42)
  p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
              y=umap_embedding[names(V(igraph.net_marker)),1]) +
    geom_edge_diagonal(edge_colour="grey", edge_alpha=0.85, edge_width=0.15,
                       arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
    geom_node_point(aes(fill=expression_perc), size=7.5, stroke=0, shape=21, alpha=0.75) +
    geom_node_text(aes(label=name), size=2.25,) +
    scale_colour_gradientn(colors=colorRampPalette(brewer.pal(9, "YlOrRd"))(100),aesthetics = "fill") +
    theme_graph(base_family = 'Helvetica') + theme(legend.position = "none") #+ 
    # annotate("text",  x=Inf, y = Inf, label = j, vjust=1, hjust=1)
  
  
  pdf(file.path(args$outdir,sprintf("Score%s/global_network_gradient_expr_perc_sample%s.pdf",args$min_chip_score,j)), 
      width = 6, height = 5.5)
  
  print(p)
  dev.off()
}


######################
## Completion token ##
######################

file.create(file.path(args$outdir,sprintf("Score%s/completed.txt",args$min_chip_score)))

##################
# Print warnings #
##################

warnings()