##############################################
##                                          ##
##  Plot_UMAP_GRN_metacells_trajectories.R  ##
##                                          ##
##############################################

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
p$add_argument('--rna_pseudobulk_cluster.sce',       type="character",                help='SingleCellExperiment (pseudobulk) based on cluster')
p$add_argument('--rna_pseudobulk_sample.sce',       type="character",                help='SingleCellExperiment (pseudobulk) based on sample')
p$add_argument('--grn_coef',       type="character",                help='GRN coefficients')
p$add_argument('--min_coef',  type="double",            default=0.25,      help='Minimal coefficient for linear regression')
p$add_argument('--max_distance',  type="integer",            default=1e5,      help='Maximum distance for a linkage between a peak and a gene')
p$add_argument('--max_pval',  type="double",            default=0.10,      help='Maximal pvalue for linear regression')
p$add_argument('--min_chip_score',  type="double",            default=0.25,      help='Minimal chip score')
p$add_argument('--TFs_to_plot',       type="character",                help='File containing TF of interest to plot')
p$add_argument('--force',       type="logical",  default=F         ,      help='Should replotting be forced?')
p$add_argument('--seed',       type="integer",  default=42         ,      help='Random seed')
p$add_argument('--tf2gene_virtual_chip',       type="character",                help='Links between tfs and genes based on virtual chipseq')
p$add_argument('--sign_connections',       type="character",                help='Metadata metacells')
p$add_argument('--outdir',       type="character",                help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

# # I/O
# args <- list()
# args$trajectory_name <- "N2P"
# args$rna_metacells.sce <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/metacells/SingleCellExperiment_metacells_nodiff.rds"
# args$rna_pseudobulk_cluster.sce <-  "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/SingleCellExperiment_pseudobulk.rds"
# args$rna_pseudobulk_sample.sce <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/SingleCellExperiment_dayofsampling.rds"
# args$min_coef <- 0.25
# args$max_distance <- 5e4
# args$max_pval <- 0.10
# args$min_chip_score <- 0.06
# args$TFs_to_plot <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/chromvar_chip/GRN_TFs_to_plot.txt"
# args$force <- TRUE
# args$seed <- 42
# args$tf2gene_virtual_chip <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/gene_regulatory_network/metacells/CISBP/TF2gene_after_virtual_chip.txt.gz"
# args$sign_connections <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/NicheNet/GRN_coef_signalling_alt5_to_add_merged.txt"
# args$grn_coef <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/gene_regulatory_network/metacells/trajectories/N2P/TFandSign/alt5/merged/global_chip_GRN_sign_alt5_coef_score0.06_merged.txt.gz"
# args$outdir <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/gene_regulatory_network/metacells/trajectories/N2P/TFandSign/alt5/merged/UMAP_fig"

#####################
## Define settings ##
#####################

args$trajectory <- sprintf("/data/homes/louisc/Project_Babraham/RNA_ATAC/trajectories/%s/%s_trajectory.txt.gz",
                           args$trajectory_name,args$trajectory_name)

if(!dir.exists(file.path(args$outdir,sprintf("Score%s",args$min_chip_score)))){
  dir.create(args$outdir, showWarnings = F)
  dir.create(file.path(args$outdir,sprintf("Score%s",args$min_chip_score)))
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
print(nrow(markers_TF_undupl))
print(head(markers_TF_undupl))

#####################
## Load signalling ##
#####################

print("Loading signalling interactions")

sign2TF <- fread(args$sign_connections)
print(dim(sign2TF))
print(head(sign2TF))

###########################
## construct node object ##
###########################

nodes <- markers_TF_undupl 

lines_to_add <- sign2TF[!duplicated(sign2TF$from),c("from","fr_max_expr")]
colnames(lines_to_add) <- c("TF","celltype")
nodes <- rbind(nodes,lines_to_add)

##############################
## Load RNA expression data ##
##############################

# SingleCellExperiment at metacell resolution
rna_metacells.sce <- readRDS(args$rna_metacells.sce)
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

GRN_coef.dt <- fread(args$grn_coef) 

##########################
## Filter TFs and genes ##
##########################

TFs_marker <- nodes$TF
print(length(TFs_marker))
print(TFs_marker)

TFs <- TFs_marker
print("TFs")
head(TFs)
print("GRN_coef.dt")
head(GRN_coef.dt)
GRN_coef.dt <- GRN_coef.dt[sender%in%TFs & target%in%TFs]
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

betas_mtx <- cor(t(logcounts(rna_metacells.sce[rownames(rna_metacells.sce)%in%TFs,])))

betas_mtx[1:6,1:6]
betas_mtx <- betas_mtx[TFs,TFs]
betas_mtx[1:6,1:6]

# Run
set.seed(args$seed)
umap_embedding <- uwot::umap(betas_mtx, spread = 0.5, min_dist = 0.5, n_neighbors=30, pca = 10)
# umap_embedding <- uwot::umap(betas_mtx, pca = 20)

print("UMAP")
# umap_embedding[,1] <- -umap_embedding[,1]
# umap_embedding[,2] <- -umap_embedding[,2]
print(dim(umap_embedding))
print(head(umap_embedding))

#########################
## Plot global network ##
#########################

print("Plot global network: celltype (marker TFs)")

# Create node and edge data.frames
node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$sender%in%TFs_marker) & (GRN_coef.dt$target%in%TFs_marker)][,c("sender","target","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]

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
  if (sum(nodes$TF==names(V(igraph.net_marker))[i])>0){
    celltypes_markers[i] <- nodes$celltype[nodes$TF==names(V(igraph.net_marker))[i]]
  }
}
stroke_colors <- opts$celltype.colors[as.numeric(celltypes_markers)]
celltypes_markers[names(V(igraph.net_marker))%in%sign2TF$from] <- 7
stroke_colors[celltypes_markers!=7] <- "#000000"
celltypes_markers <- factor(as.character(celltypes_markers),levels=as.character(1:5))

igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(celltype=celltypes_markers)

set.seed(42)

# with TF names
p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
            y=umap_embedding[names(V(igraph.net_marker)),1]) +
  geom_edge_diagonal(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
                     arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
  geom_node_point(aes(fill=celltype), size=7.5, stroke=0.5, shape=21, alpha=0.75, color = stroke_colors) +
  geom_node_text(aes(label=name), size=2.25) +
  scale_edge_colour_manual(values=c("darkblue","darkred")) +
  scale_fill_manual(values=opts$celltype.colors,na.value="white") +
  theme_graph(base_family = 'Helvetica') + theme(legend.position = "none")

pdf(file.path(args$outdir,sprintf("Score%s/global_network_markerTFs_celltype.pdf",args$min_chip_score)), 
    width = 6, height = 5.5)
print(p)
dev.off()

# blanco
p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
            y=umap_embedding[names(V(igraph.net_marker)),1]) +
  geom_edge_diagonal(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
                     arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
  geom_node_point(aes(fill=celltype), size=7.5, stroke=0.5, shape=21, alpha=0.75, color = stroke_colors) +
  scale_edge_colour_manual(values=c("darkblue","darkred")) +
  scale_fill_manual(values=opts$celltype.colors,na.value="white") +
  theme_graph(base_family = 'Helvetica') + theme(legend.position = "none")

pdf(file.path(args$outdir,sprintf("Score%s/global_network_markerTFs_celltype_blanco.pdf",args$min_chip_score)), 
    width = 6, height = 5.5)
print(p)
dev.off()

#############################################
## Plot signalling network on general plot ##
#############################################

print("Plot signalling network on general plot")

interactions <- c("activating","repressive","both")
for (int in interactions){
  print(int)

  # Create node and edge data.frames
  node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
  edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$sender%in%TFs_marker) & (GRN_coef.dt$target%in%TFs_marker)][,c("sender","target","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]

  n_pos_edges <- sum(edge_list.dt_marker$weight>0)
  n_neg_edges <- sum(edge_list.dt_marker$weight<0)

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

  # print(length(E(igraph.net_marker)$weight))
  mask_pos <- (E(igraph.tbl_marker)$weight>0) & !grepl("^Not",E(igraph.tbl_marker)$sign)
  mask_not <- (!((edge_list.dt_marker$from%in%sign2TF$from) | (edge_list.dt_marker$to%in%sign2TF$from)))
  mask_neg <- (E(igraph.tbl_marker)$weight<0) & !grepl("^Not",E(igraph.tbl_marker)$sign)

  E(igraph.tbl_marker)$sign[mask_not] <- paste0("Not signalling (n=",sum(mask_not),")")
  E(igraph.tbl_marker)$sign[mask_pos & !mask_not] <- paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")")
  E(igraph.tbl_marker)$sign[mask_neg & !mask_not] <- paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")")
  E(igraph.tbl_marker)$sign <- factor(E(igraph.tbl_marker)$sign,levels=c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
                                                                          paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
                                                                          paste0("Not signalling (n=",sum(mask_not),")")))
      
  labs_plot <- c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
                  paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
                  paste0("Not signalling (n=",sum(mask_not),")"))

                  
  print("Summary edges:")
  print(paste0("positive: ",sum(mask_pos & !mask_not)))
  print(paste0("negative: ",sum(mask_neg & !mask_not)))
  print(paste0("irrelevant: ",sum(mask_not)))

  if (int == "activating"){
    col_arrows <- c("darkred","grey","gray") 
  } else if (int == "repressive"){
    col_arrows <- c("grey","darkblue","gray") 
  } else {
    col_arrows <- c("darkred","darkblue","gray")
  }

  celltypes_markers <- rep(NA,length(names(V(igraph.net_marker))))
  for (i in 1:length(celltypes_markers)){
    if (sum(nodes$TF==names(V(igraph.net_marker))[i])>0){
      celltypes_markers[i] <- nodes$celltype[nodes$TF==names(V(igraph.net_marker))[i]]
    }
  }
  stroke_colors <- opts$celltype.colors[as.numeric(celltypes_markers)]
  celltypes_markers[names(V(igraph.net_marker))%in%sign2TF$from] <- 7
  celltypes_markers <- factor(as.character(celltypes_markers),levels=as.character(1:7))

  # only keep signalling links 
  celltypes_markers[!(names(V(igraph.net_marker))%in%unique(c(sign2TF$from,edge_list.dt_marker$to[(edge_list.dt_marker$from%in%sign2TF$from)],edge_list.dt_marker$from[(edge_list.dt_marker$to%in%sign2TF$from)])))] <- NA

  # remove activating/repressive links

  print(celltypes_markers)
  stroke_colors[(celltypes_markers!=7) | is.na(celltypes_markers)] <- "#000000"

  igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(celltype=celltypes_markers)

  celltypes_present <- sort(as.numeric(unique(celltypes_markers)))
  celltypes_present <- celltypes_present[!is.na(celltypes_present)]
  print(celltypes_present)

  if (int == "activating"){
    mask_to_front <- mask_pos & !mask_not
  } else if (int == "repressive"){
    mask_to_front <- mask_neg & !mask_not
  } else {
    mask_to_front <- (mask_pos | mask_neg) & !mask_not
  }
  igraph.tbl_marker_plot = graph_from_data_frame(d=as_data_frame(igraph.tbl_marker,what="edges") %>%
    arrange(c(mask_to_front)+1),
    vertices=as_data_frame(igraph.tbl_marker,what="vertices"))
  
  set.seed(42)
  p <- ggraph(igraph.tbl_marker_plot, x=umap_embedding[names(V(igraph.tbl_marker_plot)),2],
              y=umap_embedding[names(V(igraph.tbl_marker_plot)),1]) +
    geom_edge_diagonal(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
                    arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
    geom_node_point(aes(fill=celltype), size=7.5, stroke=0, shape=21, alpha=0.75, color = stroke_colors) +
    geom_node_text(aes(label=name), size=2.25,) +
    scale_edge_colour_manual(values=col_arrows,
                              breaks=waiver(),
                              labels=labs_plot) +
    scale_fill_manual(values=opts$celltype.colors[celltypes_present],na.value="gray") +
    theme_graph(base_family = 'Helvetica') + theme(legend.position = "none")
      
  pdf(file.path(args$outdir,sprintf("Score%s/global_network_markerTFs_signalling_%s.pdf",args$min_chip_score,int)), 
      width = 6, height = 5.5)
  print(p)
  dev.off()

  p <- ggraph(igraph.tbl_marker_plot, x=umap_embedding[names(V(igraph.tbl_marker_plot)),2],
              y=umap_embedding[names(V(igraph.tbl_marker_plot)),1]) +
    geom_edge_diagonal(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
                    arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
    geom_node_point(aes(fill=celltype), size=7.5, stroke=0, shape=21, alpha=0.75, color = stroke_colors) +
    scale_edge_colour_manual(values=col_arrows,
                              breaks=waiver(),
                              labels=labs_plot) +
    scale_fill_manual(values=opts$celltype.colors[celltypes_present],na.value="gray") +
    theme_graph(base_family = 'Helvetica') + theme(legend.position = "none")
      
  pdf(file.path(args$outdir,sprintf("Score%s/global_network_markerTFs_signalling_%s_blanco.pdf",args$min_chip_score,int)), 
      width = 6, height = 5.5)
  print(p)
  dev.off()
}

####################################################
## Plot signalling target network on general plot ##
####################################################

print("Plot to signalling target network on general plot")

interactions <- c("activating","repressive","both")
for (int in interactions){
  print(int)

  # Create node and edge data.frames
  node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
  edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$sender%in%TFs_marker) & (GRN_coef.dt$target%in%TFs_marker)][,c("sender","target","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]

  n_pos_edges <- sum(edge_list.dt_marker$weight>0)
  n_neg_edges <- sum(edge_list.dt_marker$weight<0)

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

  # print(length(E(igraph.net_marker)$weight))
  mask_pos <- (E(igraph.tbl_marker)$weight>0) & !grepl("^Not",E(igraph.tbl_marker)$sign)
  mask_not <- (!((edge_list.dt_marker$to%in%sign2TF$from)))
  mask_neg <- (E(igraph.tbl_marker)$weight<0) & !grepl("^Not",E(igraph.tbl_marker)$sign)

  E(igraph.tbl_marker)$sign[mask_not] <- paste0("Not signalling (n=",sum(mask_not),")")
  E(igraph.tbl_marker)$sign[mask_pos & !mask_not] <- paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")")
  E(igraph.tbl_marker)$sign[mask_neg & !mask_not] <- paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")")
  E(igraph.tbl_marker)$sign <- factor(E(igraph.tbl_marker)$sign,levels=c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
                                                                          paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
                                                                          paste0("Not signalling (n=",sum(mask_not),")")))
      
  labs_plot <- c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
                  paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
                  paste0("Not signalling (n=",sum(mask_not),")"))

                  
  print("Summary edges:")
  print(paste0("positive: ",sum(mask_pos & !mask_not)))
  print(paste0("negative: ",sum(mask_neg & !mask_not)))
  print(paste0("irrelevant: ",sum(mask_not)))

  if (int == "activating"){
    col_arrows <- c("darkred","grey","gray") 
  } else if (int == "repressive"){
    col_arrows <- c("grey","darkblue","gray") 
  } else {
    col_arrows <- c("darkred","darkblue","gray")
  }

  celltypes_markers <- rep(NA,length(names(V(igraph.net_marker))))
  for (i in 1:length(celltypes_markers)){
    if (sum(nodes$TF==names(V(igraph.net_marker))[i])>0){
      celltypes_markers[i] <- nodes$celltype[nodes$TF==names(V(igraph.net_marker))[i]]
    }
  }
  stroke_colors <- opts$celltype.colors[as.numeric(celltypes_markers)]
  celltypes_markers[names(V(igraph.net_marker))%in%sign2TF$from] <- 7
  celltypes_markers <- factor(as.character(celltypes_markers),levels=as.character(1:7))

  # only keep signalling links 
  celltypes_markers[!(names(V(igraph.net_marker))%in%unique(c(sign2TF$from,edge_list.dt_marker$to[(edge_list.dt_marker$from%in%sign2TF$from)],edge_list.dt_marker$from[(edge_list.dt_marker$to%in%sign2TF$from)])))] <- NA

  # remove activating/repressive links
  if (int == "activating"){
    set_to_remove <- unique(c(edge_list.dt_marker$from[mask_neg | mask_not],edge_list.dt_marker$to[mask_neg | mask_not]))
    set_to_remove <- set_to_remove[!(set_to_remove%in%unique(c(edge_list.dt_marker$from[mask_pos & !mask_not],edge_list.dt_marker$to[mask_pos & !mask_not])))]
    celltypes_markers[names(V(igraph.net_marker))%in%set_to_remove] <- NA
  } else if (int == "repressive"){
    set_to_remove <- unique(c(edge_list.dt_marker$from[mask_pos | mask_not],edge_list.dt_marker$to[mask_pos | mask_not]))
    set_to_remove <- set_to_remove[!(set_to_remove%in%unique(c(edge_list.dt_marker$from[mask_neg & !mask_not],edge_list.dt_marker$to[mask_neg & !mask_not])))]
    celltypes_markers[names(V(igraph.net_marker))%in%set_to_remove] <- NA
  } else {
    set_to_remove <- unique(c(edge_list.dt_marker$from[mask_not],edge_list.dt_marker$to[mask_not]))
    set_to_remove <- set_to_remove[!(set_to_remove%in%unique(c(edge_list.dt_marker$from[(mask_pos | mask_neg) & !mask_not],edge_list.dt_marker$to[(mask_pos | mask_neg) & !mask_not])))]
    celltypes_markers[names(V(igraph.net_marker))%in%set_to_remove] <- NA
  }
  print(celltypes_markers)
  stroke_colors[(celltypes_markers!=7) | is.na(celltypes_markers)] <- "#000000"

  igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(celltype=celltypes_markers)

  celltypes_present <- sort(as.numeric(unique(celltypes_markers)))
  celltypes_present <- celltypes_present[!is.na(celltypes_present)]
  print(celltypes_present)

  if (int == "activating"){
    mask_to_front <- mask_pos & !mask_not
  } else if (int == "repressive"){
    mask_to_front <- mask_neg & !mask_not
  } else {
    mask_to_front <- (mask_pos | mask_neg) & !mask_not
  }
  igraph.tbl_marker_plot = graph_from_data_frame(d=as_data_frame(igraph.tbl_marker,what="edges") %>%
    arrange(c(mask_to_front)+1),
    vertices=as_data_frame(igraph.tbl_marker,what="vertices"))
  
  set.seed(42)
  p <- ggraph(igraph.tbl_marker_plot, x=umap_embedding[names(V(igraph.tbl_marker_plot)),2],
              y=umap_embedding[names(V(igraph.tbl_marker_plot)),1]) +
    geom_edge_diagonal(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
                    arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
    geom_node_point(aes(fill=celltype), size=7.5, stroke=0, shape=21, alpha=0.75, color = stroke_colors) +
    geom_node_text(aes(label=name), size=2.25,) +
    scale_edge_colour_manual(values=col_arrows,
                              breaks=waiver(),
                              labels=labs_plot) +
    scale_fill_manual(values=opts$celltype.colors[celltypes_present],na.value="gray") +
    theme_graph(base_family = 'Helvetica') + theme(legend.position = "none")
      
  pdf(file.path(args$outdir,sprintf("Score%s/global_network_markerTFs_signalling_target_%s.pdf",args$min_chip_score,int)), 
      width = 6, height = 5.5)
  print(p)
  dev.off()

  p <- ggraph(igraph.tbl_marker_plot, x=umap_embedding[names(V(igraph.tbl_marker_plot)),2],
            y=umap_embedding[names(V(igraph.tbl_marker_plot)),1]) +
  geom_edge_diagonal(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
                  arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
  geom_node_point(aes(fill=celltype), size=7.5, stroke=0, shape=21, alpha=0.75, color = stroke_colors) +
  scale_edge_colour_manual(values=col_arrows,
                            breaks=waiver(),
                            labels=labs_plot) +
  scale_fill_manual(values=opts$celltype.colors[celltypes_present],na.value="gray") +
  theme_graph(base_family = 'Helvetica') + theme(legend.position = "none")
      
  pdf(file.path(args$outdir,sprintf("Score%s/global_network_markerTFs_signalling_target_%s_blanco.pdf",args$min_chip_score,int)), 
      width = 6, height = 5.5)
  print(p)
  dev.off()
}

####################################################
## Plot signalling sender network on general plot ##
####################################################

print("Plot signalling sender network on general plot")

interactions <- c("activating","repressive","both")
for (int in interactions){
  print(int)

  # Create node and edge data.frames
  node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
  edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$sender%in%TFs_marker) & (GRN_coef.dt$target%in%TFs_marker)][,c("sender","target","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]

  n_pos_edges <- sum(edge_list.dt_marker$weight>0)
  n_neg_edges <- sum(edge_list.dt_marker$weight<0)

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

  # print(length(E(igraph.net_marker)$weight))
  mask_pos <- (E(igraph.tbl_marker)$weight>0) & !grepl("^Not",E(igraph.tbl_marker)$sign)
  mask_not <- (!((edge_list.dt_marker$from%in%sign2TF$from)))
  mask_neg <- (E(igraph.tbl_marker)$weight<0) & !grepl("^Not",E(igraph.tbl_marker)$sign)

  E(igraph.tbl_marker)$sign[mask_not] <- paste0("Not signalling (n=",sum(mask_not),")")
  E(igraph.tbl_marker)$sign[mask_pos & !mask_not] <- paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")")
  E(igraph.tbl_marker)$sign[mask_neg & !mask_not] <- paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")")
  E(igraph.tbl_marker)$sign <- factor(E(igraph.tbl_marker)$sign,levels=c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
                                                                          paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
                                                                          paste0("Not signalling (n=",sum(mask_not),")")))
      
  labs_plot <- c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
                  paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
                  paste0("Not signalling (n=",sum(mask_not),")"))

                  
  print("Summary edges:")
  print(paste0("positive: ",sum(mask_pos & !mask_not)))
  print(paste0("negative: ",sum(mask_neg & !mask_not)))
  print(paste0("irrelevant: ",sum(mask_not)))

  if (int == "activating"){
    col_arrows <- c("darkred","grey","gray") 
  } else if (int == "repressive"){
    col_arrows <- c("grey","darkblue","gray") 
  } else {
    col_arrows <- c("darkred","darkblue","gray")
  }

  celltypes_markers <- rep(NA,length(names(V(igraph.net_marker))))
  for (i in 1:length(celltypes_markers)){
    if (sum(nodes$TF==names(V(igraph.net_marker))[i])>0){
      celltypes_markers[i] <- nodes$celltype[nodes$TF==names(V(igraph.net_marker))[i]]
    }
  }
  stroke_colors <- opts$celltype.colors[as.numeric(celltypes_markers)]
  celltypes_markers[names(V(igraph.net_marker))%in%sign2TF$from] <- 7
  celltypes_markers <- factor(as.character(celltypes_markers),levels=as.character(1:7))

  # only keep signalling links 
  celltypes_markers[!(names(V(igraph.net_marker))%in%unique(c(sign2TF$from,edge_list.dt_marker$to[(edge_list.dt_marker$from%in%sign2TF$from)])))] <- NA

  # remove activating/repressive links
  if (int == "activating"){
    set_to_remove <- unique(c(edge_list.dt_marker$from[mask_neg | mask_not],edge_list.dt_marker$to[mask_neg | mask_not]))
    set_to_remove <- set_to_remove[!(set_to_remove%in%unique(c(edge_list.dt_marker$from[mask_pos & !mask_not],edge_list.dt_marker$to[mask_pos & !mask_not])))]
    celltypes_markers[names(V(igraph.net_marker))%in%set_to_remove] <- NA
  } else if (int == "repressive"){
    set_to_remove <- unique(c(edge_list.dt_marker$from[mask_pos | mask_not],edge_list.dt_marker$to[mask_pos | mask_not]))
    set_to_remove <- set_to_remove[!(set_to_remove%in%unique(c(edge_list.dt_marker$from[mask_neg & !mask_not],edge_list.dt_marker$to[mask_neg & !mask_not])))]
    celltypes_markers[names(V(igraph.net_marker))%in%set_to_remove] <- NA
  } else {
    set_to_remove <- unique(c(edge_list.dt_marker$from[mask_not],edge_list.dt_marker$to[mask_not]))
    set_to_remove <- set_to_remove[!(set_to_remove%in%unique(c(edge_list.dt_marker$from[(mask_pos | mask_neg) & !mask_not],edge_list.dt_marker$to[(mask_pos | mask_neg) & !mask_not])))]
    celltypes_markers[names(V(igraph.net_marker))%in%set_to_remove] <- NA
  }
  print(celltypes_markers)
  stroke_colors[(celltypes_markers!=7) | is.na(celltypes_markers)] <- "#000000"

  igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(celltype=celltypes_markers)

  celltypes_present <- sort(as.numeric(unique(celltypes_markers)))
  celltypes_present <- celltypes_present[!is.na(celltypes_present)]
  print(celltypes_present)

  if (int == "activating"){
    mask_to_front <- mask_pos & !mask_not
  } else if (int == "repressive"){
    mask_to_front <- mask_neg & !mask_not
  } else {
    mask_to_front <- (mask_pos | mask_neg) & !mask_not
  }
  igraph.tbl_marker_plot = graph_from_data_frame(d=as_data_frame(igraph.tbl_marker,what="edges") %>%
    arrange(c(mask_to_front)+1),
    vertices=as_data_frame(igraph.tbl_marker,what="vertices"))
  
  set.seed(42)
  p <- ggraph(igraph.tbl_marker_plot, x=umap_embedding[names(V(igraph.tbl_marker_plot)),2],
              y=umap_embedding[names(V(igraph.tbl_marker_plot)),1]) +
    geom_edge_diagonal(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
                    arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
    geom_node_point(aes(fill=celltype), size=7.5, stroke=0, shape=21, alpha=0.75, color = stroke_colors) +
    geom_node_text(aes(label=name), size=2.25,) +
    scale_edge_colour_manual(values=col_arrows,
                              breaks=waiver(),
                              labels=labs_plot) +
    scale_fill_manual(values=opts$celltype.colors[celltypes_present],na.value="gray") +
    theme_graph(base_family = 'Helvetica') + theme(legend.position = "none")
      
  pdf(file.path(args$outdir,sprintf("Score%s/global_network_markerTFs_signalling_sender_%s.pdf",args$min_chip_score,int)), 
      width = 6, height = 5.5)
  print(p)
  dev.off()

    p <- ggraph(igraph.tbl_marker_plot, x=umap_embedding[names(V(igraph.tbl_marker_plot)),2],
              y=umap_embedding[names(V(igraph.tbl_marker_plot)),1]) +
    geom_edge_diagonal(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
                    arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
    geom_node_point(aes(fill=celltype), size=7.5, stroke=0, shape=21, alpha=0.75, color = stroke_colors) +
    scale_edge_colour_manual(values=col_arrows,
                              breaks=waiver(),
                              labels=labs_plot) +
    scale_fill_manual(values=opts$celltype.colors[celltypes_present],na.value="gray") +
    theme_graph(base_family = 'Helvetica') + theme(legend.position = "none")
      
  pdf(file.path(args$outdir,sprintf("Score%s/global_network_markerTFs_signalling_sender_%s_blanco.pdf",args$min_chip_score,int)), 
      width = 6, height = 5.5)
  print(p)
  dev.off()
}



####################################################
## Plot signalling sender network on general plot ##
####################################################

print("Plot signalling network per cluster combo on general plot")

combos <- rbind(c(1,1),c(1,2),c(2,2),c(2,3),c(3,3),c(3,4),c(4,4),c(4,5),c(5,5),
                c("4+5","4+5"),
                c(1,"all"),c(2,"all"),c(3,"all"),c(4,"all"),c(5,"all"))

for (row_combo in 1:nrow(combos)){
  combo <- combos[row_combo,]
  print(paste0("Cluster ",combo[1]," vs cluster ",combo[2]))
  
  if (sum(grepl("all",combo))>0){
    interactions <- c("activating","repressive")
  } else {
    interactions <- c("activating")
  }

  for (int in interactions){
    print(int)

    if (combo[1]=="4+5"){
      TFs_marker_i <- TFs_marker[TFs_marker%in%nodes$TF[nodes$celltype%in%c(4,5)]]
    } else {
      TFs_marker_i <- TFs_marker[TFs_marker%in%nodes$TF[nodes$celltype==combo[1]]]
    }

    if (combo[2]=="all"){
      TFs_marker_j <- TFs_marker
    } else if (combo[2]=="4+5"){
      TFs_marker_j <- TFs_marker[TFs_marker%in%nodes$TF[nodes$celltype%in%c(4,5)]]
    } else {
      TFs_marker_j <- TFs_marker[TFs_marker%in%nodes$TF[nodes$celltype%in%combo[2]]]
    }

    # Create node and edge data.frames
    node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
    edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$sender%in%TFs_marker) & (GRN_coef.dt$target%in%TFs_marker)][,c("sender","target","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]

    n_pos_edges <- sum(edge_list.dt_marker$weight>0)
    n_neg_edges <- sum(edge_list.dt_marker$weight<0)

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

    # print(length(E(igraph.net_marker)$weight))
    mask_pos <- (E(igraph.tbl_marker)$weight>0) & !grepl("^Not",E(igraph.tbl_marker)$sign)
    mask_not <- (!((edge_list.dt_marker$from%in%TFs_marker_i) & (edge_list.dt_marker$to%in%TFs_marker_j))) | 
      (edge_list.dt_marker$from%in%sign2TF$from) | (edge_list.dt_marker$to%in%sign2TF$from)
    mask_neg <- (E(igraph.tbl_marker)$weight<0) & !grepl("^Not",E(igraph.tbl_marker)$sign)

    E(igraph.tbl_marker)$sign[mask_not] <- paste0("Not signalling (n=",sum(mask_not),")")
    E(igraph.tbl_marker)$sign[mask_pos & !mask_not] <- paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")")
    E(igraph.tbl_marker)$sign[mask_neg & !mask_not] <- paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")")
    E(igraph.tbl_marker)$sign <- factor(E(igraph.tbl_marker)$sign,levels=c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
                                                                            paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
                                                                            paste0("Not signalling (n=",sum(mask_not),")")))
        
    labs_plot <- c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
                    paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
                    paste0("Not signalling (n=",sum(mask_not),")"))

    print("Summary edges:")
    print(paste0("positive: ",sum(mask_pos & !mask_not)))
    print(paste0("negative: ",sum(mask_neg & !mask_not)))
    print(paste0("irrelevant: ",sum(mask_not)))

    if (int == "activating"){
      col_arrows <- c("darkred","grey","gray") 
    } else if (int == "repressive"){
      col_arrows <- c("grey","darkblue","gray") 
    } else {
      col_arrows <- c("darkred","darkblue","gray")
    }

    celltypes_markers <- rep(NA,length(names(V(igraph.net_marker))))
    for (i in 1:length(celltypes_markers)){
      if (sum(nodes$TF==names(V(igraph.net_marker))[i])>0){
        celltypes_markers[i] <- nodes$celltype[nodes$TF==names(V(igraph.net_marker))[i]]
      }
    }
    stroke_colors <- opts$celltype.colors[as.numeric(celltypes_markers)]
    celltypes_markers[names(V(igraph.net_marker))%in%sign2TF$from] <- 7
    celltypes_markers <- factor(as.character(celltypes_markers),levels=as.character(1:7))

    # only keep relevant links 
    celltypes_markers[!(names(V(igraph.net_marker))%in%unique(c(TFs_marker_i,TFs_marker_j)))] <- NA
    celltypes_markers[names(V(igraph.net_marker))%in%sign2TF$from] <- NA

    # remove activating/repressive links
    if (int == "activating"){
      set_to_remove <- unique(c(edge_list.dt_marker$from[mask_neg | mask_not],edge_list.dt_marker$to[mask_neg | mask_not]))
      set_to_remove <- set_to_remove[!(set_to_remove%in%unique(c(edge_list.dt_marker$from[mask_pos & !mask_not],edge_list.dt_marker$to[mask_pos & !mask_not])))]
      celltypes_markers[names(V(igraph.net_marker))%in%set_to_remove] <- NA
    } else if (int == "repressive"){
      set_to_remove <- unique(c(edge_list.dt_marker$from[mask_pos | mask_not],edge_list.dt_marker$to[mask_pos | mask_not]))
      set_to_remove <- set_to_remove[!(set_to_remove%in%unique(c(edge_list.dt_marker$from[mask_neg & !mask_not],edge_list.dt_marker$to[mask_neg & !mask_not])))]
      celltypes_markers[names(V(igraph.net_marker))%in%set_to_remove] <- NA
    } else {
      set_to_remove <- unique(c(edge_list.dt_marker$from[mask_not],edge_list.dt_marker$to[mask_not]))
      set_to_remove <- set_to_remove[!(set_to_remove%in%unique(c(edge_list.dt_marker$from[(mask_pos | mask_neg) & !mask_not],edge_list.dt_marker$to[(mask_pos | mask_neg) & !mask_not])))]
      celltypes_markers[names(V(igraph.net_marker))%in%set_to_remove] <- NA
    }
    print(celltypes_markers)
    stroke_colors[(celltypes_markers!=7) | is.na(celltypes_markers)] <- "#000000"

    igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(celltype=celltypes_markers)

    celltypes_present <- sort(as.numeric(unique(celltypes_markers)))
    celltypes_present <- celltypes_present[!is.na(celltypes_present)]
    print(celltypes_present)

    if (int == "activating"){
      mask_to_front <- mask_pos & !mask_not
    } else if (int == "repressive"){
      mask_to_front <- mask_neg & !mask_not
    } else {
      mask_to_front <- (mask_pos | mask_neg) & !mask_not
    }
    igraph.tbl_marker_plot = graph_from_data_frame(d=as_data_frame(igraph.tbl_marker,what="edges") %>%
      arrange(c(mask_to_front)+1),
      vertices=as_data_frame(igraph.tbl_marker,what="vertices"))
    
    set.seed(42)
    p <- ggraph(igraph.tbl_marker_plot, x=umap_embedding[names(V(igraph.tbl_marker_plot)),2],
                y=umap_embedding[names(V(igraph.tbl_marker_plot)),1]) +
      geom_edge_diagonal(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
                      arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
      geom_node_point(aes(fill=celltype), size=7.5, stroke=0, shape=21, alpha=0.75, color = stroke_colors) +
      geom_node_text(aes(label=name), size=2.25,) +
      scale_edge_colour_manual(values=col_arrows,
                                breaks=waiver(),
                                labels=labs_plot) +
      scale_fill_manual(values=opts$celltype.colors[celltypes_present],na.value="gray") +
      theme_graph(base_family = 'Helvetica') + theme(legend.position = "none")
        
    pdf(file.path(args$outdir,sprintf("Score%s/global_network_markerTFs_cluster_%s_to_%s_%s.pdf",args$min_chip_score,combo[1],combo[2],int)), 
        width = 6, height = 5.5)
    print(p)
    dev.off()

        p <- ggraph(igraph.tbl_marker_plot, x=umap_embedding[names(V(igraph.tbl_marker_plot)),2],
                y=umap_embedding[names(V(igraph.tbl_marker_plot)),1]) +
      geom_edge_diagonal(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
                      arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
      geom_node_point(aes(fill=celltype), size=7.5, stroke=0, shape=21, alpha=0.75, color = stroke_colors) +
      scale_edge_colour_manual(values=col_arrows,
                                breaks=waiver(),
                                labels=labs_plot) +
      scale_fill_manual(values=opts$celltype.colors[celltypes_present],na.value="gray") +
      theme_graph(base_family = 'Helvetica') + theme(legend.position = "none")
        
    pdf(file.path(args$outdir,sprintf("Score%s/global_network_markerTFs_cluster_%s_to_%s_%s_blanco.pdf",args$min_chip_score,combo[1],combo[2],int)), 
        width = 6, height = 5.5)
    print(p)
    dev.off()
  }
}


#################################################
## Plot expression per cluster on general plot ##
#################################################

for (j in 1:5){
  print(paste0("Plot expression per cluster ",j))
  
  # Create node and edge data.frames
  node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
  edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$sender%in%TFs_marker) & (GRN_coef.dt$target%in%TFs_marker)][,c("sender","target","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]
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
    expression_perc[i] <- (logcounts(rna_pseudobulk.sce[V(igraph.tbl_marker)$tf[i],as.character(j)])-min(logcounts(rna_pseudobulk.sce[V(igraph.tbl_marker)$tf[i],1:5])))/
      (max(logcounts(rna_pseudobulk.sce[V(igraph.tbl_marker)$tf[i],1:5]))-min(logcounts(rna_pseudobulk.sce[V(igraph.tbl_marker)$tf[i],1:5])))
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
  
  # set.seed(42)
  # p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
  #             y=umap_embedding[names(V(igraph.net_marker)),1]) +
  #   geom_edge_diagonal(edge_colour="grey", edge_alpha=0.85, edge_width=0.15,
  #                      arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
  #   geom_node_point(aes(fill=expression), size=7.5, stroke=0, shape=21, alpha=0.75) +
  #   geom_node_text(aes(label=name), size=2.25,) +
  #   scale_colour_gradientn(colors=colorRampPalette(brewer.pal(9, "YlOrRd"))(100),aesthetics = "fill") +
  #   theme_graph(base_family = 'Helvetica') + theme(legend.position = "none") #+ 
  # #annotate("text",  x=Inf, y = Inf, label = cluster_names[j], vjust=1, hjust=1)
  
  # pdf(file.path(args$outdir,sprintf("Score%s/global_network_gradient_expr_cluster%s.pdf",args$min_chip_score,j)), 
  #     width = 6, height = 5.5)

  # print(p)
  # dev.off()
  
  set.seed(42)
  p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
              y=umap_embedding[names(V(igraph.net_marker)),1]) +
    geom_edge_diagonal(edge_colour="grey", edge_alpha=0.85, edge_width=0.15,
                       arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
    geom_node_point(aes(fill=z_score), size=7.5, stroke=0, shape=21, alpha=0.75) +
    geom_node_text(aes(label=name), size=2.25,) +
    scale_colour_gradientn(colors=colorRampPalette(brewer.pal(9, "YlOrRd"))(100),aesthetics = "fill") +
    theme_graph(base_family = 'Helvetica') + theme(legend.position = "none") #+ 
  
  pdf(file.path(args$outdir,sprintf("Score%s/global_network_gradient_zscore_cluster%s.pdf",args$min_chip_score,j)), 
      width = 6, height = 5.5)
  print(p)
  dev.off()
  
    p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
              y=umap_embedding[names(V(igraph.net_marker)),1]) +
    geom_edge_diagonal(edge_colour="grey", edge_alpha=0.85, edge_width=0.15,
                       arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
    geom_node_point(aes(fill=z_score), size=7.5, stroke=0, shape=21, alpha=0.75) +
    scale_colour_gradientn(colors=colorRampPalette(brewer.pal(9, "YlOrRd"))(100),aesthetics = "fill") +
    theme_graph(base_family = 'Helvetica') + theme(legend.position = "none") #+ 
  
  pdf(file.path(args$outdir,sprintf("Score%s/global_network_gradient_zscore_cluster%s_blanco.pdf",args$min_chip_score,j)), 
      width = 6, height = 5.5)
  print(p)
  dev.off()

  # set.seed(42)
  # p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
  #             y=umap_embedding[names(V(igraph.net_marker)),1]) +
  #   geom_edge_diagonal(edge_colour="grey", edge_alpha=0.85, edge_width=0.15,
  #                      arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
  #   geom_node_point(aes(fill=z_score_red), size=7.5, stroke=0, shape=21, alpha=0.75) +
  #   geom_node_text(aes(label=name), size=2.25,) +
  #   scale_colour_gradientn(colors=colorRampPalette(brewer.pal(9, "YlOrRd"))(100),aesthetics = "fill") +
  #   theme_graph(base_family = 'Helvetica') + theme(legend.position = "none") #+ 
  # #annotate("text",  x=Inf, y = Inf, label = cluster_names[j], vjust=1, hjust=1)
  
  # pdf(file.path(args$outdir,sprintf("Score%s/global_network_gradient_zscore_red_cluster%s.pdf",args$min_chip_score,j)), 
  #     width = 6, height = 5.5)
  # print(p)
  # dev.off()
  
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

  p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
              y=umap_embedding[names(V(igraph.net_marker)),1]) +
    geom_edge_diagonal(edge_colour="grey", edge_alpha=0.85, edge_width=0.15,
                       arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
    geom_node_point(aes(fill=expression_perc), size=7.5, stroke=0, shape=21, alpha=0.75) +
    scale_colour_gradientn(colors=colorRampPalette(brewer.pal(9, "YlOrRd"))(100),aesthetics = "fill") +
    theme_graph(base_family = 'Helvetica') + theme(legend.position = "none") #+ 
    #annotate("text",  x=Inf, y = Inf, label = cluster_names[j], vjust=1, hjust=1)
  
  
  pdf(file.path(args$outdir,sprintf("Score%s/global_network_gradient_expr_perc_cluster%s_blanco.pdf",args$min_chip_score,j)), 
      width = 6, height = 5.5)
  
  print(p)
  dev.off()
}


#########################################################
## Plot expression per day of sampling on general plot ##
#########################################################

day_of_sampling <- c("d0","d1","d3","d5","d7","d10","d14","d18")
for (j in day_of_sampling){
  print(paste0("Plot expression per day of sampling ",j))
  
  # Create node and edge data.frames
  node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
  edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$sender%in%TFs_marker) & (GRN_coef.dt$target%in%TFs_marker)][,c("sender","target","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]
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
    expression_values[i] <- logcounts(rna_pseudobulk_sample.sce[V(igraph.tbl_marker)$tf[i],as.character(j)])
  }
  
  z_score_values <- z_score(expression_values)
  
  expression_perc <- rep(NA,length(V(igraph.tbl_marker)$tf))
  for (i in 1:length(V(igraph.tbl_marker))){
    expression_perc[i] <- (logcounts(rna_pseudobulk_sample.sce[V(igraph.tbl_marker)$tf[i],as.character(j)])-min(logcounts(rna_pseudobulk_sample.sce[V(igraph.tbl_marker)$tf[i],])))/
      (max(logcounts(rna_pseudobulk_sample.sce[V(igraph.tbl_marker)$tf[i],]))-min(logcounts(rna_pseudobulk_sample.sce[V(igraph.tbl_marker)$tf[i],])))
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
  
  # set.seed(42)
  # p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
  #             y=umap_embedding[names(V(igraph.net_marker)),1]) +
  #   geom_edge_diagonal(edge_colour="grey", edge_alpha=0.85, edge_width=0.15,
  #                      arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
  #   geom_node_point(aes(fill=expression), size=7.5, stroke=0, shape=21, alpha=0.75) +
  #   geom_node_text(aes(label=name), size=2.25,) +
  #   scale_colour_gradientn(colors=colorRampPalette(brewer.pal(9, "YlOrRd"))(100),aesthetics = "fill") +
  #   theme_graph(base_family = 'Helvetica') + theme(legend.position = "none") #+ 
  # #annotate("text",  x=Inf, y = Inf, label = cluster_names[j], vjust=1, hjust=1)
  
  # pdf(file.path(args$outdir,sprintf("Score%s/global_network_gradient_expr_cluster%s.pdf",args$min_chip_score,j)), 
  #     width = 6, height = 5.5)

  # print(p)
  # dev.off()
  
  set.seed(42)
  p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
              y=umap_embedding[names(V(igraph.net_marker)),1]) +
    geom_edge_diagonal(edge_colour="grey", edge_alpha=0.85, edge_width=0.15,
                       arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
    geom_node_point(aes(fill=z_score), size=7.5, stroke=0, shape=21, alpha=0.75) +
    geom_node_text(aes(label=name), size=2.25,) +
    scale_colour_gradientn(colors=colorRampPalette(brewer.pal(9, "YlOrRd"))(100),aesthetics = "fill") +
    theme_graph(base_family = 'Helvetica') + theme(legend.position = "none") #+ 
  
  pdf(file.path(args$outdir,sprintf("Score%s/global_network_gradient_zscore_sample_%s.pdf",args$min_chip_score,j)), 
      width = 6, height = 5.5)
  print(p)
  dev.off()

  p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
              y=umap_embedding[names(V(igraph.net_marker)),1]) +
    geom_edge_diagonal(edge_colour="grey", edge_alpha=0.85, edge_width=0.15,
                       arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
    geom_node_point(aes(fill=z_score), size=7.5, stroke=0, shape=21, alpha=0.75) +
    scale_colour_gradientn(colors=colorRampPalette(brewer.pal(9, "YlOrRd"))(100),aesthetics = "fill") +
    theme_graph(base_family = 'Helvetica') + theme(legend.position = "none") #+ 
  
  pdf(file.path(args$outdir,sprintf("Score%s/global_network_gradient_zscore_sample_%s_blanco.pdf",args$min_chip_score,j)), 
      width = 6, height = 5.5)
  print(p)
  dev.off()
  
  # set.seed(42)
  # p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
  #             y=umap_embedding[names(V(igraph.net_marker)),1]) +
  #   geom_edge_diagonal(edge_colour="grey", edge_alpha=0.85, edge_width=0.15,
  #                      arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
  #   geom_node_point(aes(fill=z_score_red), size=7.5, stroke=0, shape=21, alpha=0.75) +
  #   geom_node_text(aes(label=name), size=2.25,) +
  #   scale_colour_gradientn(colors=colorRampPalette(brewer.pal(9, "YlOrRd"))(100),aesthetics = "fill") +
  #   theme_graph(base_family = 'Helvetica') + theme(legend.position = "none") #+ 
  # #annotate("text",  x=Inf, y = Inf, label = cluster_names[j], vjust=1, hjust=1)
  
  # pdf(file.path(args$outdir,sprintf("Score%s/global_network_gradient_zscore_red_cluster%s.pdf",args$min_chip_score,j)), 
  #     width = 6, height = 5.5)
  # print(p)
  # dev.off()
  
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
  
  
  pdf(file.path(args$outdir,sprintf("Score%s/global_network_gradient_expr_perc_sample_%s.pdf",args$min_chip_score,j)), 
      width = 6, height = 5.5)
  
  print(p)
  dev.off()

  p <- ggraph(igraph.tbl_marker, x=umap_embedding[names(V(igraph.net_marker)),2],
            y=umap_embedding[names(V(igraph.net_marker)),1]) +
  geom_edge_diagonal(edge_colour="grey", edge_alpha=0.85, edge_width=0.15,
                      arrow=arrow(length=unit(1,'mm'),angle=15), end_cap=circle(3,'mm'), start_cap=circle(3,'mm')) +
  geom_node_point(aes(fill=expression_perc), size=7.5, stroke=0, shape=21, alpha=0.75) +
  scale_colour_gradientn(colors=colorRampPalette(brewer.pal(9, "YlOrRd"))(100),aesthetics = "fill") +
  theme_graph(base_family = 'Helvetica') + theme(legend.position = "none") #+ 
  #annotate("text",  x=Inf, y = Inf, label = cluster_names[j], vjust=1, hjust=1)
  
  
  pdf(file.path(args$outdir,sprintf("Score%s/global_network_gradient_expr_perc_sample_%s_blanco.pdf",args$min_chip_score,j)), 
      width = 6, height = 5.5)
  
  print(p)
  dev.off()
}

######################
## Completion token ##
######################

file.create(file.path(args$outdir,sprintf("Score%s/completed.txt",args$min_chip_score)))

##############################
## Code for gif (run local) ##
##############################

##############################
## Code for gif (run local) ##
##############################

# library(magick)
# # Put files to combine in a directory called gif in the working directory
# setwd("C:/Users/ljcousse/Downloads/Score0.06")

# prefixes <- c("global_network_gradient_expr_perc",
#               "global_network_gradient_zscore")

# ######################
# ## Based on cluster ##
# ######################

# clusters <- c("Naïve","Early partially primed","Late partially primed","Early primed","Late primed")

# ## White background
# ###################

# for (j in prefixes){
#   files <- list.files()[grepl(paste0("^",j,"_cluster[1-5].pdf"),list.files())]
#   print(files)
#   counter <- 1
#   for (i in 1:length(files)){
#     image_i <- image_read_pdf(files[i],density=600)
#     image_i_ann <- image_annotate(image_i,clusters[counter], size = 16, color = "black",gravity = "northeast", location = "+100+100")
#     if (i == 1){
#       images <- image_i_ann
#     } else {
#       images <- c(images,image_i_ann)
#     }
#     counter <- counter+1
#   }
  
#   animation <- image_animate(images, fps = 1)
  
#   image_write(animation,paste0("C:/Users/ljcousse/Downloads/",j,"_cluster.gif"))
# }

# ## no background, white text
# ###################

# for (j in prefixes){
#   files <- list.files()[grepl(paste0("^",j,"_cluster[1-5].pdf"),list.files())]
#   print(files)
#   counter <- 1
#   for (i in 1:length(files)){
#     image_i <- image_read_pdf(files[i],density=700)
#     image_i <- image_transparent(image_i, 'white')
#     image_i_ann <- image_annotate(image_i,clusters[counter], size = 16, color = "white",gravity = "northeast", location = "+100+100")
#     if (i == 1){
#       images <- image_i_ann
#     } else {
#       images <- c(images,image_i_ann)
#     }
#     counter <- counter+1
#   }
  
#   animation <- image_animate(images, fps = 1)
  
#   image_write(animation,paste0("C:/Users/ljcousse/Downloads/",j,"_cluster_nb.gif"))
#   gc()
# }

# #####################
# ## Based on sample ##
# #####################

# samples <- c("d0","d1","d3","d5","d7","d10","d14","d18")

# ## White background
# ###################

# for (j in prefixes){
#   files <- list.files()[grepl(paste0("^",j,"_sample_d[0-9]+.pdf"),list.files())]
#   files <- files[sort(as.numeric(gsub(".pdf","",gsub(".*_sample_d","",files))),index.return=T)$ix]
#   print(files)
#   counter <- 1
#   for (i in 1:length(files)){
#     image_i <- image_read_pdf(files[i],density=600)
#     image_i_ann <- image_annotate(image_i,samples[counter], size = 20, color = "black",gravity = "northeast", location = "+100+100")
#     if (i == 1){
#       images <- image_i_ann
#     } else {
#       images <- c(images,image_i_ann)
#     }
#     counter <- counter+1
#   }
  
#   animation <- image_animate(images, fps = 1)
  
#   image_write(animation,paste0("C:/Users/ljcousse/Downloads/",j,"_sample.gif"))
# }

# ## no background, white text
# ###################

# for (j in prefixes){
#   files <- list.files()[grepl(paste0("^",j,"_sample_d[0-9]+.pdf"),list.files())]
#   files <- files[sort(as.numeric(gsub(".pdf","",gsub(".*_sample_d","",files))),index.return=T)$ix]
#   print(files)
#   counter <- 1
#   for (i in 1:length(files)){
#     image_i <- image_read_pdf(files[i],density=700)
#     image_i <- image_transparent(image_i, 'white')
#     image_i_ann <- image_annotate(image_i,samples[counter], size = 20, color = "white",gravity = "northeast", location = "+100+100")
#     if (i == 1){
#       images <- image_i_ann
#     } else {
#       images <- c(images,image_i_ann)
#     }
#     counter <- counter+1
#   }
  
#   animation <- image_animate(images, fps = 1)
  
#   image_write(animation,paste0("C:/Users/ljcousse/Downloads/",j,"_sample_nb.gif"))
#   gc()
# }

# ################
# ## Transition ##
# ################

# prefixes <- c("global_network_markerTFs")

# ## White background
# ###################

# for (j in prefixes){
#   files <- list.files()[grepl(paste0("^",j,"_cluster_[0-9\\+]+_to_[0-9\\+]+_activating.pdf"),list.files())]
#   files <- files[-grep("to_5",files)]
#   print(files)
#   counter <- 1
#   samples <- gsub("_"," ",gsub("_act.*.pdf$","",gsub("^.*cluster_","",files)))
#   for (i in 1:length(files)){
#     image_i <- image_read_pdf(files[i],density=600)
#     image_i_ann <- image_annotate(image_i,samples[counter], size = 20, color = "black",gravity = "northeast", location = "+100+100")
#     if (i == 1){
#       images <- image_i_ann
#     } else {
#       images <- c(images,image_i_ann)
#     }
#     counter <- counter+1
#   }
  
#   animation <- image_animate(images, fps = 1)
  
#   image_write(animation,paste0("C:/Users/ljcousse/Downloads/",j,"_sample.gif"))
#   gc()
# }


# ## no background, white text
# ###################

# for (j in prefixes){
#   files <- list.files()[grepl(paste0("^",j,"_cluster_[0-9\\+]+_to_[0-9\\+]+_activating.pdf"),list.files())]
#   files <- files[-grep("to_5",files)]
#   print(files)
#   counter <- 1
#   samples <- gsub("_"," ",gsub("_act.*.pdf$","",gsub("^.*cluster_","",files)))
#   for (i in 1:length(files)){
#     image_i <- image_read_pdf(files[i],density=600)
#     image_i <- image_transparent(image_i, 'white')
#     image_i_ann <- image_annotate(image_i,samples[counter], size = 20, color = "white",gravity = "northeast", location = "+100+100")
#     if (i == 1){
#       images <- image_i_ann
#     } else {
#       images <- c(images,image_i_ann)
#     }
#     counter <- counter+1
#   }
  
#   animation <- image_animate(images, fps = 1)
  
#   image_write(animation,paste0("C:/Users/ljcousse/Downloads/",j,"_sample_nb.gif"))
#   gc()
# }