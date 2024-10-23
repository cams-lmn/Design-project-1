##########################################
##                                      ##
##  Plot_GRN_metacells_trajectories.R  ##
##                                      ##
##########################################

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
p$add_argument('--max_pval',  type="double",            default=0.25,      help='Maximal pvalue for linear regression')
p$add_argument('--min_chip_score',  type="double",            default=0.25,      help='Minimal chip score')
p$add_argument('--TFs_to_plot',       type="character",                help='File containing TF of interest to plot')
p$add_argument('--force',       type="logical",  default=F         ,      help='Should replotting be forced?')
p$add_argument('--outdir',       type="character",                help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

# # args$rna_cells.sce <- io$rna.sce
# args <- list()
# io$rna_metacells.sce <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/trajectories/N2P/N2P_SingleCellExperiment.rds"
# io$rna_pseudobulk.sce <- file.path(io$basedir,"results/rna/pseudobulk/celltype/SingleCellExperiment_pseudobulk.rds")
# io$trajectory_file <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/metacell_trajectory.txt.gz")
# io$trajectory <- "N2P"
# io$metacell_metadata <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/metacells_metadata.txt.gz")
# io$grn_coef <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/gene_regulatory_networks/metacells/trajectories/N2P/global_chip_GRN_coef.txt.gz"
# io$outdir <-  "/data/homes/louisc/Project_Babraham/RNA_ATAC/gene_regulatory_networks/metacells/trajectories/N2P/fig"
# 
# # Options
# opts <- list()
# opts$min_coef <- 0.25

#####################
## Define settings ##
#####################

args$trajectory <- sprintf("/data/homes/louisc/Project_Babraham/RNA_ATAC/trajectories/%s/%s_trajectory.txt.gz",
                           args$trajectory_name,args$trajectory_name)

if(!dir.exists(file.path(args$outdir))){
  dir.create(args$outdir, showWarnings = F)
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

# #####################
# ## Load trajectory ##
# #####################
# 
# metacell_metadata.dt <- fread(args$metacell_metadata)
# 
# trajectory.dt <- fread(args$trajectory) %>%
#   setnames(c("metacell","FA1","FA2")) %>%
#   merge(metacell_metadata.dt[,c("metacell","celltype")])


# ##############################
# ## Load marker gene and TFs ##
# ##############################
# 
# # TFs <- fread(paste0(io$basedir,"/results/rna_atac/gene_regulatory_networks/TFs.txt"))[[1]]
# marker_TFs_all.dt <- fread(args$rna.atlas.marker_TFs.all) %>% 
#   .[celltype%in%celltypes.to.plot] %>%
#   .[,gene:=toupper(gene)]
# 
# marker_TFs_filt.dt <- marker_TFs_all.dt %>% 
#   .[score>=opts$min_tf_score]
# 
# print(unique(marker_TFs_filt.dt$gene))
# 
# marker_TFs_all.dt <- marker_TFs_all.dt[celltype!="Caudal_Mesoderm"]

##############################
## Load RNA expression data ##
##############################

# SingleCellExperiment at cellular resolution
# rna_cells.sce <- readRDS(args$rna_cells.sce)[,trajectory.dt$cell]

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

GRN_coef.dt <- fread(args$grn_coef) %>% 
  .[,gene:=toupper(gene)] %>% # .[gene%in%unique(marker_TFs_filt.dt$gene) & tf%in%unique(marker_TFs_filt.dt$gene)] %>%
  .[pvalue<args$max_pval & abs(beta)>=args$min_coef]

# Consider only positive links
# GRN_coef.dt <- GRN_coef.dt[beta>=0]

# GRN_coef.dt[tf=="SOX2" | gene=="T"] %>% View
# GRN_coef.dt[tf=="SOX2"] %>% View

##########################
## Filter TFs and genes ##
##########################

# Filter TFs that have few connections
# TFs <- intersect(GRN_coef.dt[,.N,by="tf"][N>=3,tf], GRN_coef.dt[,.N,by="gene"][N>=2,gene])
TFs <- union(GRN_coef.dt[,.N,by="tf"][N>=3,tf], GRN_coef.dt[,.N,by="gene"][N>=3,gene])
print("TFs")
head(TFs)
print("GRN_coef.dt")
head(GRN_coef.dt)
GRN_coef.dt <- GRN_coef.dt[tf%in%TFs & gene%in%TFs]

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

# Scale
# rna_tf_pseudobulk_scaled.mtx <- apply(rna_tf_pseudobulk.mtx,1,minmax.normalisation) %>% .[celltypes.to.plot,] %>% t 
rna_tf_pseudobulk_scaled.mtx <- apply(rna_tf_pseudobulk.mtx[,celltypes.to.plot],1,minmax.normalisation) %>% t 
rna_tf_metacells_scaled.mtx <- apply(rna_tf_metacells.mtx,1,minmax.normalisation) %>% t 

TFs_marker <- unique(c(GRN_coef.dt$tf,GRN_coef.dt$gene))
TFs_marker <- TFs_marker[TFs_marker%in%markers_TF_undupl$TF]
length(TFs_marker)

###############
## Create_KO ##
###############

edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$tf%in%TFs_marker) & 
                                     (GRN_coef.dt$gene%in%TFs_marker) &
                                     (GRN_coef.dt$beta>0)][,c("tf","gene","beta")] %>% 
  setnames(c("from","to","weight")) %>% .[!from==to]
print(edge_list.dt_marker[(edge_list.dt_marker$from=="HOXB2") | (edge_list.dt_marker$to=="HOXB2")])

##############################################
## Plot network per cluster on general plot ##
##############################################

col_arrows <- c("darkred","darkblue","gray")
res_ko <- data.frame(matrix(rep(0,length(unique(markers_TF_undupl$TF))*5),ncol=5))
rownames(res_ko) <- unique(markers_TF_undupl$TF)
counter <- 1
for (gene in unique(markers_TF_undupl$TF)){
  # gene <- "ESRRB"
  cluster_gene <-  unique(markers_TF_undupl$celltype[markers_TF_undupl$TF==gene])
  print(paste0("Gene ",gene," in cluster ",cluster_gene))
  genes_to_remove <- NULL
  passed_ko <- FALSE
  
  print("Determine KO")
  for (j in 1:5){
    print(paste0("Cluster",j))
    
    # Create node and edge data.frames
    TFs_marker_j <- TFs_marker[TFs_marker%in%markers_TF_undupl$TF[markers_TF_undupl$celltype==j]]
    TFs_marker_atleast_j <- TFs_marker[TFs_marker%in%markers_TF_undupl$TF[markers_TF_undupl$celltype>=j]]
    node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
    edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$tf%in%TFs_marker_j) & 
                                         (GRN_coef.dt$gene%in%TFs_marker_atleast_j) &
                                         (GRN_coef.dt$beta>0)][,c("tf","gene","beta")] %>% 
      setnames(c("from","to","weight")) %>% .[!from==to]

    # Create ko
    edge_list.dt_marker_ko <- edge_list.dt_marker
    if ((j == 1) | length(genes_to_remove)==0){
      genes_to_remove <- gene
    } 

    while(sum(genes_to_remove%in%edge_list.dt_marker_ko$from)>0){
      edge_list.dt_marker_ko <- edge_list.dt_marker[!(edge_list.dt_marker$from%in%genes_to_remove)]

      genes_to_remove_new <- edge_list.dt_marker$to[edge_list.dt_marker$from%in%genes_to_remove]
      genes_to_remove_new <- genes_to_remove_new[!(genes_to_remove_new%in%edge_list.dt_marker_ko$to)]
  
      genes_to_remove <- sort(unique(c(genes_to_remove,genes_to_remove_new)))
    }
    print(genes_to_remove)
    
    if (j<cluster_gene){
      res_ko[counter,j] <- length(genes_to_remove) - 1
    } else {
      res_ko[counter,j] <- length(genes_to_remove)
    }
  }
  
  for (j in 1:5){
    print(paste0("Cluster",j))
    
    # Create node and edge data.frames
    TFs_marker_j <- TFs_marker[TFs_marker%in%markers_TF_undupl$TF[markers_TF_undupl$celltype==j]]
    TFs_marker_atleast_j <- TFs_marker[TFs_marker%in%markers_TF_undupl$TF[markers_TF_undupl$celltype>=j]]
    node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
    edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$tf%in%TFs_marker_j) & 
                                         (GRN_coef.dt$gene%in%TFs_marker_atleast_j) &
                                         (GRN_coef.dt$beta>0)][,c("tf","gene","beta")] %>% 
      setnames(c("from","to","weight")) %>% .[!from==to]
    
    edge_list.dt_marker_ko <- edge_list.dt_marker
    while(sum(genes_to_remove%in%edge_list.dt_marker_ko$from)>0){
      edge_list.dt_marker_ko <- edge_list.dt_marker[!(edge_list.dt_marker$from%in%genes_to_remove)]

      genes_to_remove_new <- edge_list.dt_marker$to[edge_list.dt_marker$from%in%genes_to_remove]
      genes_to_remove_new <- genes_to_remove_new[!(genes_to_remove_new%in%edge_list.dt_marker_ko$to)]
      
      genes_to_remove <- sort(unique(c(genes_to_remove,genes_to_remove_new)))
    }
    print(genes_to_remove)
    
    if (res_ko[counter,5]<length(genes_to_remove)){
      res_ko[counter,j] <- res_ko[counter,j] + (length(genes_to_remove)-res_ko[counter,5])
    }
  }
  counter <- counter + 1
  
  print(paste0("Plot network for gene ",gene))
  
  # Create node and edge data.frames
  TFs_marker <- TFs_marker[TFs_marker%in%markers_TF_undupl$TF]
  node_list.dt_marker <- data.table(node_id=1:length(TFs_marker), node_name=TFs_marker)
  edge_list.dt_marker <- GRN_coef.dt[(GRN_coef.dt$tf%in%TFs_marker) & 
                                       (GRN_coef.dt$gene%in%TFs_marker)][,c("tf","gene","beta")] %>%
    setnames(c("from","to","weight")) %>% .[!from==to]
  
  print(genes_to_remove)
  
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
  
  mask_pos <- (E(igraph.tbl_marker)$weight>0)
  mask_neg <- (E(igraph.tbl_marker)$weight<0)
  mask_not <- (edge_list.dt_marker$from%in%c(genes_to_remove,gene)) |
    (edge_list.dt_marker$to%in%c(genes_to_remove,gene))
  # print(sum(mask_pos & mask_not))
  # print(sum(mask_neg & mask_not))
  E(igraph.tbl_marker)$sign[mask_not] <- paste0("Knocked down (n=",sum(mask_not),")")
  E(igraph.tbl_marker)$sign[mask_pos & !mask_not] <- paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")")
  E(igraph.tbl_marker)$sign[mask_neg & !mask_not] <- paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")")
  E(igraph.tbl_marker)$sign <- factor(E(igraph.tbl_marker)$sign,levels=c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
                                                                         paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
                                                                         paste0("Knocked down (n=",sum(mask_not),")")))
  print(table(E(igraph.tbl_marker)$sign))
  wh_arrows <- as.numeric(names(table(E(igraph.tbl_marker)$sign)))
  labs_plot <- c(paste0("Positive (n=",sum(mask_pos)-sum(mask_pos & mask_not),")"),
                 paste0("Negative (n=",sum(mask_neg)-sum(mask_neg & mask_not),")"),
                 paste0("Knocked down (n=",sum(mask_not),")"))
  
  
  celltypes_markers <- rep(NA,length(names(V(igraph.net_marker))))
  for (i in 1:length(celltypes_markers)){
    if (sum(markers_TF_undupl$TF==names(V(igraph.net_marker))[i])>0){
      celltypes_markers[i] <- markers_TF_undupl$celltype[markers_TF_undupl$TF==names(V(igraph.net_marker))[i]]
    }
  }
  
  celltypes_markers[names(V(igraph.net_marker))%in%genes_to_remove] <- NA
  
  igraph.tbl_marker <- igraph.tbl_marker %>% activate(nodes) %>% mutate(celltype=celltypes_markers)
  
  celltypes_present <- sort(as.numeric(unique(celltypes_markers)))
  celltypes_present <- celltypes_present[!is.na(celltypes_present)]
  print(celltypes_present)

  set.seed(42)
  p <- ggraph(igraph.tbl_marker, 'stress') +
    geom_edge_link(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.15,
                   arrow=arrow(length=unit(1.5,'mm')), end_cap=circle(1,'mm')) +
    geom_node_point(aes(fill=celltype), size=3, stroke=0, shape=21, alpha=0.75) +
    geom_node_text(aes(label=name), size=1) +
    scale_edge_colour_manual(values=col_arrows[wh_arrows],
                             breaks=waiver(),
                             labels=labs_plot[wh_arrows]) +
    scale_fill_manual(values=opts$celltype.colors[celltypes_present],na.value="gray") +
    theme_graph(base_family = 'Helvetica')
  
  pdf(file.path(args$outdir,paste0("/",gene,"_ko_markerTFs.pdf")),
      width = 5.5, height = 5.75)
  print(p)
  dev.off()
}

print(res_ko)

colnames(res_ko) <- paste0(rep("cluster",5),1:5)
write.table(res_ko,file=file.path(args$outdir,"table_ko.txt"),col.names = T,row.names = T,sep="\t",quote=F)
######################
## Completion token ##
######################

file.create(file.path(args$outdir,"completed.txt"))

