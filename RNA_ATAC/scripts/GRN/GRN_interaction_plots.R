###############################
##                           ##
##  GRN_interaction_plots.R  ##
##                           ##
###############################

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--interaction_table',    type="character",    help='Interaction table for GRN')
p$add_argument('--outdir',   type="character",    help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

# # START TEST
# ## Settings for alternative calculation of coefficients!
# args <- list()
# args$interaction_table <-  "/data/homes/louisc/Project_Babraham/RNA_ATAC/gene_regulatory_network/metacells/trajectories/N2P/alt3/marker_TF/fig/Score0.06_Coef0.25/table_interactions.txt"
# args$outdir <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/gene_regulatory_network/metacells/trajectories/N2P/alt3/marker_TF/Interaction_plots/Score0.06_Coef0.25"
# # END TEST

#####################
## Define settings ##
#####################

if(!dir.exists(args$outdir)){
  dir.create(paste(unlist(strsplit(args$outdir,"/"))[-length(unlist(strsplit(args$outdir,"/")))],collapse="/"), showWarnings = F)
  dir.create(args$outdir)
}

####################
## Color settings ##
####################

cols <- rev(brewer.pal(9, "RdBu"))
cols1 <- rev(brewer.pal(9, "RdYlBu"))
cols2 <- c("deeppink","black","turquoise1")
cols3 <- (brewer.pal(9, "YlOrRd"))

modes <- c("all","main")
for (mode in modes){
  print(paste0("Running ",mode," mode!"))
  
  ##################################
  ## Read table with interactions ##
  ##################################
  
  if (mode=="all"){
    fread(args$interaction_table) %>%
      as_tibble() %>%
      dplyr::select("from", "to", "beta", "from_cluster_all", "to_cluster_all") %>%
      na.omit() -> interactions
    print(dim(interactions))
    print(head(interactions))
    gaps_mode <- c(4,5)
  } else {
    fread(args$interaction_table) %>%
      as_tibble() %>%
      dplyr::select("from", "to", "beta", "from_cluster", "to_cluster") %>%
      dplyr::rename(from_cluster_all = from_cluster, to_cluster_all = to_cluster) %>%
      na.omit() -> interactions
    print(dim(interactions))
    print(head(interactions))
    gaps_mode <- c(2,3)
  }
  
  ##################
  ## Make TF list ##
  ##################
  
  ### List of receivers
  
  interactions %>%
    dplyr::select(to, to_cluster_all) %>%
    dplyr::distinct(to, .keep_all = T) %>%
    dplyr::rename(TF = to, cluster = to_cluster_all) -> TF.list.1
  
  ### List of senders
  
  interactions %>%
    dplyr::select(from, from_cluster_all) %>%
    dplyr::distinct(from, .keep_all=T) %>%
    dplyr::rename(TF = from, cluster = from_cluster_all) -> TF.list.2
  
  ### Full list
  
  TF.list.1 %>%
    full_join(TF.list.2) -> TF.list
  
  TF.list %>%
    group_by(cluster) %>%
    tally() %>%
    dplyr::rename(cluster.size = n) -> cluster.sizes
  
  print(head(TF.list))
  print(cluster.sizes)
  
  interactions %>%
    arrange(from_cluster_all) %>%
    pull(from_cluster_all) %>%
    unique -> TF.groups
  
  print(TF.groups)
  
  ##################################
  ## Count activating connections ##
  ##################################
  
  count.activation <- function(x){
    
    interactions %>%
      filter(beta >0) %>%
      filter(from_cluster_all == x) %>%
      group_by(to_cluster_all) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::rename(to_cluster = to_cluster_all) %>%
      left_join(cluster.sizes %>% dplyr::rename(to_cluster = cluster)) %>%
      dplyr::rename(receiver.cluster.size = cluster.size) %>%
      mutate(sender.cluster.size = (cluster.sizes %>% filter(cluster == x) %>% pull(cluster.size))) %>%
      mutate(freq = n / sum(n)) %>%
      mutate(index.1 = freq / receiver.cluster.size) %>%
      mutate(index.2 = n / (receiver.cluster.size * sender.cluster.size)) %>%
      mutate(from_cluster = x) %>%
      dplyr::relocate(from_cluster, from_cluster, to_cluster, receiver.cluster.size, sender.cluster.size, everything()) -> my.result
    
    return(my.result)
  }
  
  res_activation <- lapply(TF.groups, count.activation)  %>%
    do.call("rbind",.) %>%
    as_tibble() -> active.connection
  
  print(data.frame(active.connection))
  
  ##################################
  ## Count repressing connections ##
  ##################################
  
  count.repression <- function(x){
    
    interactions %>%
      filter(beta <0) %>%
      filter(from_cluster_all == x) %>%
      group_by(to_cluster_all) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::rename(to_cluster = to_cluster_all) %>%
      left_join(cluster.sizes %>% dplyr::rename(to_cluster = cluster)) %>%
      dplyr::rename(receiver.cluster.size = cluster.size) %>%
      mutate(sender.cluster.size = (cluster.sizes %>% filter(cluster == x) %>% pull(cluster.size))) %>%
      mutate(freq = n / sum(n)) %>%
      mutate(index.1 = freq / receiver.cluster.size) %>%
      mutate(index.2 = n / (receiver.cluster.size * sender.cluster.size)) %>%
      mutate(from_cluster = x) %>%
      dplyr::relocate(from_cluster, from_cluster, to_cluster, receiver.cluster.size, sender.cluster.size, everything()) -> my.result
    
    return(my.result)
  }
  
  lapply(TF.groups, count.repression) %>%
    do.call("rbind",.) %>%
    as_tibble() -> repressive.connection
  
  print(data.frame(repressive.connection))
  
  ##################
  ## Heatmaps all ##
  ##################
  
  ### Number of connections
  
  active.connection %>%
    mutate(n = round(n), digits = 0) %>%
    select(c(to_cluster, from_cluster, n)) %>%
    pivot_wider(names_from = "from_cluster", values_from = "n") %>%
    mutate_all(~replace(., is.na(.), 0)) %>%
    arrange(to_cluster) %>%
    column_to_rownames(var = "to_cluster") %>%
    pheatmap(
      show_rownames = T,
      show_colnames = T,
      cluster_rows = F,
      cluster_cols = F,
      gaps_col = gaps_mode,
      gaps_row = gaps_mode,
      display_numbers = T,
      number_format = "%.0f",
      number_color = "white",
      fontsize_number = 20,
      color = colorRampPalette(brewer.pal(6,"OrRd"))(100),breaks=seq(from=0, to=20, length.out=101),
      filename = paste0(args$outdir,"/transition.activating.connections.number.",mode,".clusters.pdf"), height = 6, width = 6.5
    )
  
  repressive.connection %>%
    mutate(n = round(n), digits = 0) %>%
    select(c(to_cluster, from_cluster, n)) %>%
    pivot_wider(names_from = "from_cluster", values_from = "n") %>%
    mutate_all(~replace(., is.na(.), 0)) %>%
    arrange(to_cluster) %>%
    column_to_rownames(var = "to_cluster") %>%
    pheatmap(
      show_rownames = T,
      show_colnames = T,
      cluster_rows = F,
      cluster_cols = F,
      gaps_col = gaps_mode,
      gaps_row = gaps_mode,
      display_numbers = T,
      number_format = "%.0f",
      number_color = "white",
      fontsize_number = 20,
      color = colorRampPalette(brewer.pal(6,"Blues"))(100),breaks=seq(from=0, to=20, length.out=101),
      filename = paste0(args$outdir,"/transition.repressive.connections.number.",mode,".clusters.pdf"), height = 6, width = 6.5
    )
  
  
  ### Frequency of connections
  
  active.connection %>%
    mutate(freq = freq *100) %>%
    select(c(to_cluster, from_cluster, freq)) %>%
    pivot_wider(names_from = "from_cluster", values_from = "freq") %>%
    mutate_all(~replace(., is.na(.), 0)) %>%
    arrange(to_cluster) %>%
    column_to_rownames(var = "to_cluster") %>%
    pheatmap(
      show_rownames = T,
      show_colnames = T,
      cluster_rows = F,
      cluster_cols = F,
      gaps_col = gaps_mode,
      gaps_row = gaps_mode,
      display_numbers = T,
      number_color = "white",
      fontsize_number = 20,
      color = colorRampPalette(brewer.pal(6,"OrRd"))(100),breaks=seq(from=0, to=20, length.out=101),
      filename = paste0(args$outdir,"/transition.activating.connections.freq.",mode,".clusters.pdf"), height = 6, width = 6.5
    )
  
  repressive.connection %>%
    mutate(freq = freq *100) %>%
    select(c(to_cluster, from_cluster, freq)) %>%
    pivot_wider(names_from = "from_cluster", values_from = "freq") %>%
    mutate_all(~replace(., is.na(.), 0)) %>%
    arrange(to_cluster) %>%
    column_to_rownames(var = "to_cluster") %>%
    pheatmap(
      show_rownames = T,
      show_colnames = T,
      cluster_rows = F,
      cluster_cols = F,
      gaps_col = gaps_mode,
      gaps_row = gaps_mode,
      display_numbers = T,
      number_color = "white",
      fontsize_number = 20,
      color = colorRampPalette(brewer.pal(6,"Blues"))(100),breaks=seq(from=0, to=20, length.out=101),
      filename = paste0(args$outdir,"/transition.repressive.connections.freq.",mode,".clusters.pdf"), height = 6, width = 6.5
    )
  
  
  ### Index 1: number of connections normalised to the number of outgoing connections and receiver cluster size - THIS SEEMS THE MOST RELIABLE
  
  active.connection %>%
    mutate(index.1 = index.1 *100) %>%
    select(c(to_cluster, from_cluster, index.1)) %>%
    pivot_wider(names_from = "from_cluster", values_from = "index.1") %>%
    mutate_all(~replace(., is.na(.), 0)) %>%
    arrange(to_cluster) %>%
    column_to_rownames(var = "to_cluster") %>%
    pheatmap(
      show_rownames = T,
      show_colnames = T,
      cluster_rows = F,
      cluster_cols = F,
      gaps_col = gaps_mode,
      gaps_row = gaps_mode,
      display_numbers = T,
      number_color = "white",
      fontsize_number = 20,
      color = colorRampPalette(brewer.pal(6,"OrRd"))(100),breaks=seq(from=0, to=3, length.out=101),
      filename = paste0(args$outdir,"/transition.activating.connections.index.1.",mode,".clusters.pdf"), height = 6, width = 6.5
    )
  
  repressive.connection %>%
    mutate(index.1 = index.1 *100) %>%
    select(c(to_cluster, from_cluster, index.1)) %>%
    pivot_wider(names_from = "from_cluster", values_from = "index.1") %>%
    mutate_all(~replace(., is.na(.), 0)) %>%
    arrange(to_cluster) %>%
    column_to_rownames(var = "to_cluster") %>%
    pheatmap(
      show_rownames = T,
      show_colnames = T,
      cluster_rows = F,
      cluster_cols = F,
      gaps_col = gaps_mode,
      gaps_row = gaps_mode,
      display_numbers = T,
      number_color = "white",
      fontsize_number = 20,
      color = colorRampPalette(brewer.pal(6,"Blues"))(100),breaks=seq(from=0, to=3, length.out=101),
      filename = paste0(args$outdir,"/transition.repressive.connections.index.1.",mode,".clusters.pdf"), height = 6, width = 6.5
    )
  
  
  ### Index 2: number of connections normalised to the number of outgoing connections and receiver cluster size. This index might be skewed because of too many receivers in cluster 3.
  
  active.connection %>%
    mutate(index.2 = index.2 *100) %>%
    select(c(to_cluster, from_cluster, index.2)) %>%
    pivot_wider(names_from = "from_cluster", values_from = "index.2") %>%
    mutate_all(~replace(., is.na(.), 0)) %>%
    arrange(to_cluster) %>%
    column_to_rownames(var = "to_cluster") %>%
    pheatmap(
      show_rownames = T,
      show_colnames = T,
      cluster_rows = F,
      cluster_cols = F,
      gaps_col = gaps_mode,
      gaps_row = gaps_mode,
      display_numbers = T,
      number_color = "white",
      fontsize_number = 20,
      color = colorRampPalette(brewer.pal(6,"OrRd"))(100),breaks=seq(from=0, to=10, length.out=101),
      filename = paste0(args$outdir,"/transition.activating.connections.index.2.",mode,".clusters.pdf"), height = 6, width = 6.5
    )
  
  repressive.connection %>%
    mutate(index.2 = index.2 *100) %>%
    select(c(to_cluster, from_cluster, index.2)) %>%
    pivot_wider(names_from = "from_cluster", values_from = "index.2") %>%
    mutate_all(~replace(., is.na(.), 0)) %>%
    arrange(to_cluster) %>%
    column_to_rownames(var = "to_cluster") %>%
    pheatmap(
      show_rownames = T,
      show_colnames = T,
      cluster_rows = F,
      cluster_cols = F,
      gaps_col = gaps_mode,
      gaps_row = gaps_mode,
      display_numbers = T,
      number_color = "white",
      fontsize_number = 20,
      color = colorRampPalette(brewer.pal(6,"Blues"))(100),breaks=seq(from=0, to=10, length.out=101),
      filename = paste0(args$outdir,"/transition.repressive.connections.index.2.",mode,".clusters.pdf"), height = 6, width = 6.5
    )
}

######################
## Completion token ##
######################

print(warnings())

file.create(file.path(args$outdir,"completed.txt"))


