########################
##                    ##
##  Metacells_stats.R ##
##                    ##
########################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA_ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--samples',    type="character",  nargs="+", help='Samples')
p$add_argument('--metadata',    type="character",  help='Cell metadata file')
p$add_argument('--metacell',    type="character",  help='Metacell metadata file')
p$add_argument('--metacell_prefilter',    type="character",  help='Metacell metadata file pre filtering')
p$add_argument('--mofa_model',    type="character",  help='Mofa model')
p$add_argument('--seed',     type="integer",    default=42,     help='Random seed')
p$add_argument('--batch_correction',        type="logical", default=FALSE,               help='Suffix')
p$add_argument('--incl_samples',        type="character",                               help='Suffix')
p$add_argument('--outdir',          type="character",               help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

args$filter_differentiated <- FALSE

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) 

###############
## Load MOFA ##
###############

if (tools::file_ext(args$mofa_model)=="rds") {
  MOFAobject <- readRDS(args$mofa_model)
  
  # Add sample metadata
  samples_metadata(MOFAobject) <- sample_metadata
} else if (file_ext(args$mofa_model)=="hdf5") {
  MOFAobject <- load_model(args$mofa_model, load_data = F)
  
  # Add sample metadata
  cells <- as.character(unname(unlist(MOFA2::samples_names(MOFAobject))))
  sample_metadata_to_mofa <- copy(sample_metadata) %>%
    setnames("cell","sample") %>%
    .[sample%in%cells] %>% setkey(sample) %>% .[cells]
  stopifnot(all(cells==sample_metadata_to_mofa$cell))
  samples_metadata(MOFAobject) <- sample_metadata_to_mofa
}

######################
## Batch correction ##
######################


# Select factors to use 
print(sprintf("%s dimensions used.",get_dimensions(MOFAobject)[["K"]]))
factors.to.use <- 1:get_dimensions(MOFAobject)[["K"]]

# Extract factors
Z <- get_factors(MOFAobject, factors=factors.to.use)[[1]]

Z_corrected <- Z

##########
## UMAP ##
##########

# Run
set.seed(args$seed)
umap_embedding <- uwot::umap(Z_corrected, n_neighbors=200, min_dist=0.05, metric="cosine")

# Plot
to.plot <- umap_embedding %>% as.data.table %>%
  .[,sample:=rownames(Z_corrected)] %>%
  merge(MOFAobject@samples_metadata[,c("sample","batch")] %>% as.data.table)

####################
## Plot metacells ##
####################

if (file.exists(args$metacell)){
  metacell2cell.dt <- fread(args$metacell_prefilter)
  
  to.plot$type <- "cell"
  sum(to.plot$sample%in%metacell2cell.dt$metacell)
  to.plot_metacell <- to.plot[to.plot$sample%in%metacell2cell.dt$metacell,]
  to.plot_metacell$type <- "metacell"
  to.plot_metacell <- rbind(to.plot,to.plot_metacell)
  
  p <- ggplot(to.plot_metacell, aes(x=V1, y=V2, fill=type)) +
    geom_point(size=1, shape=21, stroke=0.1) +
    guides(fill = guide_legend(override.aes = list(size=2))) +
    theme_classic() +
    scale_fill_manual(values=c("gray","red")) +
    ggplot_theme_NoAxes() 
  
  output_plot(p,sprintf("%s/pdf/mofa_umap_metacell_%s_prefilter",args$outdir,args$incl_samples), width=7, height=5,UMAP=TRUE)
  
  metacell2cell.dt <- fread(args$metacell)
  
  to.plot$type <- "cell"
  sum(to.plot$sample%in%metacell2cell.dt$metacell)
  to.plot_metacell <- to.plot[to.plot$sample%in%metacell2cell.dt$metacell,]
  to.plot_metacell$type <- "metacell"
  to.plot_metacell <- rbind(to.plot,to.plot_metacell)
  
  p <- ggplot(to.plot_metacell, aes(x=V1, y=V2, fill=type)) +
    geom_point(size=1, shape=21, stroke=0.1) +
    guides(fill = guide_legend(override.aes = list(size=2))) +
    theme_classic() +
    scale_fill_manual(values=c("gray","red")) +
    ggplot_theme_NoAxes() 
  
  output_plot(p,sprintf("%s/pdf/mofa_umap_metacell_%s",args$outdir,args$incl_samples), width=7, height=5,UMAP=TRUE)
} else {
  print("First produce metacells before you plot them ;)")
}


#####################
## Check fractions ##
#####################

table_metacells_sample <- NULL
table_metacells_cluster <- NULL
if (!is.null(args$metadata)){
  metacell2cell.dt <- fread(args$metacell_prefilter)
  
  to.plot$type <- "cell"
  print(sum(to.plot$sample%in%metacell2cell.dt$metacell))
  to.plot_metacell <- to.plot[to.plot$sample%in%metacell2cell.dt$metacell,]
  to.plot_metacell$type <- "metacell"
  to.plot_metacell <- rbind(to.plot,to.plot_metacell)
  
  metadata_clusters <- fread(args$metadata)
  metadata_clusters$cluster <- factor(metadata_clusters$cluster,levels=1:6)
  
  to.plot_metacell <- to.plot_metacell %>% merge(metadata_clusters)
  to.plot_metacell$batch <- as.factor(to.plot_metacell$batch)
  
  to.plot_metacell <- to.plot_metacell %>% dplyr::rename(cell=sample) %>% merge(metacell2cell.dt)
  to.plot_metacell <- to.plot_metacell %>% dplyr::rename(sample=cell) 
  
  table_metacells_sample <- rbind(table_metacells_sample,
                                  table(to.plot_metacell$batch[to.plot_metacell$type=="cell"])[args$samples])
  table_metacells_sample <- rbind(table_metacells_sample,round(table_metacells_sample[1,]/sum(table_metacells_sample[1,])*100,digits=1)[args$samples])
  table_metacells_sample <- rbind(table_metacells_sample,
                                  table(to.plot_metacell$batch[to.plot_metacell$type=="metacell"])[args$samples])
  table_metacells_sample <- rbind(table_metacells_sample,round(table_metacells_sample[3,]/sum(table_metacells_sample[3,])*100,digits=1)[args$samples])
  
  table_metacells_cluster <- rbind(table_metacells_cluster, 
                                   table(to.plot_metacell$cluster[to.plot_metacell$type=="cell"]))
  table_metacells_cluster <- rbind(table_metacells_cluster,round(table_metacells_cluster[1,]/sum(table_metacells_cluster[1,])*100,digits=1))
  table_metacells_cluster <- rbind(table_metacells_cluster, 
                                   table(to.plot_metacell$cluster[to.plot_metacell$type=="metacell"]))
  table_metacells_cluster <- rbind(table_metacells_cluster,round(table_metacells_cluster[3,]/sum(table_metacells_cluster[3,])*100,digits=1))
  
  metacell2cell.dt <- fread(args$metacell)
  
  print(sum(to.plot_metacell$metacell%in%metacell2cell.dt$metacell))
  to.plot_metacell <- to.plot_metacell[to.plot_metacell$metacell%in%metacell2cell.dt$metacell,]
  
  table_metacells_sample <- rbind(table_metacells_sample,
                                  table(to.plot_metacell$batch[to.plot_metacell$type=="cell"])[args$samples])
  table_metacells_sample <- rbind(table_metacells_sample,round(table_metacells_sample[5,]/sum(table_metacells_sample[5,])*100,digits=1)[args$samples])
  table_metacells_sample <- rbind(table_metacells_sample,
                                  table(to.plot_metacell$batch[to.plot_metacell$type=="metacell"])[args$samples])
  table_metacells_sample <- rbind(table_metacells_sample,round(table_metacells_sample[7,]/sum(table_metacells_sample[7,])*100,digits=1)[args$samples])
  
  table_metacells_cluster <- rbind(table_metacells_cluster, 
                                   table(to.plot_metacell$cluster[to.plot_metacell$type=="cell"]))
  table_metacells_cluster <- rbind(table_metacells_cluster,round(table_metacells_cluster[5,]/sum(table_metacells_cluster[5,])*100,digits=1))
  table_metacells_cluster <- rbind(table_metacells_cluster, 
                                   table(to.plot_metacell$cluster[to.plot_metacell$type=="metacell"]))
  table_metacells_cluster <- rbind(table_metacells_cluster,round(table_metacells_cluster[7,]/sum(table_metacells_cluster[7,])*100,digits=1))
  
}

write.table(table_metacells_sample,sprintf("%s/table_metacells_per_sample_%s.txt",args$outdir,args$incl_samples), 
            sep="\t", col.names = F, row.names = F, quote = F)
write.table(table_metacells_cluster,sprintf("%s/table_metacells_per_cluster_%s.txt",args$outdir,args$incl_samples), 
            sep="\t", col.names = F, row.names = F, quote = F)


###########################
## Load metacell results ##
###########################

cell2metacell.dt <- args$metacell_prefilter %>% map(~ fread(.)) %>% rbindlist
cell2metacell.dt <- cell2metacell.dt[cell2metacell.dt$metacell%in%metacell2cell.dt$metacell,]

print(sprintf("Number of metacells: %s", length(unique(cell2metacell.dt$metacell))))

print(summary(as.numeric(table(cell2metacell.dt$metacell))))
pdf(paste0(args$outdir,"/distribution_cells_per_metacell_",args$incl_samples,".pdf"))
hist(as.numeric(table(cell2metacell.dt$metacell)),breaks=50,main="",xlab="Amount cells per metacell")
dev.off()
