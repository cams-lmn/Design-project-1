################################################
##                                            ##
##   Plot_exprs_vs_access_BrowserTrack_genes  ##
##                                            ##
################################################

source("/data/homes/louisc/Project_Babraham/RNA_ATAC/scripts/Settings.R")

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',    type="character",    help='ArchR directory')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--genome',           type="character", default="hg38",      help='Genome')
p$add_argument('--group_by',     type="character",    help='Metadata column to group by')
p$add_argument('--norm_method',     type="character", default="ReadsInTSS",    help='Normalisation method')
p$add_argument('--min_cells',     type="integer", default=100,    help='Minimum number of cells per celltype')
p$add_argument('--tile_size',     type="integer", default=100,    help='Tile size')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
p$add_argument('--outdir',    type="character",    help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

# ## START TEST ##
# args <- list()
# args$group_by <- "cluster"
# args$norm_method <- c("ReadsInTSS")
# args$tile_size <- 100
# args$min_cells <- 100
# args$genome <- "hg38"
# args$threads <- 4
# args$archr_directory <- "/data/homes/louisc/Project_Babraham/ATAC/archR/"
# ## END TEST ##

# ## START TEST ##
args <- list()
args$outdir <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/rna_vs_acc/BrowserTrack/"
args$archr_directory <- "/data/homes/louisc/Project_Babraham/ATAC/archR/"
args$metadata <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/clustering/PeakMatrix/sample_metadata_nodiff_after_clustering.txt.gz"
args$tileSize <- 100
args$tf2gene_virtual_chip <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/gene_regulatory_network/metacells/CISBP/TF2gene_after_virtual_chip.txt.gz"
args$max_distance <- 5e4
args$sce <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/SingleCellExperiment_pseudobulk.rds"
args$genome <- "hg38"
args$threads <- 1
args$grn_coef <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/gene_regulatory_network/metacells/trajectories/N2P/TFandSign/alt7/merged/global_chip_GRN_sign_alt7_coef_score0.06_merged_cleaned.txt.gz"

# args$sce <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/RNA/pseudobulk_with_replicates.rds"
# args$atac_peak_matrix <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/PeakMatrix/pseudobulk_with_replicates.rds"
# args$motif_annotation <- "CISBP"
# args$motifmatcher <- sprintf("/data/homes/louisc/Project_Babraham/ATAC/archR/Annotations/%s-Scores.rds",args$motif_annotation)
# args$motif2gene <- sprintf("/data/homes/louisc/Project_Babraham/ATAC/archR/Annotations/%s_motif2gene.txt.gz",args$motif_annotation)
# args$test <- FALSE
# args$force_rerun <- TRUE
# args$remove_clusters <- NULL
# ## END TEST ##

#####################
## Define settings ##
#####################

if(!dir.exists(args$outdir)){dir.create(args$outdir)}

########################
## Load cell metadata ##
########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE]
colnames(sample_metadata)[1] <- "cell"

stopifnot(args$group_by%in%colnames(sample_metadata))
sample_metadata <- sample_metadata[!is.na(sample_metadata[[args$group_by]])]

# Filter celltypes by minimum number of cells
sample_metadata <- sample_metadata[,N:=.N,by=c(args$group_by)] %>% .[N>=args$min_cells] %>% .[,N:=NULL]

print(table(sample_metadata[[args$group_by]]))
print(head(sample_metadata))

########################
## Load ArchR project ##
########################

# source(here::here("atac/archR/load_archR_project.R"))

setwd(args$archr_directory)

addArchRGenome(args$genome)
addArchRThreads(threads = args$threads)

ArchRProject <- loadArchRProject(args$archr_directory)[sample_metadata$cell]

############################################
## Add cluster annotation to ArchRProject ##
############################################

sum(sample_metadata$cell==rownames(getCellColData(ArchRProject)))

ArchRProject <- addCellColData(
    ArchRProject,
    data = sample_metadata$cluster,
    name = "Clusters",
    cells = sample_metadata$cell,
    force = TRUE
  )

###########################################
## Load RNA expression data - pseudobulk ##
###########################################

sce <- readRDS(args$sce)

####################################################
## Load TF2gene links based on in silico ChIP-seq ##
####################################################

tf2gene_chip.dt <- fread(args$tf2gene_virtual_chip) %>%
  .[dist<=args$max_distance]  

gr_insilico_chip_seq <- GRanges(seqnames=gsub(":.*$","",tf2gene_chip.dt$peak),
                                      ranges=IRanges(start=as.numeric(gsub("-.*$","",gsub("^.*:","",tf2gene_chip.dt$peak))),
                                                     end=as.numeric(gsub("^.*-","",gsub("^.*:","",tf2gene_chip.dt$peak)))),
                                      chip_score=tf2gene_chip.dt$chip_score,
                                      tf=tf2gene_chip.dt$tf,
                                      gene=tf2gene_chip.dt$gene)
gr_insilico_chip_seq <- gr_insilico_chip_seq[!(gr_insilico_chip_seq$gene=="")]

#########################
## Load GRN ##
#########################

GRN_coef.dt <- fread(args$grn_coef)


####################################
## Select interesting connections ##
####################################

top <- 100

list_interesting <- data.frame(tf2gene_chip.dt)
print(dim(list_interesting))
list_interesting <- list_interesting[sort(list_interesting$chip_score,decreasing=T,index.return=T)$ix,]
list_interesting <- list_interesting[!duplicated(paste(list_interesting$tf,list_interesting$gene,sep="_")),]
print(dim(list_interesting))

GRN_coef.dt_strong_repr <-  GRN_coef.dt[GRN_coef.dt$beta < -0.75,]
GRN_coef.dt_strong_act <-  GRN_coef.dt[GRN_coef.dt$beta > 0.75,]

list_interesting_repr <- list_interesting[paste(list_interesting$tf,list_interesting$gene,sep="_")%in%paste(GRN_coef.dt_strong_repr$sender,GRN_coef.dt_strong_repr$target,sep="_"),]
dim(list_interesting_repr)
if(nrow(list_interesting_repr)>top) {
  list_interesting_repr <- head(list_interesting_repr,top)
}

list_interesting_act <- list_interesting[paste(list_interesting$tf,list_interesting$gene,sep="_")%in%paste(GRN_coef.dt_strong_act$sender,GRN_coef.dt_strong_act$target,sep="_"),]
dim(list_interesting_act)
if(nrow(list_interesting_act)>top) {
  list_interesting_act <- head(list_interesting_act,top)
}

#########################
## Get gene annotation ##
#########################

# GenomicRanges
genes.gr <- getGeneAnnotation(ArchRProject)[["genes"]]
genes.gr <- genes.gr[!is.na(genes.gr$symbol)]

####################
## Browser tracks ##
####################

sender_genes <- "EGR1"
target_genes <- "KLF5"

# opts$extend.upstream <- 1e4
# opts$extend.downstream <- 1e4
# opts$tileSize <- 50

opts$tileSize <- 30

interactions <- c("activating","repressing")
for (int in interactions){
  if(!dir.exists(sprintf("%s/%s",args$outdir,int))){dir.create(sprintf("%s/%s",args$outdir,int))}
  print(int)
  if (int=="activating"){
    sender_genes <- list_interesting_act$tf
    target_genes <- list_interesting_act$gene
  } else if (int=="repressing"){
    sender_genes <- list_interesting_repr$tf
    target_genes <- list_interesting_repr$gene
  } else {
    sender_genes <- "EGR1"
    target_genes <- "KLF5"
  }

  for (i in 1:length(sender_genes)) {

    sender_i <- sender_genes[i]
    target_i <- target_genes[i]

    print(sprintf("%s to %s",sender_i,target_i))

    clusters_for_plot <- as.numeric(gsub("^C","",ArchRProject$Clusters))
    clusters_for_plot[!(clusters_for_plot%in%c(1:5))] <- NA

    ArchRProject <- addCellColData(
      ArchRProject,
      data = clusters_for_plot,
      name = "clusters_for_plot",
      cells = rownames(getCellColData(ArchRProject)),
      force = TRUE
    )

    getCellColData(ArchRProject)

    gr_insilico_chip_seq_i <- gr_insilico_chip_seq[gr_insilico_chip_seq$tf==sender_i]

    gr_insilico_chip_seq_i <- gr_insilico_chip_seq_i[gr_insilico_chip_seq_i$gene==target_i]

    to.plot.gr <- gr_insilico_chip_seq_i
    gene.length <- abs(end(to.plot.gr) - start(to.plot.gr))
    start(to.plot.gr) <- start(to.plot.gr) - 600*10
    end(to.plot.gr) <- end(to.plot.gr) + 600*10

    # Plot
    p <- plotBrowserTrack(
      title=sprintf(": %s to %s interaction",sender_i,target_i),
      ArchRProj = ArchRProject, 
      region = to.plot.gr,
      geneSymbol = target_i,
      groupBy = "clusters_for_plot", 
      tileSize = args$tileSize,
      pal = opts$celltype.colors,
      features= gr_insilico_chip_seq_i,
      plotSummary = c("bulkTrack", "featureTrack", "geneTrack"),
      sizes = c(8, 1, 1),
      verbose=TRUE
    )

    r_expr_sender_i <- 0.01+c(0.018*minmax.10pct.normalisation(logcounts(sce)[sender_i,1:5]))
    r_expr_target_i <- 0.01+c(0.018*minmax.10pct.normalisation(logcounts(sce)[target_i,1:5]))

    # to be optimized for new dimensions of plots
    c_expr_sender_i <- circleGrob(name = "my_circle",
                          x = rep(0.94,5), y = 0.87-(c(c(1:5)-1)*0.1375), r = r_expr_sender_i,
                          gp = gpar(col = opts$celltype.colors[1:5], fill=opts$celltype.colors[1:5], lwd = 1))
    c_expr_target_i <- circleGrob(name = "my_circle",
                          x = rep(0.98,5), y = 0.87-(c(c(1:5)-1)*0.1375), r = r_expr_target_i,
                          gp = gpar(col = opts$celltype.colors[1:5], fill=opts$celltype.colors[1:5], lwd = 1))

    rect_y_axis_label <- rectGrob(x=0.03,y=0.4,width=0.05,height=1,gp = gpar(col = "white", fill="white", lwd = 1))
    rect_labels_boxes <- rectGrob(x=0.96,y=0.1,width=0.07,height=0.2,gp = gpar(col = "white", fill="white", lwd = 1))

    rect_grey_boxes_clust <- rectGrob(x=0.96,y=0.1,width=0.07,height=2,gp = gpar(col = "white", fill="white", lwd = 1))

    expr_label_sender_i <- textGrob(sprintf("%s",sender_i),x=0.94,y=0.955,gp=gpar(fontsize=7))
    expr_label_target_i <- textGrob(sprintf("%s",target_i),x=0.98,y=0.955,gp=gpar(fontsize=7))


    #pdf(sprintf("%s/%s/test_%s_to_%s_BrowserTrack.pdf",args$outdir,int,sender_i,target_i), width = 9, height = 5)
    pdf(sprintf("%s/test_%s_to_%s_BrowserTrack.pdf",args$outdir,sender_i,target_i), width = 9, height = 5)
    grid::grid.draw(p)
    grid::grid.draw(rect_grey_boxes_clust)
    grid::grid.draw(rect_labels_boxes)
    grid::grid.draw(rect_y_axis_label)

    grid::grid.draw(c_expr_sender_i)
    grid::grid.draw(c_expr_target_i)


    grid::grid.draw(expr_label_sender_i)
    grid::grid.draw(expr_label_target_i)
    dev.off()
  }
}
