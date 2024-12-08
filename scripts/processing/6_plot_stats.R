####################
##                ##
##  plot_stats.R  ##
##                ##
####################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA/scripts/Settings.R")

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--sce',             type="character",                               help='SingleCellExperiment file')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--samples',       type="character",  default="all",  nargs='+',  help='Samples to plot')
p$add_argument('--outdir',          type="character",                help='Output directory')
p$add_argument("--sort_samples",  type="logical",default=TRUE, help="Should samples be sorted?")
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

# These settings are important for consistency with ArchR, which provides little flexibility to edit cell names
opts <- list()

#####################
## Parse arguments ##
#####################

# I/O
dir.create(args$outdir, showWarnings = F)

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata) %>% 
  .[pass_rnaQC==TRUE & doublet_call==FALSE & sample%in%args$samples]

# sample_metadata <- fread(io$metadata) %>%  .[(pass_rnaQC==TRUE | pass_atacQC==TRUE)]
if (args$sort_samples){
  table(sample_metadata$sample)[args$samples]
} else {
  table(sample_metadata$sample)
}


###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(args$sce, cells=sample_metadata$cell, normalise = FALSE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

################################################
## Plot distribution of UMI counts per sample ##
################################################

to.plot <- args$samples %>% map(function(i) {
  mtx <- as.matrix(counts(sce[,sce$sample==i])[,1:min(ncol(sce[,sce$sample==i]),100)])
  data.table(table(mtx[mtx<=10])) %>% setnames(c("value","N")) %>% .[,sample:=i]
}) %>% rbindlist

to.plot[N>=1e6,N:=1e6]

p <- ggbarplot(to.plot, x="value", y="N", fill="gray70") +
  facet_wrap(~sample) +
  theme(
    axis.text =  element_text(size=rel(0.65)),
    axis.title.x = element_blank()
  )

pdf(paste0(args$outdir,"/umi_distribution_per_sample.pdf"), width=8, height=5)
print(p)
dev.off()

#########################################
## Barplots number of cells per sample ##
#########################################

to.plot <- sample_metadata[,.N,by=c("sample")]
fwrite(to.plot[,c("sample","N")], file.path(args$outdir,"ncells_per_sample.txt.gz"), sep="\t", quote=F)

if (args$sort_samples){
  to.plot$sample <- as.factor(to.plot$sample)
  to.plot$sample <- factor(to.plot$sample,levels=args$samples)
}

p <- ggbarplot(to.plot, x = "sample", y = "N") +
  # scale_fill_manual(values=opts$stage.colors) +
  labs(x="", y="Number of cells (after QC)") +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.y = element_text(colour="black",size=rel(0.7)),
    axis.text.x = element_text(colour="black",size=rel(0.8), angle=40, hjust=1)
  )

pdf(paste0(args$outdir,"/ncells_per_sample.pdf"), width=8, height=5)
print(p)
dev.off()

#########################################
## Barplots number of reads per sample ##
#########################################

to.plot <- data.table(cell=colnames(sce), sample=sce$sample, stage=sce$stage, N=colSums(counts(sce))) %>%
  .[,.(N=sum(N)),by=c("sample")]
fwrite(to.plot[,c("sample","N")], file.path(args$outdir,"nreads_per_sample.txt.gz"), sep="\t", quote=F)

if (args$sort_samples){
  to.plot$sample <- as.factor(to.plot$sample)
  to.plot$sample <- factor(to.plot$sample,levels=args$samples)
}

p <- ggbarplot(to.plot, x = "sample", y = "N") +
  # scale_fill_manual(values=opts$stage.colors) +
  labs(x="", y="Number of reads (after QC)") +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.y = element_text(colour="black",size=rel(0.7)),
    axis.text.x = element_text(colour="black",size=rel(0.8), angle=40, hjust=1)
  )

pdf(paste0(args$outdir,"/nreads_per_sample.pdf"), width=8, height=5)
print(p)
dev.off()