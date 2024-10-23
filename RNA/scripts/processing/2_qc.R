############
##        ##
##  qc.R  ##
##        ##
############

source("/data/louisc/Project_Babraham/RNA/scripts/Settings.R")

#####################
## Define arguments ##
#####################

p <- ArgumentParser(description='')
p$add_argument('--metadata',       type="character",                    help='Metadata')
p$add_argument('--outputdir',       type="character",                    help='Output directory')
p$add_argument('--min_nFeature_RNA',       type="integer",                    help='Minimum number of expressed genes')
p$add_argument('--max_nFeature_RNA',       type="integer",                    help='Maximum number of expressed genes')
p$add_argument('--mitochondrial_percent_RNA',       type="integer",                    help='Maximum percentage of mitochondrial reads')
p$add_argument('--ribosomal_percent_RNA',       type="integer",                    help='Maximum percentage of ribosomal reads')
p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
p$add_argument('--cutoff_plot_nFeature_RNA',       type="integer",                    help='Cutoff maximal plotted values for expressed genes')
p$add_argument('--cutoff_plot_mito',       type="integer",                    help='Cutoff maximal plotted values for mitochondrial reads')
p$add_argument('--cutoff_plot_ribo',       type="integer",                    help='Cutoff maximal plotted values for ribosomal reads')
p$add_argument('--sort_samples',       type="logical",  default=TRUE,       help='Sort samples?')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# args$samples <- c("d0","d1","d3","d5","d7","d10","d14","d18","DE","NE","AME_E","AME_L")
# args$metadata <- "/data/louisc/Project_Babraham/RNA/metadata.txt.gz"
# args$min_nFeature_RNA <- 2000
# args$max_nFeature_RNA <- 12500
# args$mitochondrial_percent_RNA <- 35
# args$ribosomal_percent_RNA <- 30
# args$outputdir <- "/data/louisc/Project_Babraham/RNA/qc"
# args$sort_samples <- TRUE
# args$cutoff_plot_nFeature_RNA <- 15000
# args$cutoff_plot_mito <- 60
# args$cutoff_plot_ribo <- 50
## END TEST ##

#####################
## Define settings ##
#####################

# Sanity checks
# stopifnot(args$samples%in%opts$samples)

opts <- list()

###############
## Load data ##
###############

metadata <- fread(args$metadata) %>% 
  .[sample%in%args$samples] %>%
  # .[,pass_rnaQC:=nFeature_RNA>args$min_nFeature_RNA & nCount_RNA>2**args$log_nCount_RNA & mitochondrial_percent_RNA<args$mitochondrial_percent_RNA]
  .[,pass_rnaQC:=nFeature_RNA<=args$max_nFeature_RNA & nFeature_RNA>=args$min_nFeature_RNA & mitochondrial_percent_RNA<args$mitochondrial_percent_RNA & ribosomal_percent_RNA<args$ribosomal_percent_RNA]

table(metadata$pass_rnaQC)
table(metadata$pass_rnaQC,metadata$sample)

#####################
## Plot QC metrics ##
#####################

to.plot <- metadata %>% #.[pass_rnaQC==TRUE] %>%
  .[nFeature_RNA<=args$cutoff_plot_nFeature_RNA & mitochondrial_percent_RNA<=args$cutoff_plot_mito & ribosomal_percent_RNA<=args$cutoff_plot_ribo] %>% # remove massive outliers for plotting => does not make sense if you already put passQC at TRUE
  # .[,log_nFeature_RNA:=log10(nFeature_RNA)] %>%
  melt(id.vars=c("sample","cell"), measure.vars=c("nFeature_RNA","mitochondrial_percent_RNA","ribosomal_percent_RNA"))
# melt(id.vars=c("sample","cell"), measure.vars=c("nFeature_RNA"))

facet.labels <- c("nFeature_RNA" = "Num. of genes", "mitochondrial_percent_RNA" = "Mitochondrial %", "ribosomal_percent_RNA" = "Ribosomal %")

if (args$sort_samples){
  to.plot$sample <- as.factor(to.plot$sample)
  to.plot$sample <- factor(to.plot$sample,levels=args$samples)
}

## Box plot 

p <- ggplot(to.plot, aes_string(x="sample", y="value")) +
  geom_boxplot(outlier.shape=NA, coef=1) +
  facet_wrap(~variable, scales="free_y", nrow=1, labeller = as_labeller(facet.labels)) +
  scale_fill_manual(values=opts$stage.colors) +
  theme_classic() +
  theme(
    axis.text.y = element_text(colour="black",size=rel(1)),
    axis.text.x = element_text(colour="black",size=rel(0.65), angle=20, hjust=1, vjust=1),
    axis.title.x = element_blank()
  )

pdf(sprintf("%s/qc_metrics_boxplot.pdf",args$outputdir), width=9, height=5)
# pdf(sprintf("%s/qc_metrics_boxplot.pdf",args$outputdir))
print(p)
dev.off()

## histogram 

tmp <- data.table(
  variable = c("nFeature_RNA", "nFeature_RNA", "mitochondrial_percent_RNA", "ribosomal_percent_RNA"),
  value = c(args$min_nFeature_RNA, args$max_nFeature_RNA, args$mitochondrial_percent_RNA, args$ribosomal_percent_RNA)
)
# tmp <- data.table(
#     variable = c("nFeature_RNA"),
#     value = c(args$min_nFeature_RNA)
# )

p <- gghistogram(to.plot, x="value", fill="sample", bins=50) +
  # p <- ggdensity(to.plot, x="value", fill="sample") +
  geom_vline(aes(xintercept=value), linetype="dashed", data=tmp) +
  facet_wrap(~variable, scales="free", nrow=1) +
  theme(
    axis.text =  element_text(size=rel(0.8)),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size=rel(0.75))
  )

pdf(sprintf("%s/qc_metrics_histogram.pdf",args$outputdir), width=13, height=6)
print(p)
dev.off()

to.plot.ribo <- to.plot[to.plot$variable=="ribosomal_percent_RNA",]
tmp.ribo <- tmp[grepl("ribo",tmp$variable),]
p <- gghistogram(to.plot.ribo, x="value", fill="sample", bins=50) +
  # p <- ggdensity(to.plot, x="value", fill="sample") +
  geom_vline(aes(xintercept=value), linetype="dashed", data=tmp.ribo) +
  facet_wrap(~sample, scales="fixed") +
  theme(
    axis.text =  element_text(size=rel(0.8)),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size=rel(0.75))
  )

pdf(sprintf("%s/qc_metrics_histogram_ribo_per_sample.pdf",args$outputdir), width=13, height=6)
print(p)
dev.off()

to.plot.mito <- to.plot[to.plot$variable=="mitochondrial_percent_RNA",]
tmp.mito <- tmp[grepl("mito",tmp$variable),]
p <- gghistogram(to.plot.mito, x="value", fill="sample", bins=50) +
  # p <- ggdensity(to.plot, x="value", fill="sample") +
  geom_vline(aes(xintercept=value), linetype="dashed", data=tmp.mito) +
  facet_wrap(~sample, scales="fixed") +
  theme(
    axis.text =  element_text(size=rel(0.8)),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size=rel(0.75))
  )

pdf(sprintf("%s/qc_metrics_histogram_mito_per_sample.pdf",args$outputdir), width=13, height=6)
print(p)
dev.off()

to.plot.RNA <- to.plot[to.plot$variable=="nFeature_RNA",]
tmp.RNA <- tmp[grepl("nFeature_RNA",tmp$variable),]
p <- gghistogram(to.plot.RNA, x="value", fill="sample", bins=50) +
  # p <- ggdensity(to.plot, x="value", fill="sample") +
  geom_vline(aes(xintercept=value), linetype="dashed", data=tmp.RNA) +
  facet_wrap(~sample, scales="fixed") +
  theme(
    axis.text =  element_text(size=rel(0.8)),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size=rel(0.75))
  )

pdf(sprintf("%s/qc_metrics_histogram_nFeatureRNA_per_sample.pdf",args$outputdir), width=13, height=6)
print(p)
dev.off()


#########################################################
## Plot fraction of cells that pass QC for each sample ##
#########################################################

to.plot <- metadata %>%
  .[,mean(pass_rnaQC,na.rm=T),by=c("sample")]

p <- ggbarplot(to.plot, x="sample", y="V1") +
  scale_fill_manual(values=opts$stage.colors) +
  labs(x="", y="Fraction of cells that pass QC (RNA)") +
  # facet_wrap(~stage)
  theme(
    legend.position = "none",
    axis.text.y = element_text(colour="black",size=rel(0.8)),
    axis.text.x = element_text(colour="black",size=rel(0.65), angle=20, hjust=1, vjust=1),
  )

pdf(sprintf("%s/qc_metrics_barplot.pdf",args$outputdir), width=6, height=5)
print(p)
dev.off()

##########
## Save ##
##########

fwrite(metadata, paste0(args$outputdir,"/sample_metadata_after_qc.txt.gz"), quote=F, na="NA", sep="\t")
