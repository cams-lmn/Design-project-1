##################
##              ##
##  qc_archr.R  ##
##              ##
##################

source("/data/homes/louisc/Project_Babraham/ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',    type="character",    help='ArchR directory')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--metadata_rna',    type="character",    help='metadata file')
p$add_argument('--outdir',     type="character",    help='Output directory')
p$add_argument('--min_tss_enrichment',     type="integer",    default=2,   help='Minimum TSS enrichment')
p$add_argument('--max_tss_enrichment',     type="integer",    default=50,   help='Minimum TSS enrichment')
p$add_argument('--min_number_fragments',     type="integer",    default=100,    help='Maximum number of ATAC fragments')
p$add_argument('--max_number_fragments',     type="integer",    default=1e6,    help='Maximum number of ATAC fragments')
p$add_argument('--max_blacklist_ratio',     type="double",    default=0.05,    help='Maximum Blacklist Ratio')
p$add_argument('--cutoff_nFrag',     type="integer",    default=1e6,    help='Maximum number of ATAC fragments')
p$add_argument('--cutoff_TSS_enrichment',     type="integer",    default=50,    help='Maximum number of ATAC fragments')
p$add_argument('--genome',          type="character", default="hg38",      help='Genome')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
p$add_argument('--samples', type="character",  nargs='+',      help='Samples')
p$add_argument("--sort_samples",  type="logical",default=T, help="Should samples be sorted?")
args <- p$parse_args(commandArgs(TRUE))

# START TEST ##
# args <- list()
# args$samples <- c("d0","d1","d3","d5","d7","d10","d14","d18","DE","NE","AME_E","AME_L")
# args$archr_directory <- "/data/homes/louisc/Project_Babraham/ATAC/archR"
# args$metadata <- "/data/homes/louisc/Project_Babraham/ATAC/archR/sample_metadata_after_archR.txt.gz"
# args$metadata_rna <- "/data/louisc/homes/Project_Babraham/RNA/mapping/sample_metadata_after_doublets.txt.gz"
# args$min_tss_enrichment <- 3
# args$max_tss_enrichment <- 25
# args$min_number_fragments <- 6000
# args$max_number_fragments <- 1e6
# args$max_blacklist_ratio <- 0.03
# args$max_fragment_size <- 400
# args$cutoff_nFrag <- 10^(6.2)
# args$cutoff_TSS_enrichment <- 35
# args$sort_samples <- TRUE
# args$threads <- 16
# args$genome <- "hg38"
# args$outdir <- "/data/homes/louisc/Project_Babraham/ATAC/archR/qc"
## END TEST ##

#####################
## Define settings ##
#####################

# Options
opts$chr <- paste0("chr",1:22)

print("To really recreate qc plots run snakemake for create_arrow_files rule!")
warning("To really recreate qc plots run snakemake for create_arrow_files rule!")

########################
## Load cell metadata ##
########################

sample_metadata <- fread(args$metadata)

# temporary
# sample_metadata[is.na(stage),stage:=strsplit(sample,"_") %>% map_chr(1)]

########################
## Load ArchR project ##
########################

# source(here::here("atac/archR/load_archR_project.R"))

setwd(args$archr_directory)

addArchRGenome(args$genome)
addArchRThreads(threads = args$threads)

ArchRProject <- loadArchRProject(args$archr_directory)[sample_metadata$cell]
# ArchRProject <- loadArchRProject(args$archr_directory)
# ArchRProject <- ArchRProject[sample_metadata$cell[sample_metadata$cell%in%ArchRProject$cellNames]]

print(ArchRProject)

##################
## Subset ArchR ##
##################

# Subset chr for faster computations
tss.granges <- getTSS(ArchRProject)
tss.granges <- tss.granges[seqnames(tss.granges)%in%opts$chr]

print(tss.granges)

####################
## Cutoffs object ##
####################

print("Cutoff object")

features <- c("TSSEnrichment_atac", "log_nFrags_atac", "BlacklistRatio_atac")

cutoffs <- data.table(
  variable = c("TSSEnrichment_atac", "TSSEnrichment_atac", "log_nFrags_atac", "log_nFrags_atac","BlacklistRatio_atac"),
  value = c(args$min_tss_enrichment, args$max_tss_enrichment, log10(args$min_number_fragments), 
            log10(args$max_number_fragments),args$max_blacklist_ratio)
)

print(cutoffs)

#########################
## Plot TSS Enrichment ##
#########################

print("TSS enrichment")

# Calculate TSS Enrichment
data_tss.dt <- args$samples %>% map(function(i) {
  plotTSSEnrichment(
    ArchRProj = ArchRProject[ArchRProject$Sample==i,], 
    groupBy = "Sample", 
    returnDF = TRUE,
    TSS = tss.granges
  ) %>% as.data.table %>% return
}) %>% rbindlist %>% 
  setnames("group","sample") %>% 
  melt(id.vars=c("sample","x"))
fwrite(data_tss.dt, sprintf("%s/qc_TSSenrichment.txt.gz",args$outdir))

# Plot TSS Enrichment
to.plot_tss.dt <- data_tss.dt %>% 
  merge(unique(sample_metadata[,c("sample")])) %>%
  .[,.(value=mean(value)), by = c("sample","x","variable")] %>%
  .[variable=="normValue"]

print(dim(to.plot_tss.dt))
print(head(to.plot_tss.dt))

p <- ggline(to.plot_tss.dt, x="x", y="value", numeric.x.axis=T, color = "sample", palette = opts$color_scheme, plot_type="l") + 
  #scale_x_continuous(breaks=c(-2000,-1000,0,1000,2000))+
  facet_wrap(~sample, scales="fixed") +
  labs(x="Distance from TSS (bp)", y="TSS enrichment (normalised)") +
  theme(
    axis.text.y = element_text(size=rel(0.65), color="black"),
    axis.text.x = element_text(size=rel(0.6), color="black"),
    axis.title = element_text(size=rel(0.75), color="black"),
    legend.position = "none",
    legend.title = element_blank()
  ) 

output_plot(p, file.path(args$outdir,"qc_TSSenrichment"), width=6, height=5)

# Plot fragment sizes - no diff
print("Subset non differentiated cells")
print(dim(to.plot_tss.dt))
to.plot_tss.dt_nodiff <- to.plot_tss.dt[grepl("d[0-9]+",to.plot_tss.dt$sample),]
print(dim(to.plot_tss.dt_nodiff))

p <- ggline(to.plot_tss.dt_nodiff, x="x", y="value", numeric.x.axis=T, color = "sample", palette = opts$color_scheme, plot_type="l") +
  labs(x="Distance from TSS (bp)", y="TSS enrichment (normalised)") +
  theme(
    axis.text.y = element_text(size=rel(0.65), color="black"),
    axis.text.x = element_text(size=rel(0.6), color="black"),
    axis.title = element_text(size=rel(0.75), color="black"),
    legend.position = "none",
    legend.title = element_blank()
  ) + scale_x_continuous(breaks=c(-2000,-1000,0,1000,2000))

output_plot(p, file.path(args$outdir,"qc_TSSenrichment_nodiff_singleplot"), width=6, height=5)

#####################################
## Plot Fragment size distribution ##
#####################################

print("Fragment size distribution")

# Calculate fragment sizes
data_fragmentsize.dt <- args$samples %>% map(function(i) {
  plotFragmentSizes(
    ArchRProj = ArchRProject[ArchRProject$Sample==i,], 
    groupBy = "Sample", 
    returnDF = TRUE
  ) %>% as.data.table %>% return
}) %>% rbindlist %>% 
  setnames("group","sample")

fwrite(data_fragmentsize.dt, sprintf("%s/qc_FragmentSizeDistribution.txt.gz",args$outdir))

# Plot fragment sizes
to.plot_fragmentsize.dt <- data_fragmentsize.dt %>% 
  merge(unique(sample_metadata[,c("sample")])) %>%
  .[,.(fragmentPercent=mean(fragmentPercent)), by = c("sample","fragmentSize")] #%>%

print(dim(to.plot_fragmentsize.dt))
print(head(to.plot_fragmentsize.dt))

p <- ggline(to.plot_fragmentsize.dt, x="fragmentSize", y="fragmentPercent", color = "sample", palette = opts$color_scheme, plot_type="l") +
  facet_wrap(~sample, scales="fixed") +
  labs(x="Fragment Size (bp)", y="Percentage of fragments (%)") +
  theme(
    axis.text = element_text(size=rel(0.55)),
    axis.title = element_text(size=rel(0.75)),
    legend.position = "none",
    legend.title = element_blank()
  )

output_plot(p, file.path(args$outdir,"qc_FragmentSizeDistribution"), width=8, height=4)

# Plot fragment sizes - no diff
print("Subset non differentiated cells")
print(dim(to.plot_fragmentsize.dt))
to.plot_fragmentsize.dt_nodiff <- to.plot_fragmentsize.dt[grepl("d[0-9]+",to.plot_fragmentsize.dt$sample),]
print(dim(to.plot_fragmentsize.dt_nodiff))

p <- ggline(to.plot_fragmentsize.dt_nodiff, x="fragmentSize", y="fragmentPercent", color = "sample", palette = opts$color_scheme, plot_type="l") +
  labs(x="Fragment Size (bp)", y="Percentage of fragments (%)") +
  theme(
    axis.text = element_text(size=rel(0.55)),
    axis.title = element_text(size=rel(0.75)),
    legend.position = "none",
    legend.title = element_blank()
  )

output_plot(p, file.path(args$outdir,"qc_FragmentSizeDistribution_nodiff_singleplot"), width=8, height=4)

##################################
## Plot histogram of QC metrics ##
##################################

print("Histograms")

to.plot <- sample_metadata %>%
  .[!is.na(nFrags_atac)] %>%
  .[,log_nFrags_atac:=log10(nFrags_atac)] %>%
  melt(id.vars=c("sample","cell"), measure.vars=c("TSSEnrichment_atac","log_nFrags_atac","BlacklistRatio_atac"))

print(head(to.plot))

if (args$sort_samples){
  to.plot$sample <- as.factor(to.plot$sample)
  to.plot$sample <- factor(to.plot$sample,levels=args$samples)
}

p <- gghistogram(to.plot, x="value", y="..density..", fill="sample", bins=70) +
  geom_vline(aes(xintercept=value), linetype="dashed", data=cutoffs, linewidth=1.5) +
  scale_fill_manual(values=opts$color_scheme) +
  facet_wrap(~variable, scales="free") +
  theme(
    axis.text =  element_text(size=rel(0.8)),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size=rel(0.75))
)

output_plot(p,sprintf("%s/qc_metrics_histogram",args$outdir), width=8, height=5)

# per feature

for (feature_i in features){
  print(feature_i)

  to.plot_i <- to.plot[grepl(feature_i,to.plot$variable),]

  p <- gghistogram(to.plot_i, x="value", y="..density..", fill="sample", bins=40) +
    scale_fill_manual(values=opts$color_scheme) +
    facet_wrap(~sample, scales="fixed") +
    theme(
      axis.text =  element_text(size=rel(0.8)),
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      legend.position = "top",
      legend.text = element_text(size=rel(0.75))
    )

  cutoffs_i <- cutoffs$value[grepl(feature_i,cutoffs$variable)]
  print(cutoffs_i)
  p <- p  + geom_vline(xintercept=cutoffs_i, linetype="dashed", col="red", linewidth=1.5)

  output_plot(p,sprintf("%s/qc_metrics_histogram_%s_per_sample",args$outdir,feature_i), width=8, height=5)
}

# no diff - per feature

print("Subset non differentiated cells")
print(dim(to.plot))
to.plot_nodiff <- to.plot[grepl("d[0-9]+",to.plot$sample),]
print(dim(to.plot_nodiff))

for (feature_i in features){
  print(feature_i)

  to.plot_nodiff_i <- to.plot_nodiff[grepl(feature_i,to.plot_nodiff$variable),]

  p <- gghistogram(to.plot_nodiff_i, x="value", y="..density..", fill="sample", bins=40) +
    scale_fill_manual(values=opts$color_scheme) +
    facet_wrap(~sample, scales="fixed") +
    theme(
      axis.text =  element_text(size=rel(0.8)),
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      legend.position = "top",
      legend.text = element_text(size=rel(0.75))
    )

  cutoffs_i <- cutoffs$value[grepl(feature_i,cutoffs$variable)]
  print(cutoffs_i)
  p <- p  + geom_vline(xintercept=cutoffs_i, linetype="dashed", linewidth=1.5)

  output_plot(p, sprintf("%s/qc_metrics_histogram_%s_no_diff_per_sample",args$outdir,feature_i), width=8, height=5)
}

################################
## Plot boxplot of QC metrics ##
################################

print("Boxplots")

print(head(to.plot))

# per feature

for (feature_i in features){
  print(feature_i)

  to.plot_i <- to.plot[grepl(feature_i,to.plot$variable),]

  p <- ggplot(to.plot_nodiff_i, aes_string(x="sample", y="value", fill="sample")) +
    geom_boxplot(outlier.shape=NA, coef=1) +
    scale_fill_manual(values=opts$color_scheme) +
    geom_jitter(color="black", size=0.05, alpha=0.1) +
    theme_classic() +
    theme(
      axis.text =  element_text(size=rel(0.8)),
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      legend.position = "top",
      legend.text = element_text(size=rel(0.75))
    )

  cutoffs_i <- cutoffs$value[grepl(feature_i,cutoffs$variable)]
  print(cutoffs_i)
  p <- p  + geom_hline(yintercept=cutoffs_i, linetype="dashed", col="red", linewidth=1.5)

  output_plot(p, sprintf("%s/qc_metrics_boxplot_%s_per_sample",args$outdir,feature_i), width=8, height=5)
}

# no diff - per feature

print("Subset non differentiated cells")
print(dim(to.plot))
to.plot_nodiff <- to.plot[grepl("d[0-9]+",to.plot$sample),]
print(dim(to.plot_nodiff))

for (feature_i in features){
  print(feature_i)

  to.plot_nodiff_i <- to.plot_nodiff[grepl(feature_i,to.plot_nodiff$variable),]

p <- ggplot(to.plot_nodiff_i, aes_string(x="sample", y="value", fill="sample")) +
    geom_boxplot(outlier.shape=NA, coef=1) +
    scale_fill_manual(values=opts$color_scheme) +
    geom_jitter(color="black", size=0.05, alpha=0.1) +
    theme_classic() +
    theme(
      axis.text =  element_text(size=rel(0.8)),
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      legend.position = "top",
      legend.text = element_text(size=rel(0.75))
    )

  cutoffs_i <- cutoffs$value[grepl(feature_i,cutoffs$variable)]
  print(cutoffs_i)
  p <- p  + geom_hline(yintercept=cutoffs_i, linetype="dashed", col="red", linewidth=1.5)

  output_plot(p, sprintf("%s/qc_metrics_boxplot_%s_per_sample",args$outdir,feature_i), width=8, height=5)
}

#############
## Call QC ##
#############

sample_metadata %>%
  .[,pass_atacQC:=TSSEnrichment_atac>=args$min_tss_enrichment & TSSEnrichment_atac<=args$max_tss_enrichment & 
      nFrags_atac>=args$min_number_fragments & nFrags_atac<=args$max_number_fragments & 
      BlacklistRatio_atac<=args$max_blacklist_ratio] %>%
  .[is.na(pass_atacQC),pass_atacQC:=FALSE]

print(head(sample_metadata))
print(summary(sample_metadata))

print(sample_metadata[,mean(pass_atacQC,na.rm=T),by="sample"])
# print(sample_metadata[,mean(is.na(nFrags_atac)),by="sample"])

# Filter low quality cells that did not pass RNA QC
sample_metadata[pass_rnaQC==FALSE, pass_atacQC:=FALSE]

# Save
fwrite(sample_metadata, file.path(args$outdir,"sample_metadata_after_qc.txt.gz"), quote=F, na="NA", sep="\t")

###########################################
## Plot QC statistics after QC filtering ##
###########################################

print("Barplots")

# Barplot of the fraction of cells that pass QC for each sample
to.plot <- sample_metadata %>%
  .[,mean(pass_atacQC,na.rm=T),by=c("sample")]

if (args$sort_samples){
  to.plot$sample <- as.factor(to.plot$sample)
  to.plot$sample <- factor(to.plot$sample,levels=args$samples)
}

print(to.plot)

p <- ggbarplot(to.plot, x="sample", y="V1", fill="sample") +
  scale_fill_manual(values=opts$color_scheme) +
  labs(x="", y="Fraction of cells that pass ATAC QC") +
  coord_cartesian(ylim=c(0,1)) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(colour="black",size=rel(0.8)),
    axis.text.x = element_text(colour="black",size=rel(0.65), angle=20, hjust=1, vjust=1),
  )

output_plot(p, sprintf("%s/qc_metrics_barplot",args$outdir), width=6, height=5)

# Barplot of the fraction of cells that pass QC for each sample - no diff
print("Subset non differentiated cells")
print(dim(to.plot))
to.plot_nodiff <- to.plot[grepl("d[0-9]+",to.plot$sample),]
print(dim(to.plot_nodiff))

p <- ggbarplot(to.plot_nodiff, x="sample", y="V1", fill="sample") +
  scale_fill_manual(values=opts$color_scheme) +
  labs(x="", y="Fraction of cells that pass ATAC QC") +
  coord_cartesian(ylim=c(0,1)) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(colour="black",size=rel(0.8)),
    axis.text.x = element_text(colour="black",size=rel(0.65), angle=20, hjust=1, vjust=1),
  )

output_plot(p, sprintf("%s/qc_metrics_barplot_nodiff",args$outdir), width=6, height=5)

# Boxplots of QC metrics
to.plot <- sample_metadata %>%
  .[pass_atacQC==TRUE] %>%
  .[nFrags_atac<=args$cutoff_nFrag & TSSEnrichment_atac<=args$cutoff_TSS_enrichment] %>% # remove massive outliers for plotting
  .[,log_nFrags_atac:=log10(nFrags_atac)] %>%
  # melt(id.vars=c("sample","cell"), measure.vars=c("TSSEnrichment_atac","log_nFrags","BlacklistRatio_atac"))
  melt(id.vars=c("sample","cell","sample"), measure.vars=c("TSSEnrichment_atac","log_nFrags_atac"))

facet.labels <- c("log_nFrags_atac" = "Num. of fragments (log10)", "TSSEnrichment_atac" = "TSS enrichment")

if (args$sort_samples){
  to.plot$sample <- as.factor(to.plot$sample)
  to.plot$sample <- factor(to.plot$sample,levels=args$samples)
}

print(to.plot)

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

output_plot(p, sprintf("%s/qc_metrics_boxplot",args$outdir), width=9, height=5)

###############################################################
## QC statistics after QC filtering & comparison with RNAseq ##
###############################################################

print("No filtering")

# data after QC
if (args$sort_samples){
  table(sample_metadata$sample)[args$samples]
} else {
  table(sample_metadata$sample)
}

# Comparison RNA after qc 
sample_metadata_RNA <- fread(args$metadata_rna)
sample_metadata_RNA <- sample_metadata_RNA[sample_metadata_RNA$pass_rnaQC,]

print(sprintf("%s cells have both RNA expression and chromatin accessibility measurements",length(intersect(sample_metadata_RNA$cell,sample_metadata$cell))))
print(sprintf("%s cells have RNA expression, but do not have chromatin accessibility measurements",sum(!sample_metadata_RNA$cell%in%sample_metadata$cell)))
print(sprintf("%s cells have chromatin accessibility, but do not have RNA expression measurements",sum(!sample_metadata$cell%in%sample_metadata_RNA$cell)))

print("Pass ATAC QC")
sample_metadata_passATACQC <- sample_metadata[sample_metadata$pass_atacQC,]
# data after QC
if (args$sort_samples){
  table(sample_metadata_passATACQC$sample)[args$samples]
} else {
  table(sample_metadata_passATACQC$sample)
}

print(sprintf("%s cells have both RNA expression and chromatin accessibility measurements",length(intersect(sample_metadata_RNA$cell,sample_metadata_passATACQC$cell))))
print(sprintf("%s cells have RNA expression, but do not have chromatin accessibility measurements",sum(!sample_metadata_RNA$cell%in%sample_metadata_passATACQC$cell)))
print(sprintf("%s cells have chromatin accessibility, but do not have RNA expression measurements",sum(!sample_metadata_passATACQC$cell%in%sample_metadata_RNA$cell)))

# Including doublet detection
print("Pass doublet detection")
sample_metadata_doublet <- sample_metadata[sample_metadata$doublet_call==F,]
if (args$sort_samples){
  table(sample_metadata_doublet$sample)[args$samples]
} else {
  table(sample_metadata_doublet$sample)
}
sample_metadata_RNA <- sample_metadata_RNA[sample_metadata_RNA$doublet_call==F,]

print(sprintf("%s cells have both RNA expression and chromatin accessibility measurements",length(intersect(sample_metadata_RNA$cell,sample_metadata_doublet$cell))))
print(sprintf("%s cells have RNA expression, but do not have chromatin accessibility measurements",sum(!sample_metadata_RNA$cell%in%sample_metadata_doublet$cell)))
print(sprintf("%s cells have chromatin accessibility, but do not have RNA expression measurements",sum(!sample_metadata_doublet$cell%in%sample_metadata_RNA$cell)))

# Including doublet detection
print("Pass doublet detection & ATAC QC")
sample_metadata_both <- sample_metadata[sample_metadata$doublet_call==F & sample_metadata$pass_atacQC,]
if (args$sort_samples){
  table(sample_metadata_both$sample)[args$samples]
} else {
  table(sample_metadata_both$sample)
}
sample_metadata_RNA <- sample_metadata_RNA[sample_metadata_RNA$doublet_call==F,]

print(sprintf("%s cells have both RNA expression and chromatin accessibility measurements",length(intersect(sample_metadata_RNA$cell,sample_metadata_both$cell))))
print(sprintf("%s cells have RNA expression, but do not have chromatin accessibility measurements",sum(!sample_metadata_RNA$cell%in%sample_metadata_both$cell)))
print(sprintf("%s cells have chromatin accessibility, but do not have RNA expression measurements",sum(!sample_metadata_both$cell%in%sample_metadata_RNA$cell)))

#########################################################################
## QC statistics after QC filtering & comparison with RNAseq - no diff ##
#########################################################################

print("No filtering")

sample_metadata_nodiff <- sample_metadata[grepl("d[0-9]+",sample_metadata$sample),]

if (args$sort_samples){
  table(sample_metadata_nodiff$sample)[args$samples[grepl("d[0-9]+",args$samples)]]
} else {
  table(sample_metadata_nodiff$sample)
}

# Comparison RNA after qc 
sample_metadata_RNA <- fread(args$metadata_rna)
sample_metadata_nodiff_RNA <- sample_metadata_RNA[grepl("d[0-9]+",sample_metadata_RNA$sample),]
sample_metadata_nodiff_RNA <- sample_metadata_nodiff_RNA[sample_metadata_nodiff_RNA$pass_rnaQC,]

print(sprintf("%s cells have both RNA expression and chromatin accessibility measurements",length(intersect(sample_metadata_nodiff_RNA$cell,sample_metadata_nodiff$cell))))
print(sprintf("%s cells have RNA expression, but do not have chromatin accessibility measurements",sum(!sample_metadata_nodiff_RNA$cell%in%sample_metadata_nodiff$cell)))
print(sprintf("%s cells have chromatin accessibility, but do not have RNA expression measurements",sum(!sample_metadata_nodiff$cell%in%sample_metadata_nodiff_RNA$cell)))

print("Pass ATAC QC")
sample_metadata_nodiff_passATACQC <- sample_metadata_nodiff[sample_metadata_nodiff$pass_atacQC,]
# data after QC
if (args$sort_samples){
  table(sample_metadata_nodiff_passATACQC$sample)[args$samples[grepl("d[0-9]+",args$samples)]]
} else {
  table(sample_metadata_nodiff_passATACQC$sample)
}

print(sprintf("%s cells have both RNA expression and chromatin accessibility measurements",length(intersect(sample_metadata_nodiff_RNA$cell,sample_metadata_nodiff_passATACQC$cell))))
print(sprintf("%s cells have RNA expression, but do not have chromatin accessibility measurements",sum(!sample_metadata_nodiff_RNA$cell%in%sample_metadata_nodiff_passATACQC$cell)))
print(sprintf("%s cells have chromatin accessibility, but do not have RNA expression measurements",sum(!sample_metadata_nodiff_passATACQC$cell%in%sample_metadata_nodiff_RNA$cell)))

# Including doublet detection
print("Pass doublet detection")
sample_metadata_nodiff_doublet <- sample_metadata_nodiff[sample_metadata_nodiff$doublet_call==F,]
if (args$sort_samples){
  table(sample_metadata_nodiff_doublet$sample)[args$samples[grepl("d[0-9]+",args$samples)]]
} else {
  table(sample_metadata_nodiff_doublet$sample)
}
sample_metadata_nodiff_RNA <- sample_metadata_nodiff_RNA[sample_metadata_nodiff_RNA$doublet_call==F,]

print(sprintf("%s cells have both RNA expression and chromatin accessibility measurements",length(intersect(sample_metadata_nodiff_RNA$cell,sample_metadata_nodiff_doublet$cell))))
print(sprintf("%s cells have RNA expression, but do not have chromatin accessibility measurements",sum(!sample_metadata_nodiff_RNA$cell%in%sample_metadata_nodiff_doublet$cell)))
print(sprintf("%s cells have chromatin accessibility, but do not have RNA expression measurements",sum(!sample_metadata_nodiff_doublet$cell%in%sample_metadata_nodiff_RNA$cell)))

# Including doublet detection
print("Pass doublet detection & ATAC QC")
sample_metadata_nodiff_both <- sample_metadata_nodiff[sample_metadata_nodiff$doublet_call==F & sample_metadata_nodiff$pass_atacQC,]
if (args$sort_samples){
  table(sample_metadata_nodiff_both$sample)[args$samples[grepl("d[0-9]+",args$samples)]]
} else {
  table(sample_metadata_nodiff_both$sample)
}
sample_metadata_nodiff_RNA <- sample_metadata_nodiff_RNA[sample_metadata_nodiff_RNA$doublet_call==F,]

print(sprintf("%s cells have both RNA expression and chromatin accessibility measurements",length(intersect(sample_metadata_nodiff_RNA$cell,sample_metadata_nodiff_both$cell))))
print(sprintf("%s cells have RNA expression, but do not have chromatin accessibility measurements",sum(!sample_metadata_nodiff_RNA$cell%in%sample_metadata_nodiff_both$cell)))
print(sprintf("%s cells have chromatin accessibility, but do not have RNA expression measurements",sum(!sample_metadata_nodiff_both$cell%in%sample_metadata_nodiff_RNA$cell)))
