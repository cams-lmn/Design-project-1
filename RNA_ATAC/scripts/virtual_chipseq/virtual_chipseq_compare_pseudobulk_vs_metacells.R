#########################################################
##                                                     ##
##  virtual_chipseq_compare_pseudobulk_vs_metacells.R  ##
##                                                     ##
#########################################################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA_ATAC/scripts/Settings.R")

#####################
## Define settings ##
#####################

p <- ArgumentParser(description='')
p$add_argument('--motif_annotation',  type="character",              help='Motif annotation')
p$add_argument('--virtual_chip_pseudobulk_dir',  type="character", help='Virtual chipSeq directory (pseudobulk)')
p$add_argument('--virtual_chip_metacells_dir',  type="character",  help='Virtual chipSeq directory (metacells)')
p$add_argument('--min_chip_score',  type="double",  help='Minimal chip score')
p$add_argument('--outdir',          type="character",                help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

# opts$TFs <- c("Esr1", "Gata3", "Foxa1", "Stat3", "E2f1", "Myc")
opts$TFs <- ""
dir.create(args$outdir, showWarnings = F, recursive = T)

###################################
## Load virtual ChIP-seq library ##
###################################

# pseudobulk
virtual_chip_pseudobulk.mtx <- readRDS(file.path(args$virtual_chip_pseudobulk_dir,"virtual_chip_matrix.rds"))

# metacells
virtual_chip_metacells.mtx <- readRDS(file.path(args$virtual_chip_metacells_dir,"virtual_chip_matrix.rds"))

################
## Parse data ##
################

opts$TFs <- opts$TFs[opts$TFs%in%colnames(virtual_chip_pseudobulk.mtx)]
opts$TFs <- opts$TFs[opts$TFs%in%colnames(virtual_chip_metacells.mtx)]

stopifnot(rownames(virtual_chip_metacells.mtx)==rownames(virtual_chip_pseudobulk.mtx))

####################################
## Compare distribution of scores ##
####################################

to.plot <- rbind(
  data.table(value=as.numeric(virtual_chip_metacells.mtx), class="metacells"),
  data.table(value=as.numeric(virtual_chip_pseudobulk.mtx), class="pseudobulk")
) %>% .[value!=0]

p <- gghistogram(to.plot, x="value", fill="class", bins=100) +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size=rel(0.75), color="black")
  )

pdf(sprintf("%s/histogram_scores.pdf",args$outdir), width=7, height=5)
print(p)
dev.off()

# Set x-axis breaks
x_breaks <- pretty(to.plot$value, n = 10)

# Histogram plot (2 subplots via facet)
hist_plot <- ggplot(to.plot, aes(x = value)) +
  geom_histogram(bins = 150, fill = "gray60") +
  facet_wrap(~class, ncol = 2) +
  scale_x_continuous(breaks = x_breaks) +
  labs(title = "Histogram of Virtual ChIP-seq Scores",
       x = "Virtual ChIP-seq score", y = "Count") +
  theme_minimal()

# Density plot (2 subplots via facet)
density_plot <- ggplot(to.plot, aes(x = value)) +
  geom_density(fill = "gray60", alpha = 0.5) +
  facet_wrap(~class, ncol = 2) +
  scale_x_continuous(breaks = x_breaks) +
  labs(title = "Density Plot of Virtual ChIP-seq Scores",
       x = "Virtual ChIP-seq score", y = "Density") +
  theme_minimal()

pdf(sprintf("%s/histogram_density_plots.pdf", args$outdir), width = 10, height = 8)
print(hist_plot)
print(density_plot)
dev.off()

for (i in opts$TF){

  ########################################
  ## Scatterplot of distribution scores ##
  ########################################
  
  to.plot <- data.table(metacells = virtual_chip_metacells.mtx[,i], 
                        pseudobulk = virtual_chip_pseudobulk.mtx[,i]) %>%
    .[abs(pseudobulk)>0.01]

  p <- ggscatter(to.plot, x="pseudobulk", y="metacells", size=0.5) +
    geom_abline(slope=1, intercept=0) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_vline(xintercept=0, linetype="dashed") +
    coord_cartesian(ylim=c(-0.75,0.75), xlim=c(-0.75,0.75)) +
    labs(x="Pseudoblk score", y="Metacells score") +
    theme(
      axis.text = element_text(size=rel(0.75))
    )
  
  pdf(sprintf("%s/%s_comparison.pdf",args$outdir,i), width=7, height=5)
  print(p)
  dev.off()

  ####################################
  ## Compare distribution of scores ##
  ####################################
  
  to.plot <- rbind(
    data.table(value = virtual_chip_metacells.mtx[,i], class="metacells"),
    data.table(value = virtual_chip_pseudobulk.mtx[,i], class="pseudobulk")
  ) %>% .[value!=0]
  
  p <- gghistogram(to.plot, x="value", fill="class", bins=100) +
    theme(
      legend.title = element_blank(),
      axis.text = element_text(size=rel(0.75), color="black")
    )
  
  pdf(sprintf("%s/%s_histogram_scores.pdf",args$outdir,i), width=7, height=5)
  print(p)
  dev.off()
}

###################################################
## Plot score thresholds vs number binding sites ##
###################################################

virtual_chip_pseudobulk.mtx <- virtual_chip_pseudobulk.mtx[,opts$TFs]
virtual_chip_metacells.mtx <- virtual_chip_metacells.mtx[,opts$TFs]

seq.ranges <- seq(0.01,1,by=0.01)

to.plot <- seq.ranges %>% map(function(j) {
  data.table(
    tf = opts$TFs,
    min_score = j,
    log2_N = log2(colSums(virtual_chip_metacells.mtx>=j)+1)
  )
}) %>% rbindlist

p.lineplot <- ggplot(to.plot[log2_N>0], aes_string(x="min_score", y="log2_N", color="tf")) +
  geom_line(size=1) +
  labs(y="Number of predicted binding sites (log2)", x="Minimum in silico binding score") +
  scale_color_brewer(palette="Dark2") +
  theme_classic() +
  geom_vline(xintercept=args$min_chip_score,linetype="dashed",color="gray") +
  theme(
    axis.text = element_text(color="black"),
    legend.position = c(.9,.65),
    legend.title = element_blank()
  )

pdf(file.path(args$outdir,"lineplot_number_binding_sites_per_tf.pdf"), width=6, height=5)
print(p.lineplot)
dev.off()

