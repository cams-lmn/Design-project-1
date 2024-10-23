####################################
##                                ##
##  Plot exprs vs access scatter  ##
##                                ##
####################################

source("/data/homes/louisc/Project_Babraham/RNA_ATAC/scripts/Settings.R")
library("viridis")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce_metacell',    type="character",    help='Expression data on metacell level')
p$add_argument('--atac_metacell_Distal',    type="character",    help='Accessibility data on metacell level (Distal)')
p$add_argument('--goi_gme',          type="character",  nargs="+",    help='Genes of interest: general markers epiblast')
p$add_argument('--goi_hk',          type="character",  nargs="+",    help='Genes of interest: house keeping genes')
p$add_argument('--goi_nc_regex',          type="character",  nargs="+",    help='Regex for genes of interest: negative control')
p$add_argument('--outdir',   type="character",    help='Output directory')
p$add_argument('--markers_TF',  type="character", help='TF marker file')
args <- p$parse_args(commandArgs(TRUE))

# START TEST
args <- list()
args$sce_metacell <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/metacells/SingleCellExperiment_metacells_nodiff.rds"
args$atac_metacell_Distal <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/metacells/GeneScoreMatrix_distal_summarized_experiment_metacells_nodiff.rds"
args$goi_gme <- c("POU5F1", "NANOG", "TDGF1", "SALL4", "GDF3", "KHDC3L")
args$goi_hk <- c("ACTB", "LDHA", "GAPDH", "RPL19", "HMBS", "TBP", "UBC", "HPRT1", "PGK1")
args$goi_nc <- NULL
args$goi_nc_regex <- "^OR[1-9,A-Z]+$"
args$goi_hm <- c("SPIC", "DPPA5", "DPPA3", "KLF4", "DNMT3L", "ARGFX", "TFCP2L1", "PRDM14", "FBP1", "KLF17", "KLF5", "GDF3", "FGF4", "NODAL", "DLL3", "GBX2", "ETV4", "PODXL", 
                 "SFRP2", "ETV5", "FGF2", "HES1", "SALL4", "ZIC2", "SOX11", "SALL2", "FST", "MYC", "SPP1", "FZD7", "TCF15", "CDH2", "TCF7L1", "SALL1")
args$outdir <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/rna_vs_acc/scatter"
args$pseudotime <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/trajectories/N2P/N2P_trajectory_nodiff.txt.gz"
args$seed <- 42
# END TEST

print(args)

####################
## Initialisation ##
####################

if(!dir.exists(args$outdir)){dir.create(args$outdir)}

celltypes.to.plot <- 1:5


###########################
## function scatter plot ##
###########################

create_scatter_plots <- function(goi,name){
    if(!dir.exists(sprintf("%s/%s",args$outdir,name))){dir.create(sprintf("%s/%s",args$outdir,name))}

    for (i in goi){
        print(i)
        to.plot_i <- data.frame(metacell=colnames(sce),
                            RNA=assay(sce,"minmax")[i,],
                            ATAC=assay(atac.sce.Distal,"minmax")[i,],
                            cluster=as.factor(cluster_sce))
        #print(head(to.plot_i))
        
        g <- ggplot(data=to.plot_i,aes(x=ATAC,y=RNA,fill=cluster)) +
            labs(x=sprintf("%s accessibility",i), y=sprintf("%s expression",i))  + 
            geom_hline(yintercept=0.5, linetype=2, color="gray") +
            geom_vline(xintercept=0.5, linetype=2, color="gray") +
            geom_point(size=2,shape = 21) +
            ylim(0,1) + xlim(0,1) + 
            scale_fill_manual(values=opts$celltype.colors) +
            theme_classic()
        
        output_plot(g,sprintf("%s/%s/%s_scatter",args$outdir,name,i),width=5,height=5)
    }
}

##########################
## Load RNA & ATAC data ##
##########################

print("RNA & ATAC data")

sce <- readRDS(args$sce_metacell)
atac.sce.Distal <- readRDS(args$atac_metacell_Distal)

print(sce)
print(atac.sce.Distal)

assay(atac.sce.Distal,"logcounts") <- log2(1e6*(sweep(assay(atac.sce.Distal),2,colSums(assay(atac.sce.Distal)),"/"))+1)

if (sum(colnames(assay(sce,"logcounts"))==colnames(assay(atac.sce.Distal,"logcounts")))!=ncol(sce)){
  stop("RNAseq & ATACseq not metacell info not in same order")
}

########################
## Minmax expr & ATAC ##
########################

jpeg(sprintf("%s/hist_exprs_q95_gw.jpg",args$outdir))
hist(apply(assay(sce,"logcounts"),1,quantile,probs=0.95,na.rm=T))
dev.off()

jpeg(sprintf("%s/hist_acces_q95_gw.jpg",args$outdir))
hist(apply(assay(atac.sce.Distal,"logcounts"),1,quantile,probs=0.95,na.rm=T))
dev.off()

jpeg(sprintf("%s/hist_exprs_q5_gw.jpg",args$outdir))
hist(apply(assay(sce,"logcounts"),1,quantile,probs=0.05,na.rm=T))
dev.off()

jpeg(sprintf("%s/hist_acces_q5_gw.jpg",args$outdir))
hist(apply(assay(atac.sce.Distal,"logcounts"),1,quantile,probs=0.05,na.rm=T))
dev.off()

########################
## Minmax expr & ATAC ##
########################

# assay(sce,"minmax") <- minmax.normalisation(logcounts(sce))
# assay(atac.sce.Distal,"minmax") <- minmax.normalisation(assay(atac.sce.Distal,"logcounts"))

# assay(sce,"minmax") <- minmax.10pct.normalisation(logcounts(sce))
# assay(atac.sce.Distal,"minmax") <- minmax.10pct.normalisation(assay(atac.sce.Distal,"logcounts"))

assay(sce,"minmax") <- minmax.fb.normalisation(logcounts(sce),lb=0,ub=10)
assay(atac.sce.Distal,"minmax") <- minmax.fb.normalisation(assay(atac.sce.Distal,"logcounts"),lb=1,ub=8)

#####################
## Load marker TFs ##
#####################

# TF markers with clusters
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
## Load Pseudotime ##
#####################

print("Pseudotime")

pseudotime.dt <- fread(args$pseudotime)
pseudotime.dt <- data.frame(pseudotime.dt)

pseudotime.dt$PT <- pseudotime.dt$PC1+pseudotime.dt$DC1*100
pseudotime.dt <- pseudotime.dt %>% setorder(PT)

print(dim(pseudotime.dt))
print(str(pseudotime.dt))
print(head(pseudotime.dt))

rownames(pseudotime.dt) <- pseudotime.dt$cell

sum(pseudotime.dt[colnames(sce),]$cell==colnames(sce))

cluster_sce <- pseudotime.dt[colnames(sce),]$cluster

##############################
## Scatter plots: marker TF ##
##############################

print("Marker TF")
name <- "markerTF"
create_scatter_plots(markers_TF_undupl$TF,name)

##################################
## Scatter plots: genes heatmap ##
##################################

print("Genes heatmap")
name <- "genes_heatmap"
create_scatter_plots(args$goi_hm,name)

#############################################
## Scatter plots: general markers epiblast ##
#############################################

print("General markers epiblast")
name <- "general_markers_epiblast"
create_scatter_plots(args$goi_gme,name)

#######################################
## Scatter plots: housekeeping genes ##
#######################################

print("House keeping genes")
name <- "housekeeping_genes"
create_scatter_plots(args$goi_hk,name)

#####################################
## Scatter plots: negative control ##
#####################################

print("Negative control")
name <- "negative_control"

set.seed(args$seed)
args$goi_nc <- rownames(sce)[grepl(args$goi_nc_regex,rownames(sce))]
length(args$goi_nc)
args$goi_nc <- args$goi_nc[args$goi_nc%in%rownames(atac.sce.Distal)]
length(args$goi_nc)

if(length(args$goi_nc)>20){
    args$goi_nc <- sample(args$goi_nc,20)
}

create_scatter_plots(args$goi_nc,name)


