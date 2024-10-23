###########################
##                       ##
##  GSA_visualisation.R  ##
##                       ##
###########################

source("/data/homes/louisc/Project_Babraham/RNA_ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--GSA_res',          type="character",    help='Gene set analysis results')
p$add_argument('--outdir',          type="character",    help='Output directory')
p$add_argument('--matrix',          type="character",  nargs="+",    help='Matrix to use')
args <- p$parse_args(commandArgs(TRUE))

args <- list()
args$matrix <- "RNA"
args$soi_type <- "Gene"
args$GSA_res <- sprintf("/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/%s/GSA/%s_enrichResults_%s_Clusters_Maria.txt",args$matrix,args$matrix,args$soi_type)
args$outdir <- sprintf("/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/%s/GSA/",args$matrix)

####################
## Initialisation ##
####################

######################
## Load GSA results ##
######################

GSA_results <- read.table(args$GSA_res,sep="\t",header=T)

# Pardon my Dutch
GSA_results$SHOW_IN_FIGURE <- GSA_results$SHOW_IN_FIGURE=='WAAR'

######################
## Load GSA results ##
######################

to_plot <- GSA_results[GSA_results$SHOW_IN_FIGURE,-c(which(colnames(GSA_results)%in%c("link","overlapId","userId")))]

# take into account duplicated descriptions for different gene set classes
duplicated_descriptions <- NULL
duplicated_descriptions <- unique(paste(to_plot$description,to_plot[,colnames(to_plot)==paste0(args$soi_type,"_cluster")])[which(duplicated(paste(to_plot$description,to_plot[,colnames(to_plot)==paste0(args$soi_type,"_cluster")])))])
if (length(duplicated_descriptions)>0){
    for (i in duplicated_descriptions){
        to_plot$description[paste(to_plot$description,to_plot[,colnames(to_plot)==paste0(args$soi_type,"_cluster")])==i] <- paste0(to_plot$description[paste(to_plot$description,to_plot[,colnames(to_plot)==paste0(args$soi_type,"_cluster")])==i],
                                                                                                                                   " (",to_plot$GOClass[paste(to_plot$description,to_plot[,colnames(to_plot)==paste0(args$soi_type,"_cluster")])==i],
                                                                                                                                   ")")
    }
}

to_plot$GOClass <- factor(to_plot$GOClass,levels=c("MF","BP","CC","KEGG"))
to_plot$order_column <- to_plot$enrichmentRatio+as.numeric(to_plot$GOClass)*100

to_plot$description <- paste(to_plot$description,to_plot[,colnames(to_plot)==paste0(args$soi_type,"_cluster")])
to_plot$description <- factor(to_plot$description,levels=unique(to_plot$description[sort(to_plot$order_column,decreasing=T,index.return=T)$ix]))

p <- ggplot(data=to_plot, aes(x=enrichmentRatio,y=description,fill=GOClass)) +
        geom_bar(stat="identity") + 
        facet_wrap(as.formula(paste0("~ ",args$soi_type,"_cluster")),scales="free_y") + 
        scale_y_discrete(labels=gsub(" .$","",as.character(to_plot$description)))

pdf(sprintf("%s%s_enrichResults_%s_Clusters.pdf",args$outdir,args$matrix,args$soi_type),height = 6, width = 16)
print(p)
# barplot(GSA_results$enrichmentRatio,names.arg=GSA_results$description,horiz=TRUE)
dev.off()