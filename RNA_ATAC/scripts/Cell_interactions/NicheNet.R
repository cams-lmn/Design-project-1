################
##            ##
##  NicheNet  ##
##            ##
################

source("/data/homes/louisc/Project_Babraham/RNA_ATAC/scripts/Settings.R")
library(nichenetr)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(annotate)
library(org.Hs.eg.db)
library(openxlsx)
detach("package:MOFA2")

######################
## Define arguments ##
######################

# p <- ArgumentParser(description='')
# p$add_argument('--DEG_overview',    type="character",    help='Differential analysis overview')
# p$add_argument('--TF',    type="character",    help='Transcription factor file')
# p$add_argument('--seed',            type="integer",     default=42,             help='Random seed')
# p$add_argument('--matrix',          type="character",  nargs="+",    help='Matrix to use')
# p$add_argument('--mofa_metadata',    type="character",    help='MOFA metadata')
# p$add_argument('--factors',    type="character",    help='MOFA factors')
# p$add_argument('--outdir',   type="character",    help='Output directory')
# args <- p$parse_args(commandArgs(TRUE))

# # START TEST
args <- list()
args$outdir <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/NicheNet"
args$sce <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/SingleCellExperiment_pseudobulk.rds"
args$cluster_ids <- 1:5
args$min_expr <- 4
args$DEG_overview <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/RNA/DE_res/DEG_overview.txt"
args$DE_gene_clusters <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/RNA/DE_res/genes_gene_clustering_ordered.txt"
args$threads <- 10
args$CellPhoneDB_data <- "/data/homes/louisc/Project_Babraham/pathways/CellPhoneDB_analysis/CellPhone_res_stat_LIGANDS_RECEPTORS.xlsx"
# # END TEST

#####################
## Define settings ##
#####################

dir.create(args$outdir, showWarnings = F, recursive = T)

args$groups <- combn(args$cluster_ids,2)

markers_TF <- NULL
for (i in 1:5){
  marker_tf_i <- read.table(sprintf("/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/RNA/marker_TFs/Marker_TF_cluster%s.txt",i),header=T,sep="\t")
  markers_TF <- rbind(markers_TF,cbind(marker_tf_i$Gene,rep(i,length(marker_tf_i$Gene))))
}
markers_TF <- data.frame(markers_TF)
colnames(markers_TF) <- c("TF","celltype")
markers_TF_undupl <- markers_TF[!duplicated(markers_TF$TF),]

###############
## Functions ##
###############

plot_ligand_parallel <- function(i){
    ligand <- ligands_oi[i]
    # print(ligand)
    ligands_all = ligand # this can be a list of multiple ligands if required
    targets_all = gsub("\\.","-",colnames(vis_ligand_target))[vis_ligand_target[gsub("\\.","-",rownames(vis_ligand_target))==ligand]!=0]

    active_signaling_network = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligands_all, targets_all = targets_all, weighted_networks = weighted_networks)

    if (nrow(active_signaling_network$sig)==0){
        print("No signaling info, copying gr info")
        active_signaling_network$sig <- active_signaling_network$gr
    }
    # For better visualization of edge weigths: normalize edge weights to make them comparable between signaling and gene regulatory interactions
    active_signaling_network_min_max = active_signaling_network
    active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
    active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)

    graph_min_max = diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_min_max, ligands_all = ligands_all, targets_all = targets_all, sig_color = "indianred", gr_color = "steelblue")
    return(graph_min_max)
}

##############
## Load sce ##
##############

sce <- readRDS(args$sce)

#########################
## Initialize NicheNet ##
#########################

if(!file.exists(paste0(args$outdir,"/lr_network_human.rds"))){
    options(timeout = 600)
    organism <- "human"

    lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
    lr_network <- lr_network %>% distinct(from, to)

    ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
    
    weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))

    ligand_tf_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_tf_matrix_nsga2r_final.rds"))

    sig_network <- readRDS(url("https://zenodo.org/record/7074291/files/signaling_network_human_21122021.rds"))

    gr_network <- readRDS(url("https://zenodo.org/record/7074291/files/gr_network_human_21122021.rds"))

    saveRDS(lr_network,file=paste0(args$outdir,"/lr_network_human.rds"))
    saveRDS(ligand_target_matrix,file=paste0(args$outdir,"/ligand_target_matrix.rds"))
    saveRDS(weighted_networks,file=paste0(args$outdir,"/weighted_networks.rds"))
    saveRDS(ligand_tf_matrix,file=paste0(args$outdir,"/ligand_tf_matrix.rds"))
    saveRDS(sig_network,file=paste0(args$outdir,"/sig_network.rds"))
    saveRDS(gr_network,file=paste0(args$outdir,"/gr_network.rds"))
} else {
    lr_network <- readRDS(paste0(args$outdir,"/lr_network_human.rds"))
    ligand_target_matrix <- readRDS(paste0(args$outdir,"/ligand_target_matrix.rds"))
    weighted_networks <- readRDS(paste0(args$outdir,"/weighted_networks.rds"))
    ligand_tf_matrix <- readRDS(paste0(args$outdir,"/ligand_tf_matrix.rds"))
    sig_network <- readRDS(paste0(args$outdir,"/sig_network.rds"))
    gr_network <- readRDS(paste0(args$outdir,"/gr_network.rds"))
}

print("NicheNet network")
print(lr_network)
print("NicheNet ligand target matrix")
print(ligand_target_matrix[1:5,1:5])

######################
## Define gene sets ##
#####################

## Load necessary objects  
############################

DEG_overview <- read.table(file=args$DEG_overview, header=T, sep="\t", quote="")

genes_clusters <- read.table(args$DE_gene_clusters)
colnames(genes_clusters) <- "Gene"

ix_genes_clusters <- sort(genes_clusters$Gene,index.return=T)$ix
DE_gene_clusters <- merge(genes_clusters,DEG_overview)[,c("Gene","gene_cluster")][sort(ix_genes_clusters,index.return=T)$ix,]
dim(DE_gene_clusters)
head(DE_gene_clusters)

## Marker genes
############################

list_Markergenes <- list(list(DEG_overview$Gene[DEG_overview$clust1_spec & (!is.na(DEG_overview$clust1_spec))],1),
                     list(DEG_overview$Gene[DEG_overview$clust2_spec & (!is.na(DEG_overview$clust2_spec))],2),
                     list(DEG_overview$Gene[DEG_overview$clust3_spec & (!is.na(DEG_overview$clust3_spec))],3),
                     list(DEG_overview$Gene[DEG_overview$clust4_spec & (!is.na(DEG_overview$clust4_spec))],4),
                     list(DEG_overview$Gene[DEG_overview$clust5_spec & (!is.na(DEG_overview$clust5_spec))],5),
                     list(markers_TF$TF[markers_TF$celltype=="1"],1),
                     list(markers_TF$TF[markers_TF$celltype=="2"],2),
                     list(markers_TF$TF[markers_TF$celltype=="3"],3),
                     list(markers_TF$TF[markers_TF$celltype=="4"],4),
                     list(markers_TF$TF[markers_TF$celltype=="5"],5),
                     list(DE_gene_clusters$Gene[DE_gene_clusters$gene_cluster=="A"],1),
                     list(DE_gene_clusters$Gene[DE_gene_clusters$gene_cluster=="B"],1:2),
                     list(DE_gene_clusters$Gene[DE_gene_clusters$gene_cluster=="C"],1:3),
                     list(DE_gene_clusters$Gene[DE_gene_clusters$gene_cluster=="D"],3),
                     list(DE_gene_clusters$Gene[DE_gene_clusters$gene_cluster=="E"],3:5),
                     list(DE_gene_clusters$Gene[DE_gene_clusters$gene_cluster=="F"],4:5),
                     list(DE_gene_clusters$Gene[DE_gene_clusters$gene_cluster=="G"],c(1,4,5)))

names(list_Markergenes) <- c("markerGenesCluster1","markerGenesCluster2","markerGenesCluster3","markerGenesCluster4","markerGenesCluster5",
                         "markerTFsCluster1","markerTFsCluster2","markerTFsCluster3","markerTFsCluster4","markerTFsCluster5",
                         "geneClusterA","geneClusterB","geneClusterC","geneClusterD","geneClusterE","geneClusterF","geneClusterG")

list_Markergenes[["genesCluster1"]] <- list(c(list_Markergenes[["geneClusterA"]][[1]],list_Markergenes[["geneClusterB"]][[2]]),1)
list_Markergenes[["genesCluster2"]] <- list(c(list_Markergenes[["geneClusterC"]][[1]],list_Markergenes[["geneClusterD"]][[2]]),2)
list_Markergenes[["genesCluster3"]] <- list(c(list_Markergenes[["geneClusterD"]][[1]],list_Markergenes[["geneClusterE"]][[2]]),3)
list_Markergenes[["genesCluster4"]] <- list(c(list_Markergenes[["geneClusterE"]][[1]],list_Markergenes[["geneClusterF"]][[2]]),4)
list_Markergenes[["genesCluster5"]] <- list(c(list_Markergenes[["geneClusterE"]][[1]],list_Markergenes[["geneClusterF"]][[2]]),5)

## DE genes
############################

c_names <- colnames(DEG_overview)
DEG_overview_bis <- DEG_overview

for (i in 1:5){
    cols_logFC <- which(grepl(i,c_names) & grepl("vs",c_names) & grepl("logFC",c_names))
    cols_FDR <-  which(grepl(i,c_names) & grepl("vs",c_names) & grepl("FDR",c_names))
    sign <- c(1,-1)[c(substr(gsub("[A-Z,a-z]+_","",c_names[cols_logFC]),1,1)==i)+1]

    DEG_overview_bis <- cbind(DEG_overview_bis,
                             (rowSums((t(t(DEG_overview_bis[,cols_logFC])*sign)>1) & (DEG_overview_bis[,cols_FDR]<0.01))>=1))
    c_names <- c(c_names,paste0("DE_clust",i))
}
colnames(DEG_overview_bis) <- c_names
for (i in 43:47){
    DEG_overview_bis[is.na(DEG_overview_bis[,i]),i] <- FALSE
}
head(DEG_overview_bis[rowSums(DEG_overview_bis[,43:46])>0,]) 

# Check total of DE genes
sum(DEG_overview_bis$DE_clust1 | DEG_overview_bis$DE_clust2 | DEG_overview_bis$DE_clust3 | DEG_overview_bis$DE_clust4 | DEG_overview_bis$DE_clust5)

# Make list 
list_DEgenes <- list(list(DEG_overview$Gene[DEG_overview_bis$DE_clust1],1),
                     list(DEG_overview$Gene[DEG_overview_bis$DE_clust2],2),
                     list(DEG_overview$Gene[DEG_overview_bis$DE_clust3],3),
                     list(DEG_overview$Gene[DEG_overview_bis$DE_clust4],4),
                     list(DEG_overview$Gene[DEG_overview_bis$DE_clust5],5))

names(list_DEgenes) <- c("DEGenesCluster1","DEGenesCluster2","DEGenesCluster3","DEGenesCluster4","DEGenesCluster5")

## Load pathway info 
############################

# # KEGG
# pathways <- list(NOTCH="hsa04330",
#                  HEDGEHOG="hsa04340")
# for (i in 1:length(pathways)){
#     pathway_name <- names(pathways)[i]
#     pathway_info <- read.table(sprintf("/data/homes/louisc/Project_Babraham/pathways/KEGG_pathways/%s.txt",pathways[[pathway_name]]), sep='\n')$V1
#     gene_line <- grep("GENE",pathway_info)
#     non_indent_lines <- which(substr(pathway_info,1,1)!=" ")
#     gene_info <- pathway_info[gene_line:(non_indent_lines[which(non_indent_lines==gene_line)+1]-1)]
#     pathways[[pathway_name]] <- gsub(";.+$","",gsub("^.+ [0-9][0-9][0-9]+  ","",gene_info))
# }
# list_pathways <- list(list(pathways[["NOTCH"]],1:5),
#                       list(pathways[["HEDGEHOG"]],1:5))
# names(list_pathways) <- c("NOTCH","HEDGEHOG")

# Reactome
tab_reactome <- read.csv("/data/homes/louisc/Project_Babraham/pathways/Reactome_pathways/interactions_summary_reactome_names_MR.csv",header=T)

list_reactome_pw <- c(tab_reactome$Reactome,tab_reactome$Reactome.,tab_reactome$Maria.s_suggestion)
list_reactome_pw <- unique(unlist(strsplit(list_reactome_pw,"/")))
list_reactome_pw <- list_reactome_pw[grepl("R-HSA",list_reactome_pw)]

reactome_pathways <- readRDS("/data/homes/louisc/Project_Babraham/pathways/Reactome_pathways/Reactome_pathways.rds")
reactome_pathway_names <- readRDS("/data/homes/louisc/Project_Babraham/pathways/Reactome_pathways/Reactome_pathway_names.rds")

list_pathways <- list()
names_list_pathways <- NULL
for (pw in list_reactome_pw){
    if (pw%in%names(reactome_pathways)){
        list_pathways <- c(list_pathways,list(list(getSYMBOL(reactome_pathways[[pw]], data='org.Hs.eg'),1:5)))
        names_list_pathways <- c(names_list_pathways,gsub(" ","",reactome_pathway_names[[pw]]))
    } else {
        print(pw)
    }
}
names(list_pathways) <- names_list_pathways
length(list_pathways)
length(list_pathways[[1]])

# CellPhoneDB
dat_cpdb_ligand <- read.xlsx(args$CellPhoneDB_data,2)
colnames(dat_cpdb_ligand) <- c("cluster","goi")
dat_cpdb_receptor <- read.xlsx(args$CellPhoneDB_data,3)
colnames(dat_cpdb_receptor) <- c("cluster","goi")

dat_cpdb <- rbind(dat_cpdb_ligand,dat_cpdb_receptor)
dat_cpdb <- data.frame(dat_cpdb)
dat_cpdb$cluster <- gsub("\\|[1-5]$","",dat_cpdb$cluster)
dim(dat_cpdb)
head(dat_cpdb)

sum(duplicated(paste(dat_cpdb$goi,dat_cpdb$cluster,sep="_")))
dat_cpdb[duplicated(paste(dat_cpdb$goi,dat_cpdb$cluster,sep="_")),]$goi
dat_cpdb <- dat_cpdb[!duplicated(paste(dat_cpdb$goi,dat_cpdb$cluster,sep="_")),]
dim(dat_cpdb)
head(dat_cpdb)

list_cpdb <- list(list(dat_cpdb$goi[dat_cpdb$cluster==1],1),
                  list(dat_cpdb$goi[dat_cpdb$cluster==2],2),
                  list(dat_cpdb$goi[dat_cpdb$cluster==3],3),
                  list(dat_cpdb$goi[dat_cpdb$cluster==4],4),
                  list(dat_cpdb$goi[dat_cpdb$cluster==5],5),
                  list(unique(dat_cpdb$goi),"X"))
names(list_cpdb) <- c("CellPhoneDBCluster1","CellPhoneDBCluster2","CellPhoneDBCluster3",
                      "CellPhoneDBCluster4","CellPhoneDBCluster5","CellPhoneDBClustersX")

# CellPhoneDB pathways

args$CellPhoneDB_data_full <- "/data/homes/louisc/Project_Babraham/pathways/CellPhoneDB_analysis/CellPhone_res_stat.xlsx"

dat_cpdb_full <- read.xlsx(args$CellPhoneDB_data_full,1)
 
cpdb_pw_list <- unique(dat_cpdb_full$classification[!is.na(dat_cpdb_full$classification)])

list_pathways_cpdb <- list()
for (pw in cpdb_pw_list){
    mask_pw <- dat_cpdb_full$classification==pw
    genes_pw <- unique(c(unlist(strsplit(dat_cpdb_full$gene_a[mask_pw],"_")),unlist(strsplit(dat_cpdb_full$gene_b[mask_pw],"_"))))
    genes_pw <- genes_pw[!is.na(genes_pw)]
    list_pathways_cpdb <- c(list_pathways_cpdb,list(list(genes_pw,1:5)))
}
names(list_pathways_cpdb) <- paste0(rep("CPDB",length(cpdb_pw_list)),gsub(" ","",cpdb_pw_list))
length(list_pathways_cpdb)
length(list_pathways_cpdb[[1]])


## Make gene sets  
####################  

genesets <- c(list_pathways,
              list_pathways_cpdb,
              list_Markergenes,
              list_DEgenes,
              list_cpdb)

# no underscores in geneset names please
sum(grepl("_",names(genesets)))
sum(grepl("/",names(genesets)))
names(genesets) <- gsub("/","-",names(genesets))
names(genesets) <- gsub("\\:","",names(genesets))
sum(grepl("/",names(genesets)))
sum(grepl(":",names(genesets)))

#######################
## NicheNet analysis ##
#######################


ligands_to_consider <- c("Expr","Expr","Pathway",
                         "Pathway","Pathway","Pathway",
                         "Pathway","Pathway","CellPhoneDB",
                         "CellPhoneDB","CellPhoneDB","Pathway",
                         "Pathway","Pathway","Pathway",
                         "CPDBPathway","CPDBPathway","CPDBPathway")#ltc
receptors_to_consider <- c("Expr","Expr","Expr",
                           "Expr","Pathway","Pathway",
                           "Pathway","Expr","CellPhoneDB",
                           "CellPhoneDB","CellPhoneDB","Pathway",
                           "Expr","Pathway","Expr",
                           "CPDBPathway","CPDBPathway","CPDBPathway") #rtc
targets_to_consider <- c("MarkerGene","MarkerTF","MarkerGene",
                         "MarkerTF","MarkerGene","MarkerTF",
                         "Pathway","Pathway","MarkerGene",
                         "MarkerTF","DEGene","DEinPathway",
                         "DEinPathway","DEGene","DEGene"
                         "MarkerGene","MarkerTF","DEGene") #ttc

ligands_to_consider <- c("CPDBPathway","CPDBPathway","CPDBPathway") #ltc
receptors_to_consider <- c("CPDBPathway","CPDBPathway","CPDBPathway") #rtc
targets_to_consider <- c("MarkerGene","MarkerTF","DEGene") #ttc

print(cbind(ligands_to_consider,receptors_to_consider,targets_to_consider))

for (i in 1:length(ligands_to_consider)){
    ltc <- ligands_to_consider[i]
    rtc <- receptors_to_consider[i]
    ttc <- targets_to_consider[i]

    if (ttc=="Expr"){
        stop("Not possible geneset of interest is equal to background")
    }

    model_performance_table <- NULL
    counter <- 1

    if (ltc == "Pathway"){
        name_sel_geneset <- names(genesets)
        name_sel_geneset <- name_sel_geneset[(!grepl("Cluster",names(genesets))) & (!grepl("CellPhoneDB",names(genesets))) & (!grepl("CPDB",names(genesets)))]
    } else if (ltc == "CellPhoneDB") {
        name_sel_geneset <- names(genesets)
        name_sel_geneset <- name_sel_geneset[grepl("CellPhoneDB",names(genesets))]
    } else if (ltc == "CPDBPathway") {
        name_sel_geneset <- names(genesets)
        name_sel_geneset <- name_sel_geneset[grepl("^CPDB",names(genesets))]
    } else {
        if (ttc == "MarkerTF"){
            name_sel_geneset <- names(genesets)
            name_sel_geneset <- name_sel_geneset[grepl("markerTF",name_sel_geneset)]
        } else if (ttc == "MarkerGene"){
            name_sel_geneset <- names(genesets)
            name_sel_geneset <- name_sel_geneset[grepl("markerGenes",name_sel_geneset)]
        }
    }

    dir.create(paste0(args$outdir,"/",ltc,"_ligands_",rtc,"_receptors_",ttc,"_targets"))
    for (name_geneset in name_sel_geneset){
        print(name_geneset)
        
        if (grepl("Pathway",ltc)){
            clusters_for_geneset <- c(1:5)
        } else if (ltc == "CellPhoneDB") {
            clusters_for_geneset <- genesets[[name_geneset]][[2]]
        } else {
            if (ttc == "MarkerTF"){
                clusters_for_geneset <- c(genesets[[name_geneset]][[2]],"X")
            } else if (ttc=="MarkerGene"){
                clusters_for_geneset <- c(genesets[[name_geneset]][[2]],"X")
            }             
        }

        print(sprintf("%s will be analysed for clusters: %s",name_geneset,paste(clusters_for_geneset,collapse=",")))

        for (cluster in clusters_for_geneset){
            print(sprintf("Analysis for %s in cluster %s",name_geneset,cluster))
            dir.create(paste0(args$outdir,"/",ltc,"_ligands_",rtc,"_receptors_",ttc,"_targets/",name_geneset,"_cluster",cluster))

            ## Expressed genes 
            #####################

            if (cluster=="X"){
                expressed_genes_sender <- rownames(sce)[rowMeans(logcounts(sce)[sce$celltype%in%c(1:5),])>=args$min_expr]
                expressed_genes_receiver <- rownames(sce)[rowMeans(logcounts(sce)[sce$celltype%in%c(1:5),])>=args$min_expr]
            } else {
                sce_sender <- sce[,sce$celltype %in% cluster]
                sce_receiver <- sce[,sce$celltype %in% cluster]
                expressed_genes_sender <- rownames(sce_sender)[logcounts(sce_sender)>=args$min_expr]
                expressed_genes_receiver <- rownames(sce_receiver)[logcounts(sce_receiver)>=args$min_expr]
            }

            ## Receiver & sender genes
            #######################

            if (ltc=="Expr"){
                sender_geneset_oi_raw <- expressed_genes_sender
                print("Sender set: all expressed genes")
            } else {
                sender_geneset_oi_raw <- genesets[[name_geneset]][[1]]
                print(paste0("Sender set: ",name_geneset))
            }

            if (ttc=="Expr"){
                receiver_geneset_oi_raw <- expressed_genes_receiver
                print("Receiver set: all expressed genes")
            } else {
                if (!(ltc=="Expr")){
                    if (ttc=="MarkerGene") {
                        if (cluster%in%c(1:5)){                        
                            receiver_geneset_oi_raw <- genesets[[paste0("markerGenesCluster",cluster)]][[1]]
                        } else {
                            receiver_geneset_oi_raw <- NULL
                            for (mgc in 1:5){
                                receiver_geneset_oi_raw <- c(receiver_geneset_oi_raw, 
                                                             genesets[[paste0("markerGenesCluster",mgc)]][[1]])
                            }
                            receiver_geneset_oi_raw <- unique(receiver_geneset_oi_raw)
                        }
                        print("Receiver set: marker genes")
                    } else if (ttc=="MarkerTF") {
                        if (cluster%in%c(1:5)){                        
                            receiver_geneset_oi_raw <- genesets[[paste0("markerTFsCluster",cluster)]][[1]]
                        } else {
                            receiver_geneset_oi_raw <- NULL
                            for (mgc in 1:5){
                                receiver_geneset_oi_raw <- c(receiver_geneset_oi_raw, 
                                                             genesets[[paste0("markerTFsCluster",mgc)]][[1]])
                            }
                            receiver_geneset_oi_raw <- unique(receiver_geneset_oi_raw)
                        }
                        print("Receiver set: marker TFs")
                    } else if (ttc=="DEGene") {
                        if (cluster%in%c(1:5)){                        
                            receiver_geneset_oi_raw <- genesets[[paste0("DEGenesCluster",cluster)]][[1]]
                        } else {
                            receiver_geneset_oi_raw <- NULL
                            for (mgc in 1:5){
                                receiver_geneset_oi_raw <- c(receiver_geneset_oi_raw, 
                                                             genesets[[paste0("DEGenesCluster",mgc)]][[1]])
                            }
                            receiver_geneset_oi_raw <- unique(receiver_geneset_oi_raw)
                        }
                        print("Receiver set: DE genes")
                    } else if (ttc=="DEinPathway") {
                        receiver_geneset_oi_raw <- genesets[[paste0("DEGenesCluster",cluster)]][[1]][genesets[[paste0("DEGenesCluster",cluster)]][[1]]%in%genesets[[name_geneset]][[1]]]
                        print(paste0("Receiver set: overlap of DE genes and ",name_geneset))
                    } else {
                        receiver_geneset_oi_raw <- genesets[[name_geneset]][[1]]
                        print(paste0("Receiver set: ",name_geneset))
                    }
                } else {
                    receiver_geneset_oi_raw <- genesets[[name_geneset]][[1]]
                    print(paste0("Receiver set: ",name_geneset))
                }
            }

            sender_geneset_oi_raw <- sender_geneset_oi_raw[sender_geneset_oi_raw%in%expressed_genes_sender]
            receiver_geneset_oi_raw <- receiver_geneset_oi_raw[receiver_geneset_oi_raw%in%expressed_genes_receiver]

            if(length(receiver_geneset_oi_raw)==0){
                write.table("No expressed genes in target set",col.names=F,row.names=F,quote=F,
                            file=paste0(args$outdir,"/",ltc,"_ligands_",rtc,"_receptors_",ttc,"_targets/",name_geneset,"_cluster",cluster,"/result_",name_geneset,"_cluster",cluster,".txt"))
                next()
            }

            # only consider genes also present in the NicheNet model - this excludes genes from the gene list for which the official HGNC symbol was not used by Puram et al.
            receiver_geneset_oi <- receiver_geneset_oi_raw %>% .[. %in% rownames(ligand_target_matrix)] 
            sender_geneset_oi <- sender_geneset_oi_raw %>% .[. %in% rownames(ligand_target_matrix)] # only consider genes also present in the NicheNet model - this excludes genes from the gene list for which the official HGNC symbol was not used by Puram et al.

            ## Determine genes of interest & background 
            ##############################################               

            background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
            head(background_expressed_genes)

            ## Potential ligands 
            #######################

            ligands = lr_network %>% pull(from) %>% unique()
            expressed_ligands = intersect(ligands,sender_geneset_oi)
         
            ## Potential Receptors
            #######################

            receptors = lr_network %>% pull(to) %>% unique()
            if (rtc=="Expr"){
                expressed_receptors = intersect(receptors,expressed_genes_receiver)
                print(paste0("Receptor set: all expressed genes"))
            } else if (rtc==ltc){
                expressed_receptors = intersect(receptors,sender_geneset_oi)
                print(paste0("Receptor set: ",name_geneset))
            } else {
                stop("currently not implemented")
            }

            ## Network
            #######################

            lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors)

            # Exception to few cconnections 
            if (nrow(lr_network_expressed)<2){
                write.table("Contrast skipped due to too few network connections",col.names=F,row.names=F,quote=F,
                            file=paste0(args$outdir,"/",ltc,"_ligands_",rtc,"_receptors_",ttc,"_targets/",name_geneset,"_cluster",cluster,"/result_",name_geneset,"_cluster",cluster,".txt"))
                next()
            }

            # print("Expressed network")
            # print(lr_network_expressed)

            potential_ligands = lr_network_expressed %>% pull(from) %>% unique()

            # head(potential_ligands)

            ## Ligand activity analysis
            ##############################

            ligand_activities = predict_ligand_activities(geneset = receiver_geneset_oi, background_expressed_genes = background_expressed_genes, 
                                                                    ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

            if (nrow(ligand_activities)>30){
                best_upstream_ligands = ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()
            } else {
                best_upstream_ligands = ligand_activities %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()
            }

            if(length(best_upstream_ligands)<2){
                write.table("Contrast skipped due to too few ligands",col.names=F,row.names=F,quote=F,
                            file=paste0(args$outdir,"/",ltc,"_ligands_",rtc,"_receptors_",ttc,"_targets/",name_geneset,"_cluster",cluster,"/result_",name_geneset,
                                        "_cluster",cluster,".txt"))
                next()
            }

            print(head(best_upstream_ligands))

            # Plot AUPR histogram
            p_hist_lig_activity = ggplot(ligand_activities, aes(x=aupr_corrected)) + 
            geom_histogram(color="black", fill="darkorange")  + 
            # geom_density(alpha=.1, fill="orange") + 
            geom_vline(xintercept=0,color="gray",linetype="dashed") +
            geom_vline(aes(xintercept=min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))), color="red", linetype="dashed", size=1) + 
            labs(x="ligand activity (PCC)", y = "# ligands") +
            theme_classic()

            pdf(file = paste0(args$outdir,"/",ltc,"_ligands_",rtc,"_receptors_",ttc,"_targets/",name_geneset,"_cluster",cluster,"/ligand_activity_hist_",name_geneset,
                              "_cluster",cluster,".pdf"),height = 6, width = 8)
            print( p_hist_lig_activity)
            dev.off()

            ## Ligand expression 
            #######################

            vis_ligand_expression <- logcounts(sce)[best_upstream_ligands,sce$celltype%in%args$cluster_ids]
            colnames(vis_ligand_expression) <- paste0(rep("cluster",ncol(vis_ligand_expression)),colnames(vis_ligand_expression))

            # Plot ligand expression
            color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
            p_ligand_expression = make_heatmap_ggplot(vis_ligand_expression, y_name= "Prioritized ligands",x_name="Cluster", color = color[100], legend_position = "top",
                                                      x_axis_position = "top", legend_title = "Expression\n(pseudobulk)") + theme(axis.text.y = element_text(face = "italic"))                                                                           
            pdf(file = paste0(args$outdir,"/",ltc,"_ligands_",rtc,"_receptors_",ttc,"_targets/",name_geneset,"_cluster",cluster,"/ligand_expression_",name_geneset,"_cluster",
                              cluster,".pdf"), height = 6, width = 12)
            print(p_ligand_expression)
            dev.off()

            ## Area under precision recall curve 
            #######################

            ligand_aupr_matrix = ligand_activities %>% dplyr::select(aupr_corrected) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

            vis_ligand_aupr = ligand_aupr_matrix[best_upstream_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("AUPR")
            p_ligand_aupr = vis_ligand_aupr %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", 
                                                                    x_axis_position = "top", legend_title = "AUPR\n(target gene prediction ability)")

            # Plot ligand activity
            pdf(file = paste0(args$outdir,"/",ltc,"_ligands_",rtc,"_receptors_",ttc,"_targets/",name_geneset,"_cluster",cluster,"/ligand_AUPR_",name_geneset,"_cluster",
                              cluster,".pdf"), height = 6, width = 12)
            print(p_ligand_aupr)
            dev.off()
            
            ## Plot ligand-target 
            ########################

            active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = receiver_geneset_oi, 
                                                                             ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
            active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, 
                                                                             cutoff = 0.25)

            order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
            order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
            rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
            colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

            if(length(order_ligands)<2){
                write.table("Contrast skipped due to too few ligand-target connections",col.names=F,row.names=F,quote=F,
                            file=paste0(args$outdir,"/",ltc,"_ligands_",rtc,"_receptors_",ttc,"_targets/",name_geneset,"_cluster",cluster,"/result_",name_geneset,"_cluster",
                                        cluster,".txt"))
                next()
            }

            if (is.null(nrow(active_ligand_target_links[order_targets,order_ligands]))){
                vis_ligand_target = matrix(active_ligand_target_links[order_targets,order_ligands],nrow=1)
                colnames(vis_ligand_target) <- order_ligands
                vis_ligand_target = vis_ligand_target %>% t()
                colnames(vis_ligand_target) <- order_targets
            } else {
                vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
            } 
            p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple", legend_position = "top",
                                                                                x_axis_position = "top",legend_title = "Regulatory potential")  + 
                                                                theme(axis.text.x = element_text(face = "italic")) + 
                                                                scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
            
            p_ligand_expression_target = make_heatmap_ggplot(vis_ligand_expression[gsub("\\.","-",order_ligands),], y_name= "Prioritized ligands",x_name="Cluster", 
                                                             color = color[100], legend_position = "top", x_axis_position = "top", 
                                                             legend_title = "Expression\n(pseudobulk)") + theme(axis.text.y = element_text(face = "italic"))   

            p_ligand_expression_target = make_heatmap_ggplot(vis_ligand_expression[gsub("\\.","-",order_ligands),], y_name= "Prioritized ligands", x_name="Cluster",
                                                             color = color[100], legend_position = "top", x_axis_position = "top", legend_title = "Expression\n(pseudobulk)") + 
                                                                theme(axis.text.y = element_text(face = "italic"))                                                                          

            p_ligand_aupr_target = vis_ligand_aupr[gsub("\\.","-",order_ligands),] %>% as.matrix() %>% 
                                        make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", 
                                                            x_axis_position = "top", legend_title = "AUPR\n(target gene\nprediction ability)")

            fig_ligand_target <- plot_grid(
                p_ligand_aupr_target + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
                p_ligand_expression_target + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
                p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                align = "h",
                nrow = 1,
                rel_widths = c(ncol(vis_ligand_aupr),ncol(vis_ligand_expression), ncol(vis_ligand_target)))

            legends_target = plot_grid(
                as_ggplot(get_legend(p_ligand_aupr_target)),
                as_ggplot(get_legend(p_ligand_expression_target)),
                as_ggplot(get_legend(p_ligand_target_network)),
                nrow = 1,
                align = "h")

            pdf(file = paste0(args$outdir,"/",ltc,"_ligands_",rtc,"_receptors_",ttc,"_targets/",name_geneset,"_cluster",cluster,"/ligand_target_gene_inference_",
                              name_geneset,"_cluster",cluster,".pdf"),
                height = 6, width = 12)
            print(plot_grid(fig_ligand_target, legends_target, rel_heights = c(10,2), nrow = 2, align = "hv"))
            dev.off()

            ## Ligand - receptor 
            ##########################

            lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
            best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

            lr_network_top_df_large = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

            lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
            lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

            dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
            hclust_ligands = hclust(dist_ligands, method = "ward.D2")
            order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

            order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

            if (nrow(lr_network_top_matrix)==1){
                vis_ligand_receptor_network <- t(as.matrix(lr_network_top_matrix[,order_ligands_receptor]))
                rownames(vis_ligand_receptor_network) <- rownames(lr_network_top_matrix)
            } else {
                dist_receptors = dist(lr_network_top_matrix, method = "binary")
                hclust_receptors = hclust(dist_receptors, method = "ward.D2")
                order_receptors = hclust_receptors$labels[hclust_receptors$order]

                order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))

                vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor] 
                rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
                colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
            }
            
            p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% 
                                            make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",
                                                                legend_title = "Prior interaction potential")

            p_ligand_expression_receptors = make_heatmap_ggplot(vis_ligand_expression[gsub("\\.","-", order_ligands_receptor),], y_name= "Prioritized ligands", x_name="Cluster",
                                                                color = color[100], legend_position = "top", x_axis_position = "top", legend_title = "Expression\n(pseudobulk)") + 
                                                                    theme(axis.text.y = element_text(face = "italic"))  

            p_ligand_aupr_receptors = vis_ligand_aupr[gsub("\\.","-",order_ligands_receptor),] %>% as.matrix() %>% 
                                        make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", 
                                                            x_axis_position = "top", legend_title = "AUPR\n(target gene\nprediction ability)")
                                                                    
            fig_ligand_receptor <- plot_grid(
                p_ligand_aupr_receptors + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
                p_ligand_expression_receptors + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
                p_ligand_receptor_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                align = "h",
                nrow = 1,
                rel_widths = c(col(vis_ligand_aupr),ncol(vis_ligand_expression), ncol(vis_ligand_receptor_network)))

            legends_receptor = plot_grid(
                as_ggplot(get_legend(p_ligand_aupr_receptors)),
                as_ggplot(get_legend(p_ligand_expression_receptors)),
                as_ggplot(get_legend(p_ligand_receptor_network)),
                nrow = 1,
                align = "h")

            pdf(file = paste0(args$outdir,"/",ltc,"_ligands_",rtc,"_receptors_",ttc,"_targets/",name_geneset,"_cluster",cluster,"/ligand_receptor_gene_inference_",
                              name_geneset,"_cluster",cluster,".pdf"), height = 6, width = 12)
            print(plot_grid(fig_ligand_receptor, legends_receptor, rel_heights = c(10,2), nrow = 2, align = "hv"))
            dev.off()

            ## Plot full analysis
            ##############################
            
            if (nrow(lr_network_top_matrix)==1){
                vis_ligand_receptor_network_targets <- t(as.matrix(lr_network_top_matrix[,order_ligands]))
                rownames(vis_ligand_receptor_network_targets) <- rownames(lr_network_top_matrix)
                vis_ligand_receptor_network_targets <- vis_ligand_receptor_network_targets %>% t()
            } else {
                vis_ligand_receptor_network_targets <- vis_ligand_receptor_network %>% t()
                vis_ligand_receptor_network_targets <- vis_ligand_receptor_network_targets[order_ligands,]
            }
            p_ligand_receptor_network_full = vis_ligand_receptor_network_targets %>% 
                                                make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", 
                                                                    x_axis_position = "top", legend_title = "Prior interaction potential")
        
            rel_widths <- c(ncol(vis_ligand_aupr)+2.5,ncol(vis_ligand_expression), ncol(vis_ligand_receptor_network_targets), ncol(vis_ligand_target)+2.5)
            rel_widths <- rel_widths/sum(rel_widths)
            rel_widths[rel_widths> 0.55] <- 0.55

            fig_ligand_full <- plot_grid(
                p_ligand_aupr_target + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
                p_ligand_expression_target + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
                p_ligand_receptor_network_full + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                align = "h",
                nrow = 1,
                rel_widths = rel_widths)

            legends_full = plot_grid(
                as_ggplot(get_legend(p_ligand_aupr_target)),
                as_ggplot(get_legend(p_ligand_expression_target)),
                as_ggplot(get_legend(p_ligand_receptor_network_full)),
                as_ggplot(get_legend(p_ligand_target_network)),
                nrow = 1,
                align = "h")

            pdf(file = paste0(args$outdir,"/",ltc,"_ligands_",rtc,"_receptors_",ttc,"_targets/",name_geneset,"_cluster",cluster,"/ligand_full_plot_",name_geneset,
                              "_cluster",cluster,".pdf"),
            height = 6, width = 18)
            print(plot_grid(fig_ligand_full, legends_full, rel_heights = c(10,2), nrow = 2, align = "hv"))
            dev.off()

            ## Plot ligand-target graph 
            ##############################

            dir.create(paste0(args$outdir,"/",ltc,"_ligands_",rtc,"_receptors_",ttc,"_targets/",name_geneset,"_cluster",cluster,"/loi"))

            print("Ligands of interest")
            ligands_oi <- gsub("\\.","-",rownames(vis_ligand_target))
            print(length(ligands_oi))    
            # ligands_oi <- "LAMA2"
            
            # res <- mclapply(1:length(ligands_oi),plot_ligand_parallel,mc.preschedule=T,mc.cores=args$threads)
            # names(res) <- ligands_oi
            # for (ligand in ligands_oi){
            #     res[[ligand]]$global_attrs$value[15] <- 4
            #     render_graph(res[[ligand]]) %>% export_svg() %>% charToRaw %>% rsvg_pdf(paste0(args$outdir,"/loi/",ligand,"_to_targets.pdf"))
            # }

            # As for loop for debugging
            # for (ligand in ligands_oi){
            #     print(ligand)
            #     ligands_all = ligand # this can be a list of multiple ligands if required
            #     targets_all = gsub("\\.","-",colnames(vis_ligand_target))[vis_ligand_target[gsub("\\.","-",rownames(vis_ligand_target))==ligand]!=0]

            #     active_signaling_network = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligands_all, targets_all = targets_all, weighted_networks = weighted_networks)

            #     if (nrow(active_signaling_network$sig)==0){
            #         print("No signaling info, copying gr info")
            #         active_signaling_network$sig <- active_signaling_network$gr
            #     }
            #     # For better visualization of edge weigths: normalize edge weights to make them comparable between signaling and gene regulatory interactions
            #     active_signaling_network_min_max = active_signaling_network
            #     active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
            #     active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)

            #     graph_min_max = diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_min_max, ligands_all = ligands_all, targets_all = targets_all, sig_color = "indianred", gr_color = "steelblue")

            #     render_graph(graph_min_max) %>% export_svg() %>% charToRaw %>% rsvg_pdf(paste0(args$outdir,"/cluster",cluster,"/loi/",ligand,"_to_targets.pdf"))
            # }

            ## Enrichment of set of interest genes in targets 
            ##########################

            # change rounds and folds here, to two rounds to reduce time: normally: do multiple rounds
            k_folds = 3 # 3-fold
            n_rounds = 2 # 2 rounds

            # # background for test
            # if (cluster=="X"){
            #     background_test <- genesets[[name_geneset]][[1]] %>% .[. %in% rownames(ligand_target_matrix)]
            # } else {
            #     background_test <- markers_TF$TF[!markers_TF$celltype==cluster] %>% .[. %in% rownames(ligand_target_matrix)]
            #     # background_test <- markers_TF_undupl$TF[!markers_TF_undupl$celltype==cluster] %>% .[. %in% rownames(ligand_target_matrix)]
            # }
            # background_test[!(background_test%in%receiver_geneset_oi)] 
            # print("Background test:")
            # print(length(background_test))
            background_test <- background_expressed_genes

            model_performance <- c(length(ligands_oi),length(background_test))
            
            # enrichment
            if (length(receiver_geneset_oi)==1){
                model_performance <- as.matrix(c(model_performance,rep(NA,6)))
            } else {
                geneset_oi_predictions_list = seq(n_rounds) %>% lapply(assess_rf_class_probabilities, folds = k_folds, geneset = receiver_geneset_oi, 
                                                                       background_expressed_genes = background_test,
                                                                       ligands_oi = best_upstream_ligands, ligand_target_matrix = ligand_target_matrix)

                target_prediction_performances_cv = geneset_oi_predictions_list %>% lapply(classification_evaluation_continuous_pred_wrapper) %>% 
                    bind_rows() %>% mutate(round=seq(1:nrow(.)))

                model_performance <- c(model_performance,target_prediction_performances_cv$auroc %>% mean(),
                                       target_prediction_performances_cv$aupr %>% mean(), target_prediction_performances_cv$pearson %>% mean()) 

                if (nrow(do.call(rbind,(geneset_oi_predictions_list %>% lapply(calculate_fraction_top_predicted, quantile_cutoff = 0.95))))<4){
                    model_performance <- as.matrix(c(model_performance,rep(NA,3)))
                } else {
                    target_prediction_performances_discrete_cv = geneset_oi_predictions_list %>% lapply(calculate_fraction_top_predicted, quantile_cutoff = 0.95) %>%  
                                        bind_rows() %>% ungroup() %>% mutate(round=rep(1:length(geneset_oi_predictions_list), each = n_rounds))

                    # Fraction of geneset oi genes that belong to the top 5% predicted targets
                    # print(target_prediction_performances_discrete_cv %>% filter(true_target) %>% .$fraction_positive_predicted %>% mean())

                    # Fraction of non geneset oi genes that belong to the top 5% predicted targets
                    # print(target_prediction_performances_discrete_cv %>% filter(!true_target) %>% .$fraction_positive_predicted %>% mean())

                    # Fisher test
                    target_prediction_performances_discrete_fisher = geneset_oi_predictions_list %>% lapply(calculate_fraction_top_predicted_fisher, 
                                                                                                            quantile_cutoff = 0.95) 
                    # print(target_prediction_performances_discrete_fisher %>% unlist() %>% mean())

                    model_performance <- as.matrix(c(model_performance,
                                            target_prediction_performances_discrete_cv %>% filter(true_target) %>% .$fraction_positive_predicted %>% mean(),
                                            target_prediction_performances_discrete_cv %>% filter(!true_target) %>% .$fraction_positive_predicted %>% mean(),
                                            target_prediction_performances_discrete_fisher %>% unlist() %>% mean()))
                }
            }
            rownames(model_performance) <- c("n_loi","n_bg","av_auroc","av_aupr","av_pearson","Frac_oi_top5%pred","Frac_noi_top5%pred","Fisher_pval")
            write.table(model_performance,file=paste0(args$outdir,"/",ltc,"_ligands_",rtc,"_receptors_",ttc,"_targets/",name_geneset,"_cluster",cluster,"/model_performance_",name_geneset,"_cluster",cluster,".txt"),
                        row.names=TRUE,col.names=FALSE,sep="\t",quote=FALSE)
                        
            model_performance_table <- rbind(model_performance_table,t(model_performance))
            rownames(model_performance_table)[counter] <- paste0(name_geneset,"_cluster",cluster)
            print(model_performance_table[counter,])
            counter <- counter + 1

            ## Make full table
            ##############################

            # Targets
            full_target_table <- ligand_activities$test_ligand %>% lapply(get_weighted_ligand_target_links,geneset = receiver_geneset_oi, ligand_target_matrix = ligand_target_matrix, n=Inf) %>% bind_rows() %>% drop_na()
            colnames(full_target_table)[1] <- "test_ligand"
            full_target_table <- merge(full_target_table,ligand_activities[,c(1,4)],all.x=TRUE)
            colnames(full_target_table)[1] <- "ligand"
            # print(dim(full_target_table))
            # print(head(full_target_table))
            write.table(full_target_table,file=paste0(args$outdir,"/",ltc,"_ligands_",rtc,"_receptors_",ttc,"_targets/",name_geneset,"_cluster",cluster,
                                                      "/ligand_target_table_",name_geneset,"_cluster",cluster,".txt"),col.names=T,row.names=F,quote=F,sep="\t")

            # Receptors
            lr_network_top <- lr_network %>% filter(from %in% ligand_activities$test_ligand & to %in% expressed_receptors) %>% distinct(from,to)
            best_upstream_receptors <- lr_network_top %>% pull(to) %>% unique()
            full_receptor_table <- weighted_networks$lr_sig %>% filter(from %in% ligand_activities$test_ligand & to %in% best_upstream_receptors)
            colnames(full_receptor_table)[1] <- "test_ligand"
            full_receptor_table <- merge(full_receptor_table,ligand_activities[,c(1,4)],all.x=TRUE)
            colnames(full_receptor_table)[1:2] <- c("ligand","receptor")
            # print(dim(full_receptor_table))
            # print(head(full_receptor_table))
            write.table(full_receptor_table,file=paste0(args$outdir,"/",ltc,"_ligands_",rtc,"_receptors_",ttc,"_targets/",name_geneset,"_cluster",cluster,
                                                        "/ligand_receptor_table_",name_geneset,"_cluster",cluster,".txt"),col.names=T,row.names=F,quote=F,sep="\t")
        }
    }
    print(model_performance_table)
    write.table(model_performance_table,file=paste0(args$outdir,"/",ltc,"_ligands_",rtc,"_receptors_",ttc,"_targets/model_performance_table_",ltc,"_ligands_",rtc,"_receptors_",ttc,"_targets.txt"),col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
}


# ###################
# ## Overlap table ##
# ###################

# overlap_table <- NULL
# counter <- 1
# for (name_geneset in names(genesets)){
#     print(name_geneset)
    
#     clusters_for_geneset <- c(1:5,"X")

#     print(sprintf("%s will be analysed for clusters: %s",name_geneset,paste(clusters_for_geneset,collapse=",")))

#     for (cluster in clusters_for_geneset){

#         ## Expressed genes 
#         #####################

#         if (cluster=="X"){
#             expressed_genes <- rownames(sce)[rowMeans(logcounts(sce)[sce$celltype%in%c(1:5),])>=args$min_expr]
#         } else {
#             sce_cluster <- sce[,sce$celltype %in% cluster]
#             expressed_genes <- rownames(sce_cluster)[logcounts(sce_cluster)>=args$min_expr]
#         }

#         ## Expressed genes 
#         #####################

#         pw_genes_expr <- unlist(genesets[[name_geneset]][[1]][genesets[[name_geneset]][[1]]%in%expressed_genes])
#         if (cluster=="X"){
#             marker_genes_expr <- markers_TF_undupl$TF[markers_TF_undupl$TF%in%expressed_genes]
#         } else {
#             marker_genes_expr <- markers_TF$TF[(markers_TF$celltype==cluster) & (markers_TF$TF%in%expressed_genes)]
#         }

#         overlap <- c(sum(marker_genes_expr%in%pw_genes_expr),length(marker_genes_expr),
#                          length(pw_genes_expr),length(expressed_genes))
        
#         chisq_tab <- rbind(c(overlap[1],overlap[2]-overlap[1]),
#                            c(overlap[3],overlap[4]-overlap[3]))
        
#         overlap <- c(overlap, chisq_tab[,1]/chisq_tab[,2],
#                      c(chisq_tab[,1]/chisq_tab[,2])[1]/c(chisq_tab[,1]/chisq_tab[,2])[2],
#                      chisq.test(chisq_tab)$p.value)

#         names(overlap) <- c("n_mark_in_pw","n_mark","n_pw","n_expr","rate_mark_in_pw","rate_pw_expr","enrichment","pval")

#         overlap_table <- rbind(overlap_table,overlap)
#         rownames(overlap_table)[counter] <- paste0(name_geneset,"_cluster",cluster)
#         #print(overlap_table[counter,])
#         counter <- counter + 1
#     }
# }
# print(head(overlap_table))
# write.table(overlap_table,file=paste0(args$outdir,"/overlap_table.txt"),col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")

# ######################
# ## Prediction table ##
# ######################

# prediction_eval_table <- NULL
# counter <- 1
# for (name_geneset in names(genesets)){
#     print(name_geneset)
    
#     clusters_for_geneset <- c(1:5,"X")

#     print(sprintf("%s will be analysed for clusters: %s",name_geneset,paste(clusters_for_geneset,collapse=",")))

#     for (cluster in clusters_for_geneset){

#         ## Expressed genes 
#         #####################

#         if (cluster=="X"){
#             expressed_genes <- rownames(sce)[rowMeans(logcounts(sce)[sce$celltype%in%c(1:5),])>=args$min_expr]
#         } else {
#             sce_cluster <- sce[,sce$celltype %in% cluster]
#             expressed_genes <- rownames(sce_cluster)[logcounts(sce_cluster)>=args$min_expr]
#         }

#         ## Expressed genes 
#         #####################

#         pw_genes_expr <- unlist(genesets[[name_geneset]][[1]][genesets[[name_geneset]][[1]]%in%expressed_genes])
#         if (cluster=="X"){
#             marker_genes_expr <- markers_TF_undupl$TF[markers_TF_undupl$TF%in%expressed_genes]
#         } else {
#             marker_genes_expr <- markers_TF$TF[(markers_TF$celltype==cluster) & (markers_TF$TF%in%expressed_genes)]
#         }

#         k_folds = 3 # 3-fold
#         n_rounds = 2 # 2 rounds

#         lr_network_expressed = lr_network %>% filter(from %in% pw_genes_expr & to %in% expressed_genes)
#         potential_ligands <- lr_network_expressed %>% pull(from) %>% unique()
#         if(length(potential_ligands)<2){
#             next()
#         }

#         marker_genes_expr <- marker_genes_expr %>% .[. %in% rownames(ligand_target_matrix)]
#         background_expressed_genes <- expressed_genes %>% .[. %in% rownames(ligand_target_matrix)]

#         geneset_oi_predictions_list = seq(n_rounds) %>% lapply(assess_rf_class_probabilities, folds = k_folds, geneset = marker_genes_expr, 
#                                                                 background_expressed_genes = background_expressed_genes,
#                                                                 ligands_oi = potential_ligands, ligand_target_matrix = ligand_target_matrix)

#         target_prediction_performances_discrete_cv = geneset_oi_predictions_list %>% lapply(calculate_fraction_top_predicted, quantile_cutoff = 0.95) %>%  
#                                         bind_rows() %>% ungroup() %>% mutate(round=rep(1:length(geneset_oi_predictions_list), each = n_rounds))

#         target_prediction_performances_discrete_fisher = geneset_oi_predictions_list %>% lapply(calculate_fraction_top_predicted_fisher, 
#                                                                             quantile_cutoff = 0.95) 

#         prediction_eval <- c(length(potential_ligands),
#                              length(marker_genes_expr),
#                              length(background_expressed_genes),
#                              target_prediction_performances_discrete_cv %>% filter(true_target) %>% .$fraction_positive_predicted %>% mean(),
#                              target_prediction_performances_discrete_cv %>% filter(!true_target) %>% .$fraction_positive_predicted %>% mean(),
#                              target_prediction_performances_discrete_fisher %>% unlist() %>% mean())

#         names(prediction_eval) <- c("n_lig","n_mark","n_bg","Frac_oi_top5%pred","Frac_noi_top5%pred","pval")

#         prediction_eval_table <- rbind(prediction_eval_table,prediction_eval)
#         rownames(prediction_eval_table)[counter] <- paste0(name_geneset,"_cluster",cluster)
#         #print(prediction_eval_table[counter,])
#         counter <- counter + 1
#     }
# }
# print(head(prediction_eval_table))

# #########################
# ## NicheNet on command ##
# #########################

# name_command <- "Ntbinddom_cluster_1"
# cluster <- 1
# name_pathway <- "Nucleotide-bindingdomain,leucinerichrepeatcontainingreceptor(NLR)signalingpathways"
# receiver_geneset_oi_raw <- genesets[[name_pathway]][[1]]
# receiver_geneset_oi_raw <- receiver_geneset_oi_raw[receiver_geneset_oi_raw%in%markers_TF$TF[(markers_TF$celltype==cluster)]]
# sender_geneset_oi_raw <-  genesets[[name_pathway]][[1]]

# dir.create(paste0(args$outdir,"/on_command"))
# dir.create(paste0(args$outdir,"/on_command/",name_command))

# ## Expressed genes 
# #####################

# if (cluster=="X"){
#     expressed_genes_sender <- rownames(sce)[rowMeans(logcounts(sce)[sce$celltype%in%c(1:5),])>=args$min_expr]
#     expressed_genes_receiver <- rownames(sce)[rowMeans(logcounts(sce)[sce$celltype%in%c(1:5),])>=args$min_expr]
# } else {
#     sce_sender <- sce[,sce$celltype %in% cluster]
#     sce_receiver <- sce[,sce$celltype %in% cluster]
#     expressed_genes_sender <- rownames(sce_sender)[logcounts(sce_sender)>=args$min_expr]
#     expressed_genes_receiver <- rownames(sce_receiver)[logcounts(sce_receiver)>=args$min_expr]
# }

# ## Receiver & sender genes
# #######################

# sender_geneset_oi_raw <- sender_geneset_oi_raw[sender_geneset_oi_raw%in%expressed_genes_sender]
# receiver_geneset_oi_raw <- receiver_geneset_oi_raw[receiver_geneset_oi_raw%in%expressed_genes_receiver]

# if(length(receiver_geneset_oi_raw)==0){
#     write.table("No expressed genes in target set",col.names=F,row.names=F,quote=F,
#                 file=paste0(args$outdir,"/on_command/",name_command,"/result_",name_geneset,"_cluster",cluster,".txt"))
#     } else {

#     # only consider genes also present in the NicheNet model - this excludes genes from the gene list for which the official HGNC symbol was not used by Puram et al.
#     receiver_geneset_oi <- receiver_geneset_oi_raw %>% .[. %in% rownames(ligand_target_matrix)] 
#     sender_geneset_oi <- sender_geneset_oi_raw %>% .[. %in% rownames(ligand_target_matrix)] # only consider genes also present in the NicheNet model - this excludes genes from the gene list for which the official HGNC symbol was not used by Puram et al.

#     ## Determine genes of interest & background 
#     ##############################################               

#     # print("Geneset of interest: sender")
#     # print(length(sender_geneset_oi))
#     # print(head(sender_geneset_oi))
#     # print("Geneset of interest: receiver")
#     # print(length(receiver_geneset_oi))
#     # print(head(receiver_geneset_oi))

#     background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
#     head(background_expressed_genes)

#     ## Potential ligands 
#     #######################

#     ligands = lr_network %>% pull(from) %>% unique()
#     expressed_ligands = intersect(ligands,sender_geneset_oi)

#     ## Potential Receptors
#     #######################

#     expressed_receptors = intersect(receptors,sender_geneset_oi)
#     # expressed_receptors = intersect(receptors,expressed_genes_receiver)

#     ## Network
#     #######################

#     lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors)

#     # Exception to few cconnections 
#     if (nrow(lr_network_expressed)<2){
#         write.table("Contrast skipped due to too few network connections",col.names=F,row.names=F,quote=F,
#                     paste0(args$outdir,"/on_command/",name_command,"/result_",name_command,"_cluster",cluster,".txt"))
#     } else {

#         # print("Expressed network")
#         # print(lr_network_expressed)

#         potential_ligands = lr_network_expressed %>% pull(from) %>% unique()

#         # head(potential_ligands)

#         ## Ligand activity analysis
#         ##############################

#         ligand_activities = predict_ligand_activities(geneset = receiver_geneset_oi, background_expressed_genes = background_expressed_genes, 
#                                                                 ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

#         if (nrow(ligand_activities)>30){
#             best_upstream_ligands = ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()
#         } else {
#             best_upstream_ligands = ligand_activities %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()
#         }

#         if(length(best_upstream_ligands)<2){
#             write.table("Contrast skipped due to too few ligands",col.names=F,row.names=F,quote=F,
#                         file=paste0(args$outdir,"/on_command/",name_command,"/result_",name_command,
#                                     "_cluster",cluster,".txt"))
#         } else {

#             print(head(best_upstream_ligands))

#             # Plot AUPR histogram
#             p_hist_lig_activity = ggplot(ligand_activities, aes(x=aupr_corrected)) + 
#             geom_histogram(color="black", fill="darkorange")  + 
#             # geom_density(alpha=.1, fill="orange") + 
#             geom_vline(xintercept=0,color="gray",linetype="dashed") +
#             geom_vline(aes(xintercept=min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))), color="red", linetype="dashed", size=1) + 
#             labs(x="ligand activity (PCC)", y = "# ligands") +
#             theme_classic()

#             pdf(file = paste0(args$outdir,"/on_command/",name_command,"/ligand_activity_hist_",name_command,
#                                 "_cluster",cluster,".pdf"),height = 6, width = 8)
#             print( p_hist_lig_activity)
#             dev.off()

#             ## Ligand expression 
#             #######################

#             vis_ligand_expression <- logcounts(sce)[best_upstream_ligands,sce$celltype%in%args$cluster_ids]
#             colnames(vis_ligand_expression) <- paste0(rep("cluster",ncol(vis_ligand_expression)),colnames(vis_ligand_expression))

#             # Plot ligand expression
#             color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
#             p_ligand_expression = make_heatmap_ggplot(vis_ligand_expression, y_name= "Prioritized ligands",x_name="Cluster", color = color[100], legend_position = "top",
#                                                         x_axis_position = "top", legend_title = "Expression\n(pseudobulk)") + theme(axis.text.y = element_text(face = "italic"))                                                                           
#             pdf(file = paste0(args$outdir,"/on_command/",name_command,"/ligand_expression_",name_command,"_cluster",
#                                 cluster,".pdf"), height = 6, width = 12)
#             print(p_ligand_expression)
#             dev.off()

#             ## Area under precision recall curve 
#             #######################

#             ligand_aupr_matrix = ligand_activities %>% dplyr::select(aupr_corrected) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

#             vis_ligand_aupr = ligand_aupr_matrix[best_upstream_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("AUPR")
#             p_ligand_aupr = vis_ligand_aupr %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", 
#                                                                     x_axis_position = "top", legend_title = "AUPR\n(target gene prediction ability)")

#             # Plot ligand activity
#             pdf(file = paste0(args$outdir,"/on_command/",name_command,"/ligand_AUPR_",name_command,"_cluster",
#                                 cluster,".pdf"), height = 6, width = 12)
#             print(p_ligand_aupr)
#             dev.off()

#             ## Plot ligand-target 
#             ########################

#             active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = receiver_geneset_oi, 
#                                                                                 ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
#             active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, 
#                                                                                 cutoff = 0.25)

#             order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
#             order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
#             rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
#             colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

#             if(length(order_ligands)<2){
#                 write.table("Contrast skipped due to too few ligand-target connections",col.names=F,row.names=F,quote=F,
#                             file=paste0(args$outdir,"/on_command/",name_command,"/result_",name_command,"_cluster",
#                                         cluster,".txt"))
#             } else {

#                 if (is.null(nrow(active_ligand_target_links[order_targets,order_ligands]))){
#                     vis_ligand_target = matrix(active_ligand_target_links[order_targets,order_ligands],nrow=1)
#                     colnames(vis_ligand_target) <- order_ligands
#                     vis_ligand_target = vis_ligand_target %>% t()
#                     colnames(vis_ligand_target) <- order_targets
#                 } else {
#                     vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
#                 } 
#                 p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple", legend_position = "top",
#                                                                                     x_axis_position = "top",legend_title = "Regulatory potential")  + 
#                                                                     theme(axis.text.x = element_text(face = "italic")) + 
#                                                                     scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))

#                 p_ligand_expression_target = make_heatmap_ggplot(vis_ligand_expression[gsub("\\.","-",order_ligands),], y_name= "Prioritized ligands",x_name="Cluster", 
#                                                                     color = color[100], legend_position = "top", x_axis_position = "top", 
#                                                                     legend_title = "Expression\n(pseudobulk)") + theme(axis.text.y = element_text(face = "italic"))   

#                 p_ligand_expression_target = make_heatmap_ggplot(vis_ligand_expression[gsub("\\.","-",order_ligands),], y_name= "Prioritized ligands", x_name="Cluster",
#                                                                     color = color[100], legend_position = "top", x_axis_position = "top", legend_title = "Expression\n(pseudobulk)") + 
#                                                                     theme(axis.text.y = element_text(face = "italic"))                                                                          

#                 p_ligand_aupr_target = vis_ligand_aupr[gsub("\\.","-",order_ligands),] %>% as.matrix() %>% 
#                                             make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", 
#                                                                 x_axis_position = "top", legend_title = "AUPR\n(target gene\nprediction ability)")

#                 fig_ligand_target <- plot_grid(
#                     p_ligand_aupr_target + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
#                     p_ligand_expression_target + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
#                     p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
#                     align = "h",
#                     nrow = 1,
#                     rel_widths = c(ncol(vis_ligand_aupr),ncol(vis_ligand_expression), ncol(vis_ligand_target)))

#                 legends_target = plot_grid(
#                     as_ggplot(get_legend(p_ligand_aupr_target)),
#                     as_ggplot(get_legend(p_ligand_expression_target)),
#                     as_ggplot(get_legend(p_ligand_target_network)),
#                     nrow = 1,
#                     align = "h")

#                 pdf(file = paste0(args$outdir,"/on_command/",name_command,"/ligand_target_gene_inference_",
#                                     name_command,"_cluster",cluster,".pdf"),
#                     height = 6, width = 12)
#                 print(plot_grid(fig_ligand_target, legends_target, rel_heights = c(10,2), nrow = 2, align = "hv"))
#                 dev.off()

#                 ## Ligand - receptor 
#                 ##########################

#                 lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
#                 best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

#                 lr_network_top_df_large = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

#                 lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
#                 lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

#                 dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
#                 hclust_ligands = hclust(dist_ligands, method = "ward.D2")
#                 order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

#                 order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

#                 if (nrow(lr_network_top_matrix)==1){
#                     vis_ligand_receptor_network <- t(as.matrix(lr_network_top_matrix[,order_ligands_receptor]))
#                     rownames(vis_ligand_receptor_network) <- rownames(lr_network_top_matrix)
#                 } else {
#                     dist_receptors = dist(lr_network_top_matrix, method = "binary")
#                     hclust_receptors = hclust(dist_receptors, method = "ward.D2")
#                     order_receptors = hclust_receptors$labels[hclust_receptors$order]

#                     order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))

#                     vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor] 
#                     rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
#                     colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
#                 }

#                 p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% 
#                                                 make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",
#                                                                     legend_title = "Prior interaction potential")

#                 p_ligand_expression_receptors = make_heatmap_ggplot(vis_ligand_expression[gsub("\\.","-", order_ligands_receptor),], y_name= "Prioritized ligands", x_name="Cluster",
#                                                                     color = color[100], legend_position = "top", x_axis_position = "top", legend_title = "Expression\n(pseudobulk)") + 
#                                                                         theme(axis.text.y = element_text(face = "italic"))  

#                 p_ligand_aupr_receptors = vis_ligand_aupr[gsub("\\.","-",order_ligands_receptor),] %>% as.matrix() %>% 
#                                             make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", 
#                                                                 x_axis_position = "top", legend_title = "AUPR\n(target gene\nprediction ability)")
                                                                        
#                 fig_ligand_receptor <- plot_grid(
#                     p_ligand_aupr_receptors + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
#                     p_ligand_expression_receptors + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
#                     p_ligand_receptor_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
#                     align = "h",
#                     nrow = 1,
#                     rel_widths = c(col(vis_ligand_aupr),ncol(vis_ligand_expression), ncol(vis_ligand_receptor_network)))

#                 legends_receptor = plot_grid(
#                     as_ggplot(get_legend(p_ligand_aupr_receptors)),
#                     as_ggplot(get_legend(p_ligand_expression_receptors)),
#                     as_ggplot(get_legend(p_ligand_receptor_network)),
#                     nrow = 1,
#                     align = "h")

#                 pdf(file = paste0(args$outdir,"/on_command/",name_command,"/ligand_receptor_gene_inference_",
#                                     name_command,"_cluster",cluster,".pdf"), height = 6, width = 12)
#                 print(plot_grid(fig_ligand_receptor, legends_receptor, rel_heights = c(10,2), nrow = 2, align = "hv"))
#                 dev.off()

#                 ## Plot full analysis
#                 ##############################

#                 if (nrow(lr_network_top_matrix)==1){
#                     vis_ligand_receptor_network_targets <- t(as.matrix(lr_network_top_matrix[,order_ligands]))
#                     rownames(vis_ligand_receptor_network_targets) <- rownames(lr_network_top_matrix)
#                     vis_ligand_receptor_network_targets <- vis_ligand_receptor_network_targets %>% t()
#                 } else {
#                     vis_ligand_receptor_network_targets <- vis_ligand_receptor_network %>% t()
#                     vis_ligand_receptor_network_targets <- vis_ligand_receptor_network_targets[order_ligands,]
#                 }
#                 p_ligand_receptor_network_full = vis_ligand_receptor_network_targets %>% 
#                                                     make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", 
#                                                                         x_axis_position = "top", legend_title = "Prior interaction potential")

#                 rel_widths <- c(ncol(vis_ligand_aupr)+2.5,ncol(vis_ligand_expression), ncol(vis_ligand_receptor_network_targets), ncol(vis_ligand_target)+2.5)
#                 rel_widths <- rel_widths/sum(rel_widths)
#                 rel_widths[rel_widths> 0.55] <- 0.55

#                 fig_ligand_full <- plot_grid(
#                     p_ligand_aupr_target + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
#                     p_ligand_expression_target + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
#                     p_ligand_receptor_network_full + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
#                     p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
#                     align = "h",
#                     nrow = 1,
#                     rel_widths = rel_widths)

#                 legends_full = plot_grid(
#                     as_ggplot(get_legend(p_ligand_aupr_target)),
#                     as_ggplot(get_legend(p_ligand_expression_target)),
#                     as_ggplot(get_legend(p_ligand_receptor_network_full)),
#                     as_ggplot(get_legend(p_ligand_target_network)),
#                     nrow = 1,
#                     align = "h")

#                 pdf(file = paste0(args$outdir,"/on_command/",name_command,"/ligand_full_plot_",name_command,
#                                     "_cluster",cluster,".pdf"),
#                 height = 6, width = 18)
#                 print(plot_grid(fig_ligand_full, legends_full, rel_heights = c(10,2), nrow = 2, align = "hv"))
#                 dev.off()
#             }
#         }
#     }
# }
