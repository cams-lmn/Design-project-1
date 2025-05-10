##########################
##                      ##
##  Settings.R - RNAseq ##
##                      ##
##########################


####################
## Load libraries ##
####################

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))
# suppressPackageStartupMessages(library(scds))     # problem not in bioconducter instead next
suppressPackageStartupMessages(library(scDblFinder))
suppressPackageStartupMessages(library(reticulate))
suppressPackageStartupMessages(library(scuttle))
# suppressPackageStartupMessages(library(gganimate))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(data.table))
# suppressPackageStartupMessages(library(Rclusterpp))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(SingleR))

#############
# Functions #
#############

# load_SingleCellExperiment
#############

load_SingleCellExperiment <- function(file, normalise = FALSE, features = NULL, cells = NULL, remove_non_expressed_genes = FALSE) {
  library(SingleCellExperiment); library(scran); library(scater);
  sce <- readRDS(file)
  if (!is.null(cells)) sce <- sce[,cells]
  if (!is.null(features)) sce <- sce[features,]
  if (remove_non_expressed_genes) sce <- sce[which(Matrix::rowSums(counts(sce))>15),]
  if (normalise) sce <- logNormCounts(sce)
  
  return(sce)
}

# RegressOutMatrix
#############

RegressOutMatrix <- function(mtx, covariates = NULL, features_idx = NULL, verbose = TRUE) {

  # Check features_idx
  if (is.null(features_idx)) {
    features_idx <- 1:nrow(mtx)
  }
  if (is.character(features_idx)) {
    features_idx <- intersect(features_idx, rownames(mtx))
    if (length(features_idx) == 0) {
      stop("Cannot use features that are beyond the scope of mtx")
    }
  } else if (max(features_idx) > nrow(mtx)) {
    stop("Cannot use features that are beyond the scope of mtx")
  }

  # Check dataset dimensions
  if (nrow(covariates) != ncol(mtx)) {
    stop("Uneven number of cells between latent data and expression data")
  }

  # Subset
  mtx <- mtx[features_idx,]
  mtx.dimnames <- dimnames(mtx)

  # Create formula for regression
  vars.to.regress <- colnames(covariates)
  fmla <- paste('GENE ~', paste(vars.to.regress, collapse = '+')) %>% as.formula

  # In this code, we'll repeatedly regress different Y against the same X
  # (covariates) in order to calculate residuals.  Rather that repeatedly
  # call lm to do this, we'll avoid recalculating the QR decomposition for the
  # covariates matrix each time by reusing it after calculating it once
  regression.mat <- cbind(covariates, mtx[1,])
  colnames(regression.mat) <- c(colnames(covariates), "GENE")
  qr <- lm(fmla, data = regression.mat, qr = TRUE)$qr
  rm(regression.mat)

  # Make results matrix
  data.resid <- matrix(
    nrow = nrow(mtx),
    ncol = ncol(mtx)
  )

  if (verbose) pb <- txtProgressBar(char = '=', style = 3, file = stderr())

  # Extract residuals from each feature by using the pre-computed QR decomposition
  for (i in 1:length(features_idx)) {
    regression.mat <- cbind(covariates, mtx[features_idx[i], ])
    colnames(regression.mat) <- c(vars.to.regress, 'GENE')
    regression.mat <- qr.resid(qr = qr, y = mtx[features_idx[i],])  # The function qr.resid returns the residuals when fitting y to the matrix with QR decomposition.
    data.resid[i, ] <- regression.mat
    if (verbose) {
      setTxtProgressBar(pb = pb, value = i / length(features_idx))
    }
  }

  if (verbose) close(con = pb)

  dimnames(data.resid) <- mtx.dimnames

  return(data.resid)
}

# ggplot_theme_NoAxes
#############

ggplot_theme_NoAxes <- function() {
  theme(
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
}

# Pseudobulk
#############

pseudobulk_sce_fn <- function(x, assay = NULL, by, fun = c("sum", "mean", "median"), scale = FALSE) {

  # check validity of input arguments
  fun <- match.arg(fun)
  if (is.null(assay))  assay <- assayNames(x)[1]

  # store aggregation parameters &
  # nb. of cells that went into aggregation
  md <- metadata(x)
  md$agg_pars <- list(assay = assay, by = by, fun = fun, scale = scale)

  # get aggregation function
  # fun <- switch(fun, sum = "rowSums", mean = "rowMeans", median = "rowMedians")

  # drop missing factor levels & tabulate number of cells
  cd <- dplyr::mutate_if(as.data.frame(colData(x)), is.factor, droplevels)
  colData(x) <- DataFrame(cd, row.names = colnames(x),check.names = FALSE)
  md$n_cells <- table(as.data.frame(colData(x)[, by]))

  # assure 'by' colData columns are factors so that missing combinations aren't dropped
  for (i in by)
    if (!is.factor(x[[i]]))
      x[[i]] <- factor(x[[i]])

  # split cells & compute pseudo-bulks
  cs <- .split_cells(x, by)
  # pb <- .pb(x, cs, assay, fun)
  pb <- .pb(x=x, by=by, fun=fun)
  if (scale & length(by) == 2) {
    ls <- lapply(.pb(x, cs, "counts", "rowSums"), colSums)
    pb <- lapply(seq_along(pb), function(i) pb[[i]] / 1e6 * ls[[i]])
    names(pb) <- names(ls)
  }

  # construct SCE
  pb <- SingleCellExperiment(pb, metadata = md)

  # propagate 'colData' columns that are unique across 2nd 'by'
  if (length(by) == 2) {
    cd <- colData(x)
    ids <- colnames(pb)
    counts <- vapply(ids, function(u) {
      m <- as.logical(match(cd[, by[2]], u, nomatch = 0))
      vapply(cd[m, ], function(u) length(unique(u)), numeric(1))
    }, numeric(ncol(colData(x))))
    cd_keep <- apply(counts, 1, function(u) all(u == 1))
    cd_keep <- setdiff(names(which(cd_keep)), by)
    if (length(cd_keep) != 0) {
      m <- match(ids, cd[, by[2]], nomatch = 0)
      cd <- cd[m, cd_keep, drop = FALSE]
      rownames(cd) <- ids
      colData(pb) <- cd
    }
  }
  return(pb)
}

.split_cells <- function(x, by) {
  if (is(x, "SingleCellExperiment")) x <- colData(x)
  cd <- data.frame(x[by], check.names = FALSE)
  cd <- data.table(cd, cell = rownames(x)) %>% split(by = by, sorted = TRUE, flatten = FALSE)
  purrr::map_depth(cd, length(by), "cell")
}

.pb <- function(x, by, fun) {

  # compute pseudobulks
  # y <- scuttle::summarizeAssayByGroup(x, assay.type = assay, ids = (ids <- colData(x)[by]), statistics = fun, BPPARAM = BiocParallel::SerialParam())
  y <- scuttle::summarizeAssayByGroup(x, ids = colData(x)[by], statistics = fun)
  colnames(y) <- y[[by[length(by)]]]

  if (length(by) == 1)  return(assay(y))

  # reformat into one assay per 'by[1]'
  if (is.factor(ids <- y[[by[1]]]))
    ids <- droplevels(ids)
  is <- split(seq_len(ncol(y)), ids)
  ys <- map(is, ~assay(y)[, .])

  # fill in missing combinations
  for (i in seq_along(ys)) {
    fill <- setdiff(unique(y[[by[2]]]), colnames(ys[[i]]))
    if (length(fill != 0)) {
      foo <- matrix(0, nrow(x), length(fill))
      colnames(foo) <- fill
      foo <- cbind(ys[[i]], foo)
      o <- paste(sort(unique(y[[by[2]]])))
      ys[[i]] <- foo[, o]
    }
  }
  return(ys)
}

# z_score
############
z_score <- function(x){
  (x - mean(x)) / sd(x)
}


# output_plot
############

output_plot <- function(p,fn,width,height,UMAP=FALSE){
  axes_line_size <- 1.5

  pdf(paste0(fn,".pdf"), width=width, height=height)
  print(p)
  dev.off()

  png(paste0(fn,".png"), width=width, height=height, units="in", res=300)
  print(p)
  dev.off()
  
  if (UMAP){
    p_nolabs <- p +
      theme(
        legend.position = "none",
        plot.title=element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y =  element_blank(),
        axis.title.x =  element_blank(),
        axis.line=element_blank()
      ) 
      width <- 7
      height <- 7
  } else {
    p_nolabs <- p +
      theme(
        legend.position = "none",
        plot.title=element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y =  element_blank(),
        axis.title.x =  element_blank(),
        axis.line=element_line(size=axes_line_size),
        axis.ticks=element_line(size=axes_line_size,colour="black"),
        axis.ticks.length=unit(.25, "cm"),
      )  
  }

  pdf(paste0(fn,"_nolabs.pdf"), width=width, height=height)
  print(p_nolabs)
  dev.off()

  png(paste0(fn,"_nolabs.png"), width=width, height=height, units="in", res=300)
  print(p_nolabs)
  dev.off()
}


###############
# Set options # # Needs to be changes
###############

opts <- list()
opts$color_scheme <- c("#A06932","#DBB216","#EFE32A","#D7DB54","#A7B019","#9AA126",
                       "#7E7721","#5D6821","#353D8A","#CA3639","#DCADCD","#7A327E")

opts$celltype.colors <- c("#DAF7A6","#FFC300","#BCA043","#C70039","#900C3F","#EBECF0","#FFFFFF")

# opts$naieve_pluri <- c("KLF4","KLF5","TFCP2L1","DNMT3L","FGF4","KLF17")
# opts$general_pluri <- c("NANOG","POU5F1","SALL4","SOX2","TDGF1")
# opts$common_post_impl <- c("ETV4","ETV5","MYC","SOX11","FZD7","CDH2","SALL2","SFRP2","ZIC2")
# opts$housekeeping <- c("GADPH","ACTB","TBP","B2M","HPRT1","HMBS","PGK1","PPIA","GUSB","TFRC","YWHAZ","SDHA","LDHA","UBC")