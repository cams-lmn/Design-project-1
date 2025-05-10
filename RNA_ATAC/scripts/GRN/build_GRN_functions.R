##########
## conv ##
##########

myfun_conv <- function(i, cutoff = 0.84) {
  # Check if the gene is in 'genes'
  if (!i %in% genes) {
    warning(paste("Gene", i, "not found in 'genes'"))
    return(NULL)
  }

  # Print progress every 10%
  if (which(genes == i) %% round(length(genes) / 10) == 0) {
    print(paste0(round(which(genes == i) / length(genes) * 100), "% done!"))
  }
  
  # Check if 'tmp' is a data frame and has 'gene' column
  if (!is.data.frame(tmp) || !"gene" %in% colnames(tmp)) {
    stop("tmp must be a data frame with a 'gene' column")
  }
  
  # Get TFs for the gene 'i'
  tfs <- tmp[tmp$gene == i, tf]
  
  dat_return <- NULL
  plot_counter <- 0
  for (j in tfs) {
    x <- rna_tf.mtx[j, ]
    y <- rna_targets.mtx[i, ]

    tf_expression_rate <- sum(x > 0) / length(x)
    gene_expression_rate <- sum(y > 0) / length(y)

    pid_file <- sprintf("%s/x_y_%d.txt", args$outdir, Sys.getpid())

    if (tf_expression_rate >= cutoff && gene_expression_rate >= cutoff) {
      # Check if x and y are matrices and ensure proper row extraction
      if (is.matrix(rna_tf.mtx)) {
        x <- rna_tf.mtx[j, ]
      }
      if (is.matrix(rna_targets.mtx)) {
        y <- rna_targets.mtx[i, ]
      }
    
      if (args$plot_correlations && plot_counter < 10) {
        to.plot <- data.frame(
          tf = rep(j, length(rna_tf.mtx[j, ])),
          gene = rep(i, length(rna_tf.mtx[j, ])),
          x = rna_tf.mtx[j, ], y = rna_targets.mtx[i, ],
          sample = factor(gsub("#.+$", "", names(rna_targets.mtx[i, ])), 
                  levels = args$samples)
        )
        
        p1 <- ggscatter(to.plot, x = "x", y = "y", fill = "sample", size = 5, shape = 21, 
                        add = "reg.line", add.params = list(color = "black", fill = "lightgray"), conf.int = TRUE) +
          stat_cor(method = "pearson", label.x.npc = 0.5, label.y.npc = 0.05) +
          stat_regline_equation(label.x.npc = 0.5, label.y.npc = 0.1) +
          scale_fill_manual(values = opts$color_scheme[1:length(unique(to.plot$sample))]) +
          labs(x = sprintf("%s expression", j), y = sprintf("%s expression", i)) +
          theme_classic() +
          theme(axis.text = element_text(size = rel(0.7)))
        
        pdf(sprintf("%s/plots_cor_score%s/%s_vs_%s.pdf", args$outdir, args$min_chip_score, j, i), 
            width = 7, height = 5)
        suppressMessages(print(p1))
        dev.off()
        plot_counter <- plot_counter + 1
      }
    
      # Linear regression model
      lm.fit <- lm(y ~ x)
      
      if (is.na(coef(lm.fit)[2])) {
        dat_return <- rbind(dat_return, c(j, i, NA, NA, paste(1:6, collapse = "-")))
      } else {
        dat_return <- rbind(dat_return, c(j, i, 
                                          round(coef(lm.fit)[[2]], 3), 
                                          format(summary(lm.fit)$coefficients[2, 4], digits = 3),
                                          paste(clusters, collapse = "-")))
      }
      last <- nrow(dat_return)
      cat(tf_expression_rate, gene_expression_rate, dat_return[last, 3], dat_return[last, 4],"\n", file = pid_file, sep = "\t", append = TRUE)
    }
  }

  # If no data was collected, return NULL
  if (is.null(dat_return) || nrow(dat_return) == 0) {
    return(NULL)
  }
  
  colnames(dat_return) <- c("tf", "gene", "beta", "pvalue", "clusters_for_calc")
  dat_return <- data.frame(dat_return, stringsAsFactors = F)
  return(dat_return)
}