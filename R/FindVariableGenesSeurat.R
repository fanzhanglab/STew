require(Seurat)

FindVariableGenesSeurat <- function(data, x.low.cutoff = 0.1, x.high.cutoff = 8,
                                     y.cutoff = 1, y.high.cutoff = Inf, num.bin = 0,
                                     binning.method = "equal_width", sort.results = TRUE,
                                     display.progress = TRUE, ...)
{

  genes.use <- rownames(data)
  if (!inherits(data, "dgCMatrix")) {
    data <- as(as.matrix(data), "dgCMatrix")
  }
  ## (1) get means and variances
  gene.mean <- FastExpMean(data, display.progress)
  names(gene.mean) <- genes.use
  gene.dispersion <- FastLogVMR(data, display.progress)
  names(gene.dispersion) <- genes.use

  gene.dispersion[is.na(x = gene.dispersion)] <- 0
  gene.mean[is.na(x = gene.mean)] <- 0

  mv.df <- data.frame(gene.mean, gene.dispersion)
  rownames(mv.df) <- rownames(data)

  ## (OPTIONAL) do the binning correction
  if (num.bin > 0) {
    if (binning.method == "equal_width") {
      data_x_bin <- cut(x = gene.mean, breaks = num.bin)
    }
    else if (binning.method == "equal_frequency") {
      data_x_bin <- cut(x = gene.mean, breaks = c(-1, quantile(gene.mean[gene.mean >
                                                                           0], probs = seq(0, 1, length.out = num.bin))))
    }
    else {
      stop(paste0("Invalid selection: '", binning.method,
                  "' for 'binning.method'."))
    }
    names(x = data_x_bin) <- names(x = gene.mean)
    mean_y <- tapply(X = gene.dispersion, INDEX = data_x_bin,
                     FUN = mean)
    sd_y <- tapply(X = gene.dispersion, INDEX = data_x_bin,
                   FUN = sd)
    gene.dispersion.scaled <- (gene.dispersion - mean_y[as.numeric(x = data_x_bin)])/sd_y[as.numeric(x = data_x_bin)]
    gene.dispersion.scaled[is.na(x = gene.dispersion.scaled)] <- 0
    ##names(gene.dispersion.scaled) <- names(gene.mean)

    mv.df$gene.dispersion.scaled <- gene.dispersion.scaled
  }

  return(mv.df)
}


