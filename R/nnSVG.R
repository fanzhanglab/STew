#Only needs BioCParallel and BRISC

nnSVG <- function(input, spatial_coords = NULL, X = NULL, 
                  assay_name = "logcounts", 
                  n_neighbors = 10, order = "AMMD", 
                  n_threads = 1, BPPARAM = NULL, 
                  verbose = FALSE) {
  
  if (is(input, "SpatialExperiment")) {
    spe <- input
    stopifnot(assay_name %in% assayNames(spe))
    # check for rows or columns of all zero counts
    if ("counts" %in% assayNames(spe)) {
      if (sum(rowSums(counts(spe)) == 0) > 0 | sum(colSums(counts(spe)) == 0) > 0) {
        warning("Rows (genes) and/or columns (spots) containing all zero counts ", 
                "have been found. Please see examples in tutorial for code to ", 
                "filter out zeros and/or low-expressed genes to avoid errors.")
      }
    }
  }
  
  if (!is.null(X)) {
    stopifnot(nrow(X) == ncol(input))
  }
  
  if (is.null(BPPARAM)) {
    BPPARAM <- MulticoreParam(workers = n_threads)
  }
  
  # -----------------------
  # run BRISC for each gene
  # -----------------------
  
  if (is(input, "MulticoreParam")) {
    y <- assays(spe)[[assay_name]]
    coords <- spatialCoords(spe)
  } else {
    y <- input
    coords <- spatial_coords
    row_names <- rownames(input)
  }
  
  # scale coordinates proportionally
  range_all <- max(apply(coords, 2, function(col) diff(range(col))))
  coords <- apply(coords, 2, function(col) (col - min(col)) / range_all)
  
  # calculate ordering of coordinates
  order_brisc <- BRISC_order(coords, order = order, verbose = verbose)
  
  # calculate nearest neighbors
  nn_brisc <- BRISC_neighbor(coords, n.neighbors = n_neighbors, n_omp = 1, 
                             search.type = "tree", ordering = order_brisc, 
                             verbose = verbose)
  
  # run BRISC using parallelization
  ix <- seq_len(nrow(y))
  out_brisc <- bplapply(ix, function(i) {
    # fit model (intercept-only model if x is NULL)
    y_i <- y[i, ]
    suppressWarnings({
      runtime <- system.time({
        out_i <- BRISC_estimation(coords = coords, y = y_i, x = X, 
                                  cov.model = "exponential", 
                                  ordering = order_brisc, neighbor = nn_brisc, 
                                  verbose = verbose)
      })
    })
    res_i <- c(
      out_i$Theta, 
      loglik = out_i$log_likelihood, 
      runtime = runtime[["elapsed"]]
    )
    res_i
  }, BPPARAM = BPPARAM)
  
  # collapse output list into matrix
  mat_brisc <- do.call("rbind", out_brisc)
  
  # --------------------
  # calculate statistics
  # --------------------
  
  
  lc <- input
  # mean logcounts
  mat_brisc <- cbind(
    mat_brisc, 
    mean = rowMeans(lc)
  )
  # variance of logcounts
  mat_brisc <- cbind(
    mat_brisc, 
    var = rowVars(as.matrix(lc))
  )
  # spatial coefficient of variation
  mat_brisc <- cbind(
    mat_brisc, 
    spcov = sqrt(mat_brisc[, "sigma.sq"]) / mat_brisc[, "mean"]
  )
  
  
  # proportion of spatial variance out of total variance
  mat_brisc <- cbind(
    mat_brisc, 
    prop_sv = mat_brisc[, "sigma.sq"] / (mat_brisc[, "sigma.sq"] + mat_brisc[, "tau.sq"])
  )
  
  # ------------------------------------------
  # likelihood ratio (LR) statistics and tests
  # ------------------------------------------
  
  if (is(input, "SpatialExperiment")) {
    nrows <- nrow(spe)
    ncols <- ncol(spe)
  } else {
    nrows <- nrow(input)
    ncols <- ncol(input)
  }
  
  # calculate log likelihoods for nonspatial models
  
  loglik_lm <- vapply(seq_len(nrows), function(i) {
    y_i <- y[i, ]
    if (is.null(X)) {
      X <- rep(1, ncols)
    }
    # model formula without intercept to enable weighted model
    as.numeric(logLik(lm(y_i ~ X - 1)))
  }, numeric(1))
  
  mat_brisc <- cbind(
    mat_brisc, 
    loglik_lm = loglik_lm
  )
  
  # calculate LR statistics and tests (Wilks' theorem, asymptotic chi-square
  # with 2 degrees of freedom)
  
  LR_stat <- -2 * (mat_brisc[, "loglik_lm"] - mat_brisc[, "loglik"])
  
  pval <- 1 - pchisq(LR_stat, df = 2)
  padj <- p.adjust(pval, method = "BH")
  
  # rank SVGs according to LR statistics
  LR_rank <- rank(-1 * LR_stat)
  
  mat_brisc <- cbind(
    mat_brisc, 
    LR_stat = LR_stat, 
    rank = LR_rank, 
    pval = pval, 
    padj = padj
  )
  
  # --------------
  # return outputs
  # --------------
  
  if (is(input, "SpatialExperiment")) {
    # return in rowData of spe object
    stopifnot(nrow(spe) == nrow(mat_brisc))
    rowData(spe) <- cbind(rowData(spe), mat_brisc)
    spe
  } else {
    # return as numeric matrix
    stopifnot(nrow(input) == nrow(mat_brisc))
    rownames(mat_brisc) <- row_names
    mat_brisc
  }
}


