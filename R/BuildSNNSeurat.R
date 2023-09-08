#' BuildSNNSeurat()
#'
#' Build SNN matrix by inputting canonical variates of spatial-driven canonical loading
#'
#' @param k.param  Maximum number of nearest neighbors to compute. Default - 30.
#' @param prune.SNN Sets the cutoff for acceptable Jaccard index when computing the neighborhood overlap for the SNN construction. Any edges with values less than or equal to this will be set to 0 and removed from the SNN graph. Essentially sets the stringency of pruning (0 --- no pruning, 1 --- prune everything). Default - 1/15.
#' @param nn.eps Error bound: default of 0.0 implies exact nearest neighbor search
#' @param obj stCCA S3 Object
#' @param cols columns selected
#'
#' @export
#'
#'
#' @import Seurat
#' future
#' future.apply
#' RANN
BuildSNNSeurat <- function(obj, k.param = 30, prune.SNN = 1/15, nn.eps = 0, cols = NULL) {

  environment(BuildSNNSeurat) <- asNamespace("Seurat")

  future::plan(multisession, workers = 4)

  data.use <- obj$joint_emb[, cols]
  my.knn <- nn2(data = data.use, k = k.param, searchtype = "standard", eps = nn.eps)
  nn.ranked <- my.knn$nn.idx

  snn_res <- Seurat:::ComputeSNN(nn_ranked = nn.ranked, prune = prune.SNN)
  rownames(snn_res) <- row.names(data.use)
  colnames(snn_res) <- row.names(data.use)

  obj$snn <- snn_res

  return(obj)
}
