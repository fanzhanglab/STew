#' FindClusters()
#'
#' Perform cell clustering using different resolutions
#'
#' @param obj stCCA S3 Object
#' @param resolution_list List of resolutions to use for clustering
#'
#' @export
#'
#' @import future.apply
#' future
#' Seurat
FindClusters <- function(resolution_list, obj) {

  future::plan(multisession, workers = 8)
  snn <- obj$snn
  ids_cos <- Reduce(cbind, future_lapply(resolution_list, function(res_use) {
    Seurat:::RunModularityClustering(SNN = snn, modularity = 1,
                                     resolution = res_use, algorithm = 3, n.start = 10,
                                     n.iter = 10, random.seed = 0, print.output = FALSE,
                                 edge.file.name = NULL)
  }, future.seed = TRUE))
  colnames(ids_cos) <- sprintf("res_%.2f", resolution_list)
  ids_cos <- as.data.frame(ids_cos)

  obj$clusters <- ids_cos
  return(obj)
}
