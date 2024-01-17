#' spatial_adjacency_matrix()
#'
#' Creating spatial-driven adjacency matrix using cell spatial coordinates
#'
#' @param obj stCCA S3 Object
#'
#' @return A matrix that describes the adjacency relationships between cells. The entries in the matrix are binary values (0 or 1) that indicate whether cells are adjacent to each other or not.
#' @export
#'
#' @import MERINGUE
spatial_adjacency_matrix <- function(obj) {

  coordinate <- as.data.frame(obj$spatial)
  adjacency_matrix <- getSpatialNeighbors(coordinate)
  rownames(adjacency_matrix) <- rownames(coordinate)
  colnames(adjacency_matrix) <- rownames(coordinate)
  obj$adj_matrix <- adjacency_matrix
  return(obj)

}

