#' joint_embedding()
#'
#' Extracting spatial-informed canonical loading
#'
#' @param obj stCCA S3 Object
#'
#' @export
joint_embedding <- function(obj) {

  stcc <- obj$optimal_cca
  meta_cell <- stcc$v
  rownames(meta_cell) <- stcc$znames
  meta_cell <- data.frame(meta_cell)

  for ( col in 1:ncol(meta_cell)){
    colnames(meta_cell)[col] <-  sub("X", "", colnames(meta_cell)[col])
  }  # rename columns

  colnames(meta_cell)[1:ncol(meta_cell)] <- paste0("stCC",      colnames(meta_cell)[1:ncol(meta_cell)], sep="") # rename columns
  obj$joint_emb <- meta_cell
  return(obj)

}

