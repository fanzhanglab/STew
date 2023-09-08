#' STew S3 object
#'
#' @param count_exp count_exp data
#' @param spatial spatial data
#'
#' @return S3 Object
#' @export
STew_Obj = function(count_exp = NULL, spatial = NULL){
  stew = list()
  stew$count_exp <- count_exp
  stew$spatial <- spatial
  class(stew) <- 'STew_Object'
  return (stew)
}


