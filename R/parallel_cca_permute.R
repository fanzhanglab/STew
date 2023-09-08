#' Selecting tuning parameters for sparse CCA
#'
#' @references We adapted and modified the code from the CCA function in the PMA package to be compatible and scalable with single-cell multimodal integration through parallel computing.
#'
#' @param x Expression-driven adjacency matrix
#' @param z Spatial-driven adjacency matrix
#' @param penaltyxs -desc
#' @param penaltyzs -desc
#' @param niter The number of iterations to be performed. Default - 20.
#' @param v -desc
#' @param trace -desc
#' @param nperms -desc
#' @param standardize -desc
#' @param upos -desc
#' @param uneg -desc
#' @param vpos -desc
#' @param vneg -desc
#' @param outcome -desc
#' @param y -desc
#' @param cens -desc
#' @param obj stCCA S3 Object
#'
#' @export
parallel_cca_permute <- function(x, z, obj, penaltyxs=NULL, penaltyzs=NULL, niter=3,v=NULL,trace=FALSE,nperms=25, standardize=TRUE, upos=FALSE, uneg=FALSE, vpos=FALSE, vneg=FALSE, outcome=NULL, y=NULL, cens=NULL){

    u <- NULL
    call <- match.call()
    if(!is.null(penaltyxs) && !is.null(penaltyzs) && length(penaltyxs)>1 && length(penaltyzs)>1 && length(penaltyxs)!=length(penaltyzs)) stop("Penaltyxs and Penaltyzs must be same length, or one must have length 1. This is because tuning parameters are considered in pairs.")
    if(is.null(penaltyxs)) penaltyxs <- seq(.2,.8,len=5)
    if(is.null(penaltyzs)) penaltyzs <- seq(.2,.8,len=5)
    if(length(unique(penaltyxs))==1 && length(unique(penaltyzs))==1){
        out <- CCA.permute.justone(x=x,z=z,penaltyx=penaltyxs[1],penaltyz=penaltyzs[1],niter=niter,v=v,trace=trace,nperms=nperms,standardize=standardize,upos=upos,uneg=uneg,vpos=vpos,vneg=vneg,outcome=outcome,y=y,cens=cens)
    }
    if(length(penaltyxs)==1 && length(penaltyzs)>1) out <- CCA.permute.zonly(x=x,z=z,penaltyx=penaltyxs,penaltyzs=penaltyzs,niter=niter,v=v,trace=trace,nperms=nperms,standardize=standardize,upos=upos,uneg=uneg,vpos=vpos,vneg=vneg,outcome=outcome,y=y,cens=cens)
    if(length(penaltyxs)>1 && length(penaltyzs)==1) out <- CCA.permute.xonly(x=x,z=z,penaltyxs=penaltyxs,penaltyz=penaltyzs,niter=niter,v=v,trace=trace,nperms=nperms,standardize=standardize,upos=upos,uneg=uneg,vpos=vpos,vneg=vneg,outcome=outcome,y=y,cens=cens)
    if(length(penaltyzs)>1 && length(penaltyxs)>1) out <- CCA.permute.both(x=x,z=z,penaltyxs=penaltyxs,penaltyzs=penaltyzs,niter=niter,v=v,trace=trace,nperms=nperms,standardize=standardize,upos=upos,uneg=uneg,vpos=vpos,vneg=vneg,outcome=outcome,y=y,cens=cens)
    out$call <- call
    out$upos <- upos
    out$uneg <- uneg
    out$vpos <- vpos
    out$vneg <- vneg
    class(out) <- "parallel_cca_permute"
    obj$bestpenaltyx <- out$bestpenaltyx
    obj$bestpenaltyz <- out$bestpenaltyz
    obj$v.init <- out$v.init
    return(obj)
}


