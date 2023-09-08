##' cca_main()
#'
#' @description Perform Sparse CCA using the optimal tuning parameters
#'
#'
#' @param x Expression-driven adjacency matrix
#' @param z Spatial-driven adjacency matrix
#' @param penaltyx The penalty to be applied to the matrix x, i.e. the penalty that results in the canonical vector u. If typex is "standard" then the L1 bound on u is penaltyx*sqrt(ncol(x)). In this case penaltyx must be between 0 and 1 (larger L1 bound corresponds to less penalization).
#' @param penaltyz The penalty to be applied to the matrix z, i.e. the penalty that results in the canonical vector v. If typez is "standard" then the L1 bound on v is penaltyz*sqrt(ncol(z)). In this case penaltyz must be between 0 and 1 (larger L1 bound corresponds to less penalization).
#' @param K The number of u's and v's desired; that is, the number of canonical vectors to be obtained.
#' @param niter The number of iterations to be performed. Default - 20.
#' @param v -desc
#' @param trace -desc
#' @param standardize -desc
#' @param xnames -desc
#' @param znames -desc
#' @param upos -desc
#' @param uneg -desc
#' @param vpos -desc
#' @param vneg -desc
#' @param outcome -desc
#' @param y -desc
#' @param cens -desc
##' @param obj stCCA S3 Object
#'
#' @export
cca_main <- function(x, z, obj = NULL, penaltyx=NULL, penaltyz=NULL, K=1, niter=15, v=NULL, trace=TRUE, standardize=TRUE, xnames=colnames(x), znames=colnames(z), upos=FALSE, uneg=FALSE, vpos=FALSE, vneg=FALSE, outcome=NULL, y=NULL, cens=NULL){


    call <- match.call()
    if(sum(is.na(x))+sum(is.na(z)) > 0) stop("Cannot have NAs in x or z")
    if(nrow(x)!=nrow(z)) stop("x and z must have same number of rows")
    if(standardize){
        sdx <- apply(x,2,sd)
        sdz <- apply(z,2,sd)
        if(min(sdx)==0) stop("Cannot standardize because some of the columns of x have std. dev. 0")
        if(min(sdz)==0) stop("Cannot standardize because some of the columns of z have std. dev. 0")
        x <- scale(x,TRUE,sdx)
        z <- scale(z,TRUE,sdz)
    }
    if(!is.null(outcome)){
        pheno.out <- CCAPhenotypeZeroSome(x,z,y,qt=.8, cens=cens, outcome=outcome)
        x <- pheno.out$x
        z <- pheno.out$z
    }
    v <- CheckVs(v,x,z,K)
    if(is.null(penaltyx)){penaltyx <- .3#pmax(1.001,.3*sqrt(ncol(x)))/sqrt(ncol(x))
    }
    if(is.null(penaltyz)){penaltyz <- .3#pmax(1.001,.3*sqrt(ncol(z)))/sqrt(ncol(z))
    }
    if(!is.null(penaltyx)){
        if(penaltyx<0 || penaltyx>1) stop("Penaltyx must be between 0 and 1.")

    }
    if(!is.null(penaltyz)){
        if(penaltyz<0 || penaltyz>1) stop("Penaltyz must be between 0 and 1.")

    }
    out <- CCAAlgorithm(x=x,z=z,v=v,penaltyx=penaltyx,penaltyz=penaltyz,K=K,niter=niter,trace=trace,upos=upos,uneg=uneg,vpos=vpos,vneg=vneg)
    out$outcome <- outcome
    out$call <- call
    out$xnames <- xnames
    out$znames <- znames
    out$penaltyx<-penaltyx
    out$penaltyz<-penaltyz
    out$K <- K
    out$niter <- niter
    out$upos <- upos
    out$uneg <- uneg
    out$vpos <- vpos
    out$vneg <- vneg
    out$xnames <- xnames
    out$znames <- znames
    out$v.init <- v
    out$cors <- numeric(K)
    for(k in 1:K){
        if(sum(out$u[,k]!=0)>0 && sum(out$v[,k]!=0)>0) out$cors[k] <- cor(x%*%out$u[,k],z%*%out$v[,k])
    }
    class(out) <- "cca_main"
    if(!is.null(obj)) {
    obj$optimal_cca <- out
    return(obj)
    } else {
      return(out)
    }
}

