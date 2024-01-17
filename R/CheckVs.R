CheckVs <- function(v,x,z,K){ # If v is NULL, then get v as appropriate.
  ##print(list(v=v, x = x, z = z, K = K))
  library(RSpectra)
  if(!is.null(v) && !is.matrix(v)) v <- matrix(v,nrow=ncol(z))
  if(!is.null(v) && ncol(v)<K) v <- NULL
  if(!is.null(v) && ncol(v)>K) v <- matrix(v[,1:K],ncol=K)
  if(is.null(v) && ncol(z)>nrow(z) && ncol(x)>nrow(x)){
    v <- try(matrix(fastsvd(x,z)$v[,1:K],ncol=K), silent=TRUE)
    attempt <- 1
    while(("try-error" %in% class(v))  && attempt < 10){a
      v <- try(matrix(fastsvd(x,z)$v[,1:K],ncol=K), silent=TRUE)
      attempt <- attempt+1
    }
    if(attempt==10) stop("Problem computing SVD.")
  } else if (is.null(v) && (ncol(z)<=nrow(z) || ncol(x)<=nrow(x))){
    attempt <- 1
    v <- try(matrix(svds(t(x)%*%z, k = K)$v[,1:K],ncol=K), silent=TRUE)
    while(("try-error" %in% class(v)) && attempt<10){
      v <- try(matrix(svds(t(x)%*%z, k = K)$v[,1:K],ncol=K), silent=TRUE)
      attempt <- attempt+1
    }
    if(attempt==10) stop("Problem computing SVD.")
  }
  return(v)
}
