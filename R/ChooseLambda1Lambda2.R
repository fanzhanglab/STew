library(stats)
ChooseLambda1Lambda2 <- function(y,chrom=NULL){ # using as of 8/25/09
  if(is.null(chrom)) chrom <- rep(1, length(y))
  y.lo <- lowess(y, f=min(50/length(y),.2))$y
  y.lo <- soft(y.lo, quantile(abs(y.lo), .45))
  s1=sum(abs(y.lo))
  lambdas=exp(seq(log(10^(-4)), .2, len=50))
  out <- matrix(NA, nrow=length(y),ncol=length(lambdas))
  for(chr in unique(chrom)){
    out[chrom==chr,] <- t(FLSA(y[chrom==chr],lambda1=0,lambda2=lambdas))
  }
  for(i in 1:length(lambdas)) out[,i] <- soft(out[,i], lambdas[i])
  s1s <- apply(abs(out),2,sum)
  errs <- abs(s1s-s1)
  bestlam <- lambdas[which.min(errs)]
  return(bestlam)
}
