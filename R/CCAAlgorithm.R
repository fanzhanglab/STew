CCAAlgorithm <- function(x,z,v,penaltyx,penaltyz,K,niter,trace,upos,uneg,vpos,vneg){
    
    if(K>1) v.init <- v[apply(z^2,2,sum)!=0,]
    if(K==1) v.init <- v[apply(z^2,2,sum)!=0]
    
    v.init <- matrix(v.init,ncol=K)
    u=v=d=NULL
    xres <- x; zres <- z
    xres <- x[,apply(x^2,2,sum)!=0]
    zres <- z[,apply(z^2,2,sum)!=0]
    for(k in 1:K){
        if(vpos && sum(abs(v.init[v.init[,k]>0,k]))<sum(abs(v.init[v.init[,k]<0,k]))) v.init[,k] <- -v.init[,k]
        if(vneg && sum(abs(v.init[v.init[,k]<0,k]))<sum(abs(v.init[v.init[,k]>0,k]))) v.init[,k] <- -v.init[,k]
        out <- SparseCCA(xres,zres,v.init[,k],penaltyx, penaltyz,niter,trace, upos, uneg, vpos, vneg)
        coef <- out$d
        d <- c(d, coef)
        xres <- rbind(xres, sqrt(coef)*t(out$u))
        zres <- rbind(zres, -sqrt(coef)*t(out$v))
        u <- cbind(u, out$u)
        v <- cbind(v, out$v)
    }
    
    ubig <- u
    vbig <- v
    
    ubig <- matrix(0,nrow=ncol(x),ncol=K)
    ubig[apply(x^2,2,sum)!=0,] <- u
    
    
    vbig <- matrix(0,nrow=ncol(z),ncol=K)
    vbig[apply(z^2,2,sum)!=0,] <- v
    
    return(list(u=ubig,v=vbig,d=d))
}
