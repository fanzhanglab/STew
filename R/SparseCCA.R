SparseCCA <- function(x,y,v,penaltyx, penaltyz,niter,trace = FALSE, upos, uneg, vpos, vneg){
    vold <- rnorm(length(v))
    u <- rnorm(ncol(x))
    for(i in 1:niter){
        if(sum(is.na(u))>0 || sum(is.na(v))>0){
            v <- rep(0, length(v))
            vold <- v
        }
        if(sum(abs(vold-v))>1e-6){
            if(trace) cat(i,fill=F)
            # Update u #
            unew <- rep(NA, ncol(x))

            #argu <- t(x)%*%(y%*%v)
            argu <- matrix(y%*%v,nrow=1)%*%x
            if(upos) argu <- pmax(argu,0)
            if(uneg) argu <- pmin(argu,0)
            lamu <- BinarySearch(argu,penaltyx*sqrt(ncol(x)))
            su <- soft(argu,lamu)
            u <-  matrix(su/l2n(su), ncol=1)

            # Done updating u #
            # Update v #
            vnew <- rep(NA, ncol(y))

            vold <- v
            #argv <- (t(u)%*%t(x))%*%y
            argv <- matrix(x%*%u,nrow=1)%*%y
            if(vpos) argv <- pmax(argv,0)
            if(vneg) argv <- pmin(argv,0)
            lamv <- BinarySearch(argv,penaltyz*sqrt(ncol(y)))
            sv <- soft(argv, lamv)
            v <-  matrix(sv/l2n(sv),ncol=1)

            # Done updating v #
        }
    }
    if(trace) cat(fill=T)
    # Update d #
    d <-  sum((x%*%u)*(y%*%v))
    # Done updating d #
    if(sum(is.na(u))>0 || sum(is.na(v))>0){
        u <- matrix(rep(0,ncol(x)),ncol=1)
        v <- matrix(rep(0,ncol(y)),ncol=1)
        d <- 0
    }
    return(list(u=u,v=v,d=d))
}
