library(future.apply)
options(future.seed = TRUE, future.globals.maxSize = 1024 * 1024 * 1024)

CCA.permute.both=
    function(x,z,penaltyxs,penaltyzs,niter,v,trace,nperms,standardize,upos,uneg,vpos,vneg,outcome,y,cens, seed = 1234){
        set.seed(seed)
        call <- match.call()
        if(standardize){
            x <- scale(x,TRUE,TRUE)
            z <- scale(z,TRUE,TRUE)
        }
        v <- CheckVs(v,x,z,1)
        ccperms=nnonzerous.perms=nnonzerovs.perms=matrix(NA, length(penaltyxs), nperms)
        ccs=nnonzerous=nnonzerovs=numeric(length(penaltyxs))


        # Generate all combinations of indices
        loop_indices <- expand.grid(i = 1:nperms, j = 1:length(penaltyxs))

        # Process each combination in parallel
        results_list <- future_lapply(1:nrow(loop_indices), function(idx) {
            i <- loop_indices[idx, "i"]
            j <- loop_indices[idx, "j"]

            # if (trace && .Platform$OS.type != "windows") cat("\n Permutation ", i, " out of ", nperms, " - Penalty index: ", j, "\n")
            sampz <- sample(1:nrow(z))
            sampx <- sample(1:nrow(x))

            if (i == 1) {
                out <- cca_main(x, z, penaltyx = penaltyxs[j], penaltyz = penaltyzs[j], y = y, outcome = outcome, cens = cens, niter = niter, v = v, trace = FALSE, upos = upos, uneg = uneg, vpos = vpos, vneg = vneg, standardize = FALSE)
                nnonzerous[j] <- sum(out$u != 0)
                nnonzerovs[j] <- sum(out$v != 0)
                if (mean(out$u == 0) != 1 && mean(out$v == 0) != 1) {
                    ccs[j] <- cor(x %*% out$u, z %*% out$v)
                } else {
                    ccs[j] <- 0
                }
            }
            out <- cca_main(x[sampx, ], z[sampz, ], penaltyx = penaltyxs[j], penaltyz = penaltyzs[j], y = y, outcome = outcome, cens = cens, niter = niter, v = v, trace = FALSE, upos = upos, uneg = uneg, vpos = vpos, vneg = vneg, standardize = FALSE)
            nnonzerous_perm <- sum(out$u != 0)
            nnonzerovs_perm <- sum(out$v != 0)
            if (mean(out$u == 0) != 1 && mean(out$v == 0) != 1) {
                ccperm <- cor(x[sampx, ] %*% out$u, z[sampz, ] %*% out$v)
            } else {
                ccperm <- 0
            }

            return(list(j = j, i = i, ccperm = ccperm, nnonzerous_perm = nnonzerous_perm, nnonzerovs_perm = nnonzerovs_perm, ccs = ifelse(i == 1, ccs[j], NA)))
        }, future.seed = TRUE)


        for (result in results_list) {
            j <- result$j
            i <- result$i
            ccperms[j, i] <- result$ccperm
            nnonzerous.perms[j, i] <- result$nnonzerous_perm
            nnonzerovs.perms[j, i] <- result$nnonzerovs_perm
            if (i == 1) {
                ccs[j] <- result$ccs
            }
        }
        # if(trace && .Platform$OS.type=="windows") close(pb)
        cc.norm <- ftrans(ccs)
        ccperm.norm <- ftrans(ccperms)
        zstats <- (cc.norm - rowMeans(ccperm.norm))/(apply(ccperm.norm,1,sd) + .05)
        # 0.05 added to the denominator to avoid getting zstat of INFINITY
        if(trace) cat(fill=T)
        pvals <- apply(sweep(ccperms,1,ccs,"-")>=0,1,mean)
        results <- list(zstats=zstats,penaltyxs=penaltyxs, penaltyzs=penaltyzs,bestpenaltyx=penaltyxs[which.max(zstats)], bestpenaltyz=penaltyzs[which.max(zstats)], cors=ccs, corperms=ccperms, ft.cors=cc.norm,ft.corperms=rowMeans(ccperm.norm),nnonzerous=nnonzerous,nnonzerovs=nnonzerovs, nnonzerous.perm=rowMeans(nnonzerous.perms),nnonzerovs.perm=rowMeans(nnonzerovs.perms),call=call,v.init=v,pvals=pvals,nperms=nperms,pvalbestz=pvals[which.max(zstats)])
        return(results)
    }
