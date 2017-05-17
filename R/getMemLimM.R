getMemLimM <- function(m=NA, n=NA, mem=NA, factor=0.7, verbose=FALSE) {
    ## NOTE: in estimating the chunk size, we don't know ahead of time if there are NAs or not!
    ## so these calculations assume M in getMAInt is an n*n matrix
    ## if there aren't any NAs, M is a scalar and we end up underestimating memory usage (better than the other way around)

    ## n must be defined now, or this doesn't work!
    if (is.na(n)) stop('Fatal: n must be defined!')
    
    ## try to get total memory from the system if mem wasn't specified, so it works reasonably
    if (is.na(mem)) {
        mem <- getMemLim(factor=factor, verbose=verbose) # infer a reasonable default from system!
    } else {
        ## assuming mem is in GB, let's convert to bytes
        mem <- mem*1024*1024*1024
    }

    ## OLD...
    ## estimating total memory usage in bytes, for getMAInt
    ## mem(X) = m*n*4+40 # genotypes are stored as ints!
    ## mem(FHat) = m*8+40 # double vector
    ## fixed cases double X at least temporarily?
    ## X <- X-1 may also double X (but later, not at the same time)
    ## mem(M) = n*n*4+40 # as int (RcppEigen), double otherwise
    ## mem(is.na(X)) same as X
    ## mem(!is.na(X)) same as X
    ## mem(SA) = n*n*4+40 # as int (RcppEigen), double otherwise
    
    ## mem = # version for int M,SA (RcppEigen), no more FHat
    ## = 2*(m*n*4+40) + 2*(n*n*4+40)
    ## = 2*m*n*4 + 2*n*n*4 + 4*40
    ## = 8*m*n + 8*n*n + 160
    ## = 8*(m*n + n*n + 20)
    ## given fixed n, solve for m:
    ## get maximum m (number of SNPs) given n and the memory requested
    mc <- (mem/8 - 20 - n*n)/n

    ## NOTE m may be missing if X is a function, so we can't make these simplifying decisions (to balance load) without m in that case...
    if (!is.na(m)) {
        if (m < mc) {
            mc <- m # use the smaller one
        } else {
            ## should "redistribute" based on number of chunks, to lower memory even more per iteration
            mc <- ceiling( m/ceiling(m/mc) ) # this lowers mc even more, balances load better
        }
    }
    
    memAct <- 8*(mc*n + n*n + 20)
    if (verbose) message('Choice of chunk size should limit mem to about ', round( memAct/(1024*1024*1024), 2 ), ' GB')

    return(mc) # return desired value
}

