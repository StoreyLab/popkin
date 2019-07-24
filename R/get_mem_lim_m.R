# an internal constant shared by this function and get_mem_lim
GB <- 1024*1024*1024

# mem is in GB
get_mem_lim_m <- function(n, m = NA, mem = NA, factor = 0.7, verbose = FALSE) {
    # NOTE: in estimating the chunk size, we don't know ahead of time if there are NAs or not!
    # so these calculations assume M in getMAInt is an n*n matrix
    # if there aren't any NAs, M is a scalar and we end up underestimating memory usage (better than the other way around)

    # n must be defined now, or this doesn't work!
    if (missing(n))
        stop('`n` is required!')
    
    # try to get total memory from the system if mem wasn't specified, so it works reasonably
    if (is.na(mem)) {
        mem <- get_mem_lim(factor = factor, verbose = verbose) # infer a reasonable default from system!
    } else {
        # assuming mem is in GB, let's convert to bytes
        mem <- mem * GB
    }

    # OLD...
    # estimating total memory usage in bytes, for getMAInt
    # mem(X) = m*n*4+40 # genotypes are stored as ints!
    # mem(FHat) = m*8+40 # double vector
    # fixed cases double X at least temporarily?
    # X <- X-1 may also double X (but later, not at the same time)
    # mem(M) = n*n*4+40 # as int (RcppEigen), double otherwise
    # mem(is.na(X)) same as X
    # mem(!is.na(X)) same as X
    # mem(SA) = n*n*4+40 # as int (RcppEigen), double otherwise
    
    # mem = # version for int M,SA (RcppEigen), no more FHat
    # = 2*(m*n*4+40) + 2*(n*n*4+40)
    # = 2*m*n*4 + 2*n*n*4 + 4*40
    # = 8*m*n + 8*n*n + 160
    # = 8*(m*n + n*n + 20)
    # given fixed n, solve for m:
    # get maximum m (number of SNPs) given n and the memory requested
    data <- solve_m_mem_lim(
        mem = mem,
        n = n,
        m = m,
        mat_m_n = 1, # X (0.5) + ?
        mat_n_n = 1 # A + M (0.5 + 0.5)
    )
    #    m_chunk <- (mem/8 - 20 - n*n)/n

    if (verbose)
        message('Choice of chunk size should limit mem to about ', round( data$mem_chunk / GB, 2 ), ' GB')

    return(data$m_chunk) # return desired value
}

