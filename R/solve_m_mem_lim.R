# an internal constant shared by this function and get_mem_lim
GB <- 1024 * 1024 * 1024

# general function shared with other, related projects/packages
# (popkin, popkinsuppl, fit to kinship, rgls)
# mem is in GB, if passed
# matrix counts should be integers for double numerical data, or halves for integer data (this is not tested)
solve_m_mem_lim <- function(
                            n,
                            m = NA,
                            mat_m_n = 0,
                            mat_n_n = 0,
                            vec_m = 0,
                            vec_n = 0,
                            mem = NA,
                            mem_factor = 0.7
                            ) {
    # mandatory arguments
    if (missing(n))
        stop('`n` is required!')

    # check that data makes sense
    # everything has to be non-negative (a few things strictly positive)
    if (n <= 0)
        stop('`n` must be positive!  Passed ', n)
    if (!is.na(m) && m <= 0) # don't test if NA
        stop('`m` must be positive!  Passed ', m)
    if (mat_m_n < 0)
        stop('`mat_m_n` must be non-negative!  Passed ', mat_m_n)
    if (mat_n_n < 0)
        stop('`mat_n_n` must be non-negative!  Passed ', mat_n_n)
    if (vec_m < 0)
        stop('`vec_m` must be non-negative!  Passed ', vec_m)
    if (vec_n < 0)
        stop('`vec_n` must be non-negative!  Passed ', vec_n)

    if (is.na(mem)) {
        # get available memory if it's missing, times the factor as usual
        mem <- get_mem_lim(factor = mem_factor)
    } else {
        # sanity checks
        if (mem <= 0)
            stop('`mem` must be positive!  Passed ', mem)
        # assuming mem is in GB, let's convert to bytes
        mem <- mem * GB
    }
    
    # another sanity check
    # there is no problem if there are no `m`s, (below the denominator becomes zero)
    if (mat_m_n == 0 && vec_m == 0)
        stop('At least one of `mat_m_n` or `vec_m` must be non-zero!  (there is no `m` to solve for otherwise)')

    # memory overheads per type
    # library(pryr)
    # object_size( matrix( nrow = 0, ncol = 0 ) )
    # mo <- 216 # matrix overhead
    # object_size( numeric() )
    # ao <- 48  # array overhead
    mo8 <- 27 # mo/8 shortcut
    ao8 <- 6 # ao/8 shortcut
    
    # the forward formula is
    # (all double matrices/vectors)
    ## mem <- (
    ##     + mat_m_n * (m * n + mo8)
    ##     + mat_n_n * (n * n + mo8)
    ##     + vec_m * (m + ao8)
    ##     + vec_n * (n + ao8)
    ## ) * 8
    # ignoring overheads for extra clarity (easier to think about int cases)
    ## mem <- (
    ##     + mat_m_n * m * n
    ##     + mat_n_n * n * n
    ##     + vec_m * m
    ##     + vec_n * n
    ## ) * 8

    # turn n into double
    # to prevent integer overflows (experienced with `n * n` in particular, below)
    n <- n + 0.0
    
    # reversed we get (with overheads)
    m_chunk <- (
        + mem / 8
        - mat_m_n * mo8
        - vec_m * ao8
        - mat_n_n * (n * n + mo8)
        - vec_n * (n + ao8)
    ) / (
        + mat_m_n * n
        + vec_m
    )
    ## # ignoring overheads
    ## m_chunk <- (
    ##     + mem / 8
    ##     - mat_n_n * n * n
    ##     - vec_n * n
    ## ) / (
    ##     + mat_m_n * n
    ##     + vec_m
    ## )

    # if n is large and the memory too low, it's possible that m_chunk is negative
    # test for that and die if that is so
    if (m_chunk < 0)
        stop('The resulting `m_chunk` was negative!  This is because either `mat_n_n` or `vec_n` are non-zero and `n` alone is too large for the available memory (even for `m_chunk == 0`).  The solution is to free more memory (ideal) or to reduce `n` if possible.')
    
    # NOTE m may be missing if X is a function, so we can't make these simplifying decisions (to balance load) without m in that case...
    # but m is not missing for BEDMatrix and regular R matrices (most common cases by far)
    if (!is.na(m)) {
        # note both cases here return integers
        if (m < m_chunk) {
            m_chunk <- m # use the smaller one
        } else {
            # should "redistribute" based on number of chunks, to lower memory even more per iteration
            m_chunk <- ceiling( m / ceiling( m / m_chunk ) ) # this lowers m_chunk even more, balances load better
        }
    } else {
        # just round down here, to have integers in all cases
        m_chunk <- floor( m_chunk )
    }

    # actual memory in use per chunk, in bytes
    mem_chunk <- (
        + mat_m_n * (m_chunk * n + mo8)
        + mat_n_n * (n * n + mo8)
        + vec_m * (m_chunk + ao8)
        + vec_n * (n + ao8)
    ) * 8
    # ignoring overheads
    ## mem_chunk <- (
    ##     + mat_m_n * m_chunk * n
    ##     + mat_n_n * n * n
    ##     + vec_m * m_chunk
    ##     + vec_n * n
    ## ) * 8
    
    # return 
    return(
        list(
            m_chunk = m_chunk,
            mem_chunk = mem_chunk,
            mem_lim = mem # report the reading from get_mem_lim, or whatever was set by user but in bytes (rarer)
        )
    )
}
