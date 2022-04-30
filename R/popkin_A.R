#' Compute popkin's `A` and `M` matrices from genotypes
#'
#' This function returns lower-level, intermediate calculations for the main `popkin` function.
#' These are not intended for most users, but rather for researchers studying the estimator.
#' 
#' @inheritParams popkin
#'
#' @return A named list containing:
#'
#' - `A`: n-by-n matrix, for individuals `j` and `k`, of average `w_i * ( ( x_ij - 1 ) * ( x_ik - 1 ) - 1)` values across all loci `i` in `X`; if `mean_of_ratios = FALSE`, `w_i = 1`, otherwise `w_i = 1 / (p_est_i * (1 - p_est_i) )` where `p_est_i` is the reference allele frequency.
#' - `M`: n-by-n matrix of sample sizes (number of loci with non-missing individual `j` and `k` pairs, used to normalize `A`)
#'
#' @examples
#' # Construct toy data
#' X <- matrix(c(0,1,2,1,0,1,1,0,2), nrow = 3, byrow = TRUE) # genotype matrix
#'
#' # NOTE: for BED-formatted input, use BEDMatrix!
#' # "file" is path to BED file (excluding .bed extension)
#' # library(BEDMatrix)
#' # X <- BEDMatrix(file) # load genotype matrix object
#'
#' obj <- popkin_A(X) # calculate A and M from genotypes
#' A <- obj$A
#' M <- obj$M
#'
#' @seealso
#' The main [popkin()] function (a wrapper of this `popkin_A` function and [popkin_A_min_subpops()] to estimate the minimum `A` value).
#'
#' @export
popkin_A <- function(
                     X,
                     n = NA,
                     loci_on_cols = FALSE,
                     mean_of_ratios = FALSE,
                     mem_factor = 0.7,
                     mem_lim = NA,
                     m_chunk_max = 1000 # gave good performance in tests
                     ) {
    if ( missing( X ) )
        stop( 'Genotype matrix `X` is required!' )

    # internally code uses this more specific variable name
    n_ind <- n
    
    # for some more recent memory tests (internal hack)
    mem_debugging <- FALSE
    
    # determine some behaviors depending on data type
    # first validate class and set key booleans
    # NOTE BEDMatrix is also class 'matrix', so have to test in this order
    isFn <- FALSE
    if (is.function(X)) {
        isFn <- TRUE
        if ( is.na( n_ind ) )
            stop('missing number of individuals "n", which is required when X is a function.')
    } else if ('BEDMatrix' %in% class(X)) { # same as general matrix but transposed
        loci_on_cols <- TRUE # this is always imposed for this particular format!
    } else if (!is.matrix(X)) {
        stop('X has unsupported class: ', toString( class(X) ) )
    } 
    
    # extract dimensions from data (not possible for function version)
    # also get individual names (IDs)
    names_X <- NULL # default
    if (isFn) {
        m_loci <- NA # have to define as NA to pass to get_mem_lim_m below
    } else {
        if (loci_on_cols) {
            if ( !is.na( n_ind ) && n != nrow( X ) ) 
                warning('User set number of samples that does not match X dimensions (will go with latter): ', n_ind, ' != ', nrow(X))
            n_ind <- nrow(X)
            m_loci <- ncol(X)
            names_X <- rownames(X)
        } else {
            if (!is.na(n_ind) && n_ind != ncol(X)) 
                warning('User set number of samples that does not match X dimensions (will go with latter): ', n_ind, ' != ', ncol(X))
            n_ind <- ncol(X)
            m_loci <- nrow(X)
            names_X <- colnames(X)
        }
    } 
    
    # initialize desired matrices
    A <- matrix(0, nrow = n_ind, ncol = n_ind)
    M <- matrix(0, nrow = n_ind, ncol = n_ind)
    
    # transfer names from X to A if present
    # this will carry over all the way to the final kinship matrix!
    # (M need not have names at all)
    if (!is.null(names_X)) {
        colnames(A) <- names_X
        rownames(A) <- names_X
    }
    
    # infer the number of SNPs to break data into, since we're limited by memory
    # given fixed n, solve for m:
    # get maximum m (number of SNPs) given n and the memory requested
    data <- solve_m_mem_lim(
        n = n_ind,
        m = m_loci,
        mat_m_n = 1.5, # X (0.5) + is.na(X) (0.5) + ? (BEDMatrix seems to consume too much additional memory, so be extra conservative overall)
        mat_n_n = 1, # A + M (0.5 + 0.5; tmp copies?)
        vec_m = 1.5, # Vi, indexes_not_fixed (0.5)
        mem = mem_lim,
        mem_factor = mem_factor
    )
    m_chunk <- data$m_chunk
    # cap value to a nice performing value (very good speed, minimal memory)
    if ( m_chunk > m_chunk_max )
        m_chunk <- m_chunk_max
    if (mem_debugging) {
        # hack: report things to troubleshoot
        message('mem_lim: ', round(data$mem_lim / GB, 1), ' GB')
        message('mem_chunk: ', round(data$mem_chunk / GB, 1), ' GB')
        message('m_chunk: ', m_chunk)
    }
    
    # navigate chunks
    i_chunk <- 1 # start of first chunk (needed for matrix inputs only; as opposed to function inputs)
    while (TRUE) { # start an infinite loop, break inside as needed
        if (isFn) {
            Xi <- X( m_chunk ) # get next "m_chunk" SNPs
            if (is.null(Xi)) break # stop when SNPs run out (only happens for functions X, not matrices)
        } else {
            # hacky message
            if (mem_debugging)
                message('chunk started!')
            # here m is known...
            # this means all SNPs have been covered!
            if (i_chunk > m_loci)
                break
            
            # range of SNPs to extract in this chunk
            indexes_loci_chunk <- i_chunk : min(i_chunk + m_chunk - 1, m_loci)
            
            if (loci_on_cols) {
                Xi <- t(X[, indexes_loci_chunk, drop = FALSE]) # transpose for our usual setup
            } else  {
                Xi <- X[indexes_loci_chunk, , drop = FALSE]
            }
            i_chunk <- i_chunk + m_chunk # update starting point for next chunk! (overshoots at the end, that's ok)
        }

        # standard mean times half
        Vi <- rowMeans(Xi, na.rm = TRUE) / 2
        # variance estimate (length-i_chunk vector; factor of 2 or 4 cancels out in the end so it is ignored)
        Vi <- Vi * ( 1 - Vi )
        
        # in the mean_of_ratios formulation it's extra critical to handle fixed loci correctly... but let's just do the same thing both ways
        indexes_not_fixed <- Vi > 0 # these are the good cases
        # filter everything if needed (will increase memory but it's easy to see things are handled correctly)
        if (any(!indexes_not_fixed)) {
            Xi <- Xi[indexes_not_fixed,]
            Vi <- Vi[indexes_not_fixed]
        }

        # only appears in this form in the rest of the code, MOR only
        if (mean_of_ratios)
            Vi <- 1 / sqrt( Vi )
        
        # center before cross product...
        Xi <- Xi - 1
        
        # some complicated steps are only necessary when there's missing data...
        if (anyNA(Xi)) {
            M <- M + crossprod( !is.na(Xi) ) # add non-missingness counts to pair count matrix
            # need this awkward adjustment for MOR
            if (mean_of_ratios)
                A <- A - crossprod( (!is.na(Xi)) * Vi )
            
            # before applying cross product, to prevent NA errors, just set those values to zero and it works out!
            Xi[is.na(Xi)] <- 0
        } else {
            M <- M + nrow(Xi) # this is correct denominator
            # need this adjustment for MOR
            if (mean_of_ratios)
                A <- A - sum( Vi^2 )
        }
        
        if (mean_of_ratios)
            # will average into A but scaling first!
            Xi <- Xi * Vi # this works! (scales each row as needed)
        
        # cross product matrix at this SNP, add to running sum.  We'll add an extra -1 later... (this is computationally faster and maybe even more numerically stable)
        A <- A + crossprod(Xi)
        # NOTE: M and m count the same loci (including fixed loci in this case; as long as there's consistency the difference just cancels out as expected)
    }
    
    # calculate final (properly averaged) estimate!
    A <- A / M
    # add -1 only for ROM (has already been included for MOR)
    if ( !mean_of_ratios )
        A <- A - 1
    
    # return parts of interest
    return(
        list(
            A = A,
            M = M
        )
    )
}
