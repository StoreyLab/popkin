#' @useDynLib popkin
#' @importFrom Rcpp sourceCpp
NULL

# Compute A matrix from genotypes
#
# Given the biallelic genotypes of \eqn{n} individuals, this function returns the \eqn{n}-by-\eqn{n} matrix \eqn{A} that satisfies
# \deqn{E[A] = \alpha(\Phi - 1),}
# where \eqn{\Phi} is the kinship matrix and \eqn{\alpha} is a nuisance scaling factor (determined by the unknown ancestral allele frequencies of each locus).
# Thus a \eqn{\Phi} estimate can be recovered from \eqn{A} after a separate step that estimates \eqn{\alpha = -\min E[A]}{M = -min E[A]} (see \code{\link{min_mean_subpops}} for one example).
#
# The matrix X (or the vectors returned by the function X) must have values only in c(0,1,2,NA), encoded to count the number of reference alleles at the locus, or NA for missing data.
#
# @param X Genotype matrix, BEDMatrix object, or a function X(mc) that returns the genotype matrix of all individuals at mc successive loci, and NULL when no loci are left.
# @param n Number of individuals (required only when X is a function, ignored otherwise)
# @param loci_on_cols If true, X has loci on columns and individuals on rows; if false, loci are on rows and individuals on columns. Has no effect if X is a function.  If X is a BEDMatrix object, loci_on_cols=TRUE is set automatically.
# @param mem_lim Memory limit in GB used to calculate the "chunk size" (numbers of SNPs). Note memory usage is somewhat underestimated and is not controlled strictly.  Default is 1GB, except in linux it is the free memory in the system times 0.7.
#
# @return The A matrix.
#
# @examples
# # Construct toy data
# X <- matrix(c(0,1,2,1,0,1,1,0,2), nrow = 3, byrow = TRUE) # genotype matrix
#
# # NOTE: for BED-formatted input, use BEDMatrix!
# # "file" is path to BED file (excluding .bed extension)
# # library(BEDMatrix)
# # X <- BEDMatrix(file) # load genotype matrix object
#
# A <- get_A(X) # calculate A from genotypes
# 
get_A <- function(X, n_ind = NA, loci_on_cols = FALSE, mem_factor = 0.7, mem_lim = NA) {
    # for some more recent memory tests (internal hack)
    mem_debugging <- FALSE
    
    # determine some behaviors depending on data type
    # first validate class and set key booleans
    isFn <- FALSE
    if (class(X) == 'function') {
        isFn <- TRUE
        if (is.na(n_ind))
            stop('missing number of individuals "n", which is required when X is a function.')
    } else if (class(X) == 'BEDMatrix') { # same as general matrix but transposed
        loci_on_cols <- TRUE # this is always imposed for this particular format!
    } else if (!is.matrix(X)) {
        stop('X has unsupported class: ', class(X))
    } 
    
    # extract dimensions from data (not possible for function version)
    # also get individual names (IDs)
    names_X <- NULL # default
    if (isFn) {
        m_loci <- NA # have to define as NA to pass to get_mem_lim_m below
    } else {
        if (loci_on_cols) {
            if (!is.na(n_ind) && n_ind != nrow(X)) 
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
    
    # initialize desired matrix
    A <- matrix(0, nrow = n_ind, ncol = n_ind)
    M <- matrix(0, nrow = n_ind, ncol = n_ind) # normalization now varies per individual pair (this tracks NAs, so subtract from overall m below)
    
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
        mat_m_n = 1.5, # X (0.5) + ? + (BEDMatrix seems to consume too much additional memory, so be extra conservative overall)
        mat_n_n = 1, # A + M (0.5 + 0.5)
        mem = mem_lim,
        mem_factor = mem_factor
    )
    m_chunk <- data$m_chunk
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

        # before passing along to my RcppEigen code, I need to make sure the genotypes are treated by R as integers or RcppEigen dies on me
        # I'm not sure if this always works though...
        if (storage.mode(Xi) != 'integer')
            storage.mode(Xi) <- 'integer'
        
        # solve chunk using very efficient RcppEigen code!
        # it is both runtime and memory efficient!
        obj <- getMAInt(Xi)
        # increment each part
        A <- A + obj$SA # increment sum of A_{ijk}'s
        M <- M + obj$M # increment M's too
    }

    # return final estimate!
    A / M - 1
}
