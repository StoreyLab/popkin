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
get_A <- function(X, n = NA, loci_on_cols = FALSE, mem_factor = 0.7, mem_lim = NA) {
    # determine some behaviors depending on data type
    # first validate class and set key booleans
    isFn <- FALSE
    if (class(X) == 'function') {
        isFn <- TRUE
        if (is.na(n)) stop('missing number of individuals "n", which is required when X is a function.')
    } else if (class(X) == 'BEDMatrix') { # same as general matrix but transposed
        loci_on_cols <- TRUE # this is always imposed for this particular format!
    } else if (class(X) != 'matrix') {
        stop('X has unsupported class: ', class(X))
    } 
    
    # extract dimensions from data (not possible for function version)
    # also get individual names (IDs)
    namesX <- NULL # default
    if (isFn) {
        m <- NA # have to define as NA to pass to get_mem_lim_m below
    } else {
        if (loci_on_cols) {
            if (!is.na(n) && n != nrow(X)) 
                warning('User set number of samples that does not match X dimensions (will go with latter): ', n, ' != ', nrow(X))
            n <- nrow(X)
            m <- ncol(X)
            namesX <- rownames(X)
        } else {
            if (!is.na(n) && n != ncol(X)) 
                warning('User set number of samples that does not match X dimensions (will go with latter): ', n, ' != ', ncol(X))
            n <- ncol(X)
            m <- nrow(X)
            namesX <- colnames(X)
        }
    } 
    
    # initialize desired matrix
    A <- matrix(0, nrow = n, ncol = n)
    M <- matrix(0, nrow = n, ncol = n) # normalization now varies per individual pair (this tracks NAs, so subtract from overall m below)
    
    # transfer names from X to A if present
    # this will carry over all the way to the final kinship matrix!
    # (M need not have names at all)
    if (!is.null(namesX)) {
        colnames(A) <- namesX
        rownames(A) <- namesX
    }
    
    # infer the number of SNPs to break data into, since we're limited by memory
    # given fixed n, solve for m:
    # get maximum m (number of SNPs) given n and the memory requested
    data <- solve_m_mem_lim(
        n = n,
        m = m,
        mat_m_n = 1, # X (0.5) + ?
        mat_n_n = 1, # A + M (0.5 + 0.5)
        mem = mem_lim,
        mem_factor = mem_factor
    )
    mc <- data$mem_chunk

    # navigate chunks
    mci <- 1 # start of first chunk (needed for matrix inputs only; as opposed to function inputs)
    while (TRUE) { # start an infinite loop, break inside as needed
        if (isFn) {
            Xi <- X( mc ) # get next "mc" SNPs
            if (is.null(Xi)) break # stop when SNPs run out (only happens for functions X, not matrices)
        } else {
            # here m is known...
            if (mci > m) break # this means all SNPs have been covered!
            
            is <- mci:min(mci+mc-1, m) # range of SNPs to extract in this chunk
            
            if (loci_on_cols) {
                Xi <- t(X[, is, drop = FALSE]) # transpose for our usual setup
            } else  {
                Xi <- X[is, , drop = FALSE]
            }
            mci <- mci + mc # update starting point for next chunk! (overshoots at the end, that's ok)
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
    A/M - 1
}
