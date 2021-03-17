#' Estimate kinship from a genotype matrix and subpopulation assignments
#'
#' Given the biallelic genotypes of `n` individuals, this function returns the `n`-by-`n` kinship matrix such that the kinship estimate between the most distant subpopulations is zero on average (this sets the ancestral population to the most recent common ancestor population).
#'
#' The subpopulation assignments are only used to estimate the baseline kinship (the zero value).
#' If the user wants to re-estimate the kinship matrix using different subpopulation labels,
#' it suffices to rescale it using [rescale_popkin()]
#' (as opposed to starting from the genotypes again, which gives the same answer but more slowly).
#' 
#' @param X Genotype matrix, BEDMatrix object, or a function `X(m)` that returns the genotypes of all individuals at `m` successive locus blocks each time it is called, and `NULL` when no loci are left.
#' If a regular matrix, `X` must have values only in `c(0, 1, 2, NA)`, encoded to count the number of reference alleles at the locus, or `NA` for missing data.
#' @param subpops The length-`n` vector of subpopulation assignments for each individual.
#' If `NULL`, every individual is effectively treated as a different population.
#' @param n Number of individuals (required only when `X` is a function, ignored otherwise).
#' If `n` is missing but `subpops` is not, `n` is taken to be the length of `subpops`.
#' @param loci_on_cols If `TRUE`, `X` has loci on columns and individuals on rows; if `FALSE` (default), loci are on rows and individuals on columns.
#' Has no effect if `X` is a function.
#' If `X` is a BEDMatrix object, `loci_on_cols` is ignored (set automatically to `TRUE` internally).
#' @param mem_factor Proportion of available memory to use loading and processing genotypes.
#' Ignored if `mem_lim` is not `NA`.
#' @param mem_lim Memory limit in GB, used to break up genotype data into chunks for very large datasets.
#' Note memory usage is somewhat underestimated and is not controlled strictly.
#' Default in Linux and Windows is `mem_factor` times the free system memory, otherwise it is 1GB (OSX and other systems).
#' @param want_M If `TRUE`, includes the matrix `M` of non-missing pair counts in the return value, which are sample sizes that can be useful in modeling the variance of estimates.
#' Default `FALSE` is to return the kinship matrix only.
#' @param m_chunk_max Sets the maximum number of loci to process at the time.
#' Actual number of loci loaded may be lower if memory is limiting.
#'
#' @return If `want_M = FALSE`, returns the estimated `n`-by-`n` kinship matrix only.
#' If `X` has names for the individuals, they will be copied to the rows and columns of this kinship matrix.
#' If `want_M = TRUE`, a named list is returned, containing:
#'
#' - `kinship`: the estimated `n`-by-`n` kinship matrix
#' - `M`: the `n`-by-`n` matrix of non-missing pair counts (see `want_M` option).
#'
#' @examples
#' # Construct toy data
#' X <- matrix(
#'     c(0, 1, 2,
#'       1, 0, 1,
#'       1, 0, 2),
#'     nrow = 3,
#'     byrow = TRUE
#' ) # genotype matrix
#' subpops <- c(1,1,2) # subpopulation assignments for individuals
#' 
#' # NOTE: for BED-formatted input, use BEDMatrix!
#' # "file" is path to BED file (excluding .bed extension)
#' ## library(BEDMatrix)
#' ## X <- BEDMatrix(file) # load genotype matrix object
#'
#' kinship <- popkin(X, subpops) # calculate kinship from genotypes and subpopulation labels
#'
#' @seealso
#' [popkin_af()] for coancestry estimation from allele frequency matrices.
#'
#' @export
popkin <- function(
                   X,
                   subpops = NULL,
                   n = NA,
                   loci_on_cols = FALSE,
                   mem_factor = 0.7,
                   mem_lim = NA,
                   want_M = FALSE,
                   m_chunk_max = 1000
                   ) {
    # wrapper around popkin_A combined with subpopulation-based estimation of A_Emin
    if ( missing( X ) )
        stop( 'Genotype matrix `X` is required!' )
    
    # repeat some validations before the hard work... (some are repeated again inside each function, but let's do it sooner)
    # test coherence between subpops and n...
    if (!is.null(subpops)) {
        if (is.na(n)) {
            n <- length(subpops)
        } else if (n != length(subpops)) {
            # if both were specified, they better agree
            stop('the length of `subpops` (', length(subpops), ') disagreed with the input `n` (', n, ')')
        }
        # also compare to X
        if ( !is.function(X) ) {
            # NOTE BEDMatrix is also class 'matrix', so have to test in this order
            if ('BEDMatrix' %in% class(X)) {
                n2 <- nrow(X)
            } else if (is.matrix(X)) {
                if (loci_on_cols) {
                    n2 <- nrow(X)
                } else {
                    n2 <- ncol(X)
                }
            } else
                stop('genotype matrix `X` has unsupported class: ', toString( class(X) ) )
            if (n != n2)
                stop('the length of `subpops` (', n, ') disagrees with the number of individuals in the genotype matrix (', n2, ')')
        }
    }
    
    # actually run code
    
    # this is the main workhorse, estimating the numerators
    obj <- popkin_A(
        X,
        n = n,
        loci_on_cols = loci_on_cols,
        mem_factor = mem_factor,
        mem_lim = mem_lim,
        m_chunk_max = m_chunk_max
    )
    A <- obj$A
    M <- obj$M # no longer needed unless user wants it
    
    # the denominator is a simple average, a scalar shared by all individuals
    A_min <- popkin_A_min_subpops(A, subpops)
    
    # the kinship matrix is this simple ratio
    kinship <- 1 - A / A_min
    
    # figure out what to return
    if ( want_M ) {
        return(
            list(
                kinship = kinship,
                M = M
            )
        )
    } else {
        return( kinship )
    }
}
