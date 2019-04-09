#' Estimate kinship from a genotype matrix and subpopulation assignments
#'
#' Given the biallelic genotypes of \eqn{n} individuals, this function returns the \eqn{n \times n}{n-by-n} kinship matrix such that the kinship estimate between the most distant subpopulations is zero on average (this sets the ancestral population \eqn{T} to the most recent common ancestor population).
#'
#' The subpopulation assignments are only used to estimate the baseline kinship (the zero value).
#' If the user wants to re-estimate the kinship matrix using different subpopulation labels,
#' it suffices to rescale it using \code{\link{rescale_popkin}}
#' (as opposed to starting from the genotypes again, which gives the same answer less efficiently).
#' 
#' The matrix \eqn{X} must have values only in \code{c(0,1,2,NA)}, encoded to count the number of reference alleles at the locus, or \code{NA} for missing data.
#'
#' @param X Genotype matrix, BEDMatrix object, or a function \eqn{X(m)} that returns the genotypes of all individuals at \eqn{m} successive locus blocks each time it is called, and NULL when no loci are left.
#' @param subpops The length-\eqn{n} vector of subpopulation assignments for each individual.  If missing, every individual is effectively treated as a different population.
#' @param n Number of individuals (required only when \eqn{X} is a function, ignored otherwise).  If \eqn{n} is missing but \code{subpops} is not, \eqn{n} is taken to be the length of \code{subpops}.
#' @param lociOnCols If true, \eqn{X} has loci on columns and individuals on rows; if false (the default), loci are on rows and individuals on columns. Has no effect if \eqn{X} is a function.  If \eqn{X} is a BEDMatrix object, \code{lociOnCols=TRUE} is set automatically.
#' @param memLim Memory limit in GB, used to break up genotype data into chunks for very large datasets. Note memory usage is somewhat underestimated and is not controlled strictly.  Default in Linux and Windows is 70 \% of the free system memory, otherwise it is 1GB (OSX and other systems).
#'
#' @return The estimated \eqn{n \times n}{n-by-n} kinship matrix.
#' If \eqn{X} has names for the individuals, they will be copied to the rows and columns of this kinship matrix.
#'
#' @examples
#' ## Construct toy data
#' X <- matrix(c(0,1,2,1,0,1,1,0,2), nrow=3, byrow=TRUE) # genotype matrix
#' subpops <- c(1,1,2) # subpopulation assignments for individuals
#' 
#' ## NOTE: for BED-formatted input, use BEDMatrix!
#' ## "file" is path to BED file (excluding .bed extension)
#' # library(BEDMatrix)
#' # X <- BEDMatrix(file) # load genotype matrix object
#'
#' kinship <- popkin(X, subpops) # calculate kinship from genotypes and subpopulation labels
#'
#' @export
popkin <- function(X, subpops = NULL, n = NA, lociOnCols = FALSE, memLim = NA) {
    ## wrapper around get_A combined with subpopulation-based estimation of A_Emin

    ## repeat some validations before the hard work... (some are repeated again inside each function, but let's do it sooner)
    ## test coherence between subpops and n...
    if (!is.null(subpops)) {
        if (is.na(n)) {
            n <- length(subpops)
        } else if (n != length(subpops)) {
            ## if both were specified, they better agree
            stop('the length of subpops (', length(subpops), ') disagreed with the input n (', n, ')')
        }
        ## also compare to X
        if (class(X) != 'function') {
            if (class(X) == 'BEDMatrix') {
                n2 <- nrow(X)
            } else if (class(X) == 'matrix') {
                if (lociOnCols) {
                    n2 <- nrow(X)
                } else {
                    n2 <- ncol(X)
                }
            } else stop('X has unsupported class: ', class(X))
            if (n != n2) stop('the length of subpops (', n, ') disagreed with the number of individuals in the genotype matrix (', n2, ')')
        }
    }
    # actually run code
    # this is the main workhorse, estimating the numerators
    A <- get_A(X, n = n, lociOnCols = lociOnCols, memLim = memLim)
    # the denominator is a simple average, a scalar shared by all individuals
    AEMin <- min_mean_subpops(A, subpops)
    # the kinship matrix is this simple ratio
    kinship <- 1 - A / AEMin # return this matrix!
}
