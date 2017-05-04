#' Estimate kinship from a genotype matrix and subpopulation assignments
#'
#' Given the biallelic genotypes of \eqn{n} individuals, this function returns the \eqn{n}-by-\eqn{n} kinship matrix \eqn{\Phi} such that the kinship estimate between the most distant subpopulations is zero on average.
#' The subpopulation labels 
#'
#' The matrix X (or the vectors returned by the function X) must have values only in c(0,1,2,NA), encoded to count the number of reference alleles at the locus, or NA for missing data.
#'
#' \code{popkin} is a wrapper function that applies \code{getA}, \code{minAvgSubpops}, and \code{getKinshipFromA}.
#'
#' @param X Genotype matrix, class BEDMatrix object, or a function that returns the genotypes of all individuals at successive loci each time it is called, and NULL when no loci are left.
#' @param subpops The length-\eqn{n} vector of subpopulation assignments for each individual.  If missing, every individual is effectively treated as a different population.
#' @param n Number of individuals (required only when X is a function, ignored otherwise).  If n is missing but subpops is not, n is taken to be the length of the subpops vector.
#' @param lociOnCols If true, X has loci on columns and individuals on rows; if false, loci are on rows and individuals on columns. Has no effect if X is a function.  If X is a BEDMatrix object, lociOnCols=TRUE is set automatically.
#' @param memLim Memory limit in GB used to calculate the "chunk size" (numbers of SNPs). Note memory usage is somewhat underestimated and is not controlled strictly.  Default is 2GB, except in linux it is the free memory in the system times 0.7.
#' @param verbose If true, prints messages to indicate which step is being performed.
#'
#' @return The estimated \eqn{n \times n} kinship matrix.
#'
#' @examples
#' \dontrun{
#' ## This example assumes input is in BED format and is loaded using BEDMatrix
#' ## "file" is path to BED file (excluding .bed extension)
#' library(BEDMatrix)
#' X <- BEDMatrix(file) # load genotype matrix object
#' Phi <- popkin(X, subpops) # calculate kinship from genotypes and subpopulation labels "subpops"
#' }
#'
#' @export
popkin <- function(X, subpops=NULL, n=NA, lociOnCols=FALSE, memLim=NA, verbose=FALSE) {
    ## wrapper around getA combined with subpopulation-based estimation of A_Emin

    ## repeat some validations before the hard work... (some are repeated again inside each function, but let's do it sooner)
    ## test coherence between subpops and n...
    if (!is.null(subpops)) {
        if (is.na(n)) {
            n <- length(subpops)
        } else if (n != length(subpops)) {
            ## if both were specified, they better agree
            stop('Fatal: the length of subpops (', length(subpops), ') disagreed with the input n (', n, ')')
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
            } else stop('Fatal: X has unsupported class: ', class(X))
            if (n != n2) stop('Fatal: the length of subpops (', n, ') disagreed with the number of individuals in the genotype matrix (', n2, ')')
        }
    }
    ## actually run code
    if (verbose) message('Making A...')
    A <- getA(X, n=n, lociOnCols=lociOnCols, memLim=memLim)
    if (verbose) message("Estimating A_Emin using subpopulations...")
    AEMin <- minAvgSubpops(A, subpops)
    if (verbose) message("Transforming A into final kinship matrix...")
    Phi <- getKinshipFromA(A, AEMin)
}
