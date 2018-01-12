#' Estimate the individual-level pairwise FST matrix
#'
#' This function construct the individual-level pairwise \eqn{F_{ST}}{FST} matrix implied by the input kinship matrix.
#' If the input is the true kinship matrix, the return value corresponds to the true pairwise \eqn{F_{ST}}{FST} matrix.
#' On the other hand, if the input is the estimated kinship returned by \code{\link{popkin}}, then the return value is the pairwise \eqn{F_{ST}}{FST} estimates described in our paper.
#' In all cases the diagonal of the pairwise \eqn{F_{ST}}{FST} matrix is zero by definition.
#'
#' @param Phi The \eqn{n \times n}{n-by-n} kinship matrix
#'
#' @return The \eqn{n \times n}{n-by-n} pairwise \eqn{F_{ST}}{FST} matrix
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
#' ## estimate the kinship matrix "Phi" from the genotypes "X"!
#' Phi <- popkin(X, subpops) # calculate kinship from X and optional subpop labels
#'
#' ## lastly, compute pairwise FST matrix from the kinship matrix
#' pwF <- pwfst(Phi)
#'  
#' @export
pwfst <- function(Phi) {
    ## sanity check
    if (nrow(Phi) != ncol(Phi)) stop('Fatal: input kinship matrix is not square (dims: ', nrow(Phi), ' x ', ncol(Phi), ')')
    ## the below code works best with Phi scaled like a coancestry matrix
    Phi <- inbrDiag(Phi)
    ## extract inbreeding coefficients
    inbrs <- diag(Phi)
    ## construct other things for estimation
    n <- length(inbrs)
    inbrMat <- matrix(inbrs, n, n) # repeats inbreeding coefficients along columns
    ## return the following desired pairwise Fst matrix!
    ( ( inbrMat + t(inbrMat) ) / 2 - Phi ) / (1 - Phi)
}
