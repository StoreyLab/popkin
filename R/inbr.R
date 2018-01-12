#' Extract inbreeding coefficients from a kinship matrix
#'
#' The kinship matrix \eqn{\Phi^T} contains inbreeding coefficients \eqn{f_j^T} along the diagonal, present as \eqn{\phi_{jj}^T = \frac{1}{2}(1+f_j^T)}{\phi_jj^T = (1+f_j^T)/2}.  This function extracts the vector of \eqn{f_j^T} values from the input \eqn{\Phi^T}.
#' 
#' @param Phi The \eqn{n \times n}{n-by-n} kinship matrix \eqn{\Phi^T}.
#'
#' @return The length-\eqn{n} vector of inbreeding coefficients \eqn{f_j^T} for each individual \eqn{j}.
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
#' ## extract inbreeding coefficients from Phi
#' inbr <- inbr(Phi)
#' 
#' @export
inbr <- function(Phi) {
    2 * diag(Phi) - 1  # returns vector of inbreeding coefficients!
}

