#' Replace kinship diagonal with inbreeding coefficients
#'
#' The usual kinship matrix contains self-kinship values \eqn{\phi_{jj}^T = \frac{1}{2}(1+f_j^T)}{\phi_jj^T = (1+f_j^T)/2} where \eqn{f_j^T} are inbreeding coefficients.
#' This function returns a modified kinship matrix with each \eqn{\phi_{jj}^T}{\phi_jj^T} replaced with \eqn{f_j} (off-diagonal \eqn{j \ne k}{j != k} values stay the same).
#' This form produces more aesthetically pleasing visualizations, but is not appropriate for modeling (e.g. in GWAS or heritability estimation).
#'
#' @param Phi The kinship matrix with self-kinship values along the diagonal
#'
#' @return The modified kinship matrix, with inbreeding coefficients along the diagonal
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
#' ## lastly, replace diagonal of kinship matrix with inbreeding coefficients
#' PhiMod <- inbrDiag(Phi)
#'
#' @export
inbrDiag <- function(Phi) {
    ## returns same kinship matrix but with inbreeding along diagonal instead of self-kinship
    ## these are always better for plots
    diag(Phi) <- inbr(Phi) # this is the only transformation needed
    Phi # return edited matrix
}

