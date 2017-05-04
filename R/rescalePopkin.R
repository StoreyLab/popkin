#' Rescale kinship matrix to set a given kinship value to zero.
#'
#' Rescales the input kinship matrix \eqn{\Phi} so that the value \eqn{\phi_{min}} in the original kinship matrix becomes zero, using the formula
#' \deqn{\Phi' = \frac{\Phi - \phi_{min}}{1 - \phi_{min}}.}
#' If subpopulation labels 'subpops' are provided, they are used to estimate \eqn{\phi_{min}} using \code{\link{minAvgSubpops}} internally.
#' If both subpops and phiMin are provided, script stops with a fatal error.
#'
#' @param Phi An \eqn{n \times n} kinship matrix.
#' @param subpops The length-\eqn{n} vector of subpopulation assignments for each individual.
#' @param phiMin A scalar kinship value to define the new zero kinship.
#'
#' @return The rescaled \eqn{n \times n} kinship matrix, with the desired level of relatedness set to zero.
#'
#' @examples
#' \dontrun{
#' ## suppose we first estimate kinship without subpopulations, which will be more biased
#' ## This example assumes input is in BED format and is loaded using BEDMatrix
#' ## "file" is path to BED file (excluding .bed extension)
#' library(BEDMatrix)
#' X <- BEDMatrix(file) # load genotype matrix object
#' Phi <- popkin(X) # calculate kinship from genotypes, WITHOUT subpopulation labels "subpops"
#' ## then we visualize this matrix, figure out a reasonable subpopulation partition "subpops"
#'
#' ## now we can adjust the kinship matrix!
#' Phi2 <- rescalePopkin(Phi, subpops) # direct way, recommended
#' ## prev is faster but otherwise equivalent to re-estimating Phi from scratch with subpops:
#' ## Phi2 <- popkin(X, subpops) 
#'
#' ## can also manually set the level of relatedness phiMin we want to be zero:
#' Phi2 <- rescalePopkin(Phi, phiMin=phiMin)
#' }
#'
#' @export
rescalePopkin <- function(Phi, subpops, phiMin) {
    ## validate inputs
    if (missing(Phi)) {
        stop('Fatal: you must provide a kinship matrix "Phi" to rescale!')
    } else if (class(Phi) != 'matrix') {
        stop('Fatal: input kinship matrix "Phi" must be class "matrix"!')
    }
    if (!missing(subpops)) {
        if (!missing(phiMin)) {
            stop('Fatal: provided both "subpops" and "phiMin"!  Please provide only one.')
        } else {
            phiMin <- minAvgSubpops(Phi, subpops)
        }
    } else if (missing(phiMin)) stop('Fatal: did not provide either of "subpops" or "phiMin".  Please provide exactly one of these.')
    ## finally, perform a simple IBD rescaling
    Phi <- (Phi - phiMin)/(1 - phiMin) # return this matrix!
}
