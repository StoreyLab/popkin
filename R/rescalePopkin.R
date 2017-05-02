#' Rescale kinship matrix to set a given kinship value to zero.
#'
#' Rescales the input kinship matrix \eqn{\Phi} so that the value \eqn{\phi_{min}} in the original kinship matrix becomes zero, using the formula
#' \deqn{\Phi' = \frac{\Phi - \phi_{min}}{1 - \phi_{min}}.}
#'
#' @param Phi An \eqn{n \times n} kinship matrix.
#' @param phiMin A scalar kinship value to define the new zero kinship.
#'
#' @return The rescaled \eqn{n \times n} kinship matrix, with the desired level of relatedness set to zero.
#'
#' @examples
#' \dontrun{
#' ## suppose first we estimate the kinship matrix without subpopulations, which is more likely to be biased
#' ## This example assumes input is in BED format and is loaded using BEDMatrix
#' ## "file" is path to BED file (excluding .bed extension)
#' library(BEDMatrix)
#' X <- BEDMatrix(file) # load genotype matrix object
#' Phi <- popkin(X) # calculate kinship from genotypes, WITHOUT subpopulation labels "subpops"
#' ## then we visualize this matrix, figure out a reasonable subpopulation partition "subpops"
#' phiMin <- minAvgSubpops(Phi, subpops) # now we have a more reasonable estimate of the minimum kinship
#' Phi <- rescalePopkin(Phi, phiMin) # set that level of relatedness to zero!
#' }
#'
#' @export
rescalePopkin <- function(Phi, phiMin) {
    ## simple IBD rescaling
    Phi <- (Phi - phiMin)/(1 - phiMin) # return this matrix!
}
