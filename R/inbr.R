#' Extract inbreeding coefficients from a kinship matrix
#'
#' The kinship matrix contains transformed inbreeding coefficients along the diagonal.
#' This function extracts the vector of inbreeding values from the input kinship matrix, by transforming the diagonal using the formula `2 * x - 1`.
#' 
#' @param kinship The `n`-by-`n` kinship matrix.
#'
#' @return The length-`n` vector of inbreeding coefficient for each individual.
#'
#' @examples
#' #########
#' # illustrate the main transformation on a 2x2 kinship matrix:
#' # same inbreeding values for both individuals
#' inbr <- 0.2
#' # corresponding self kinship (diagonal values) for both individuals
#' kinship_self <- (1 + inbr)/2
#' # actual kinship matrix (zero kinship between individuals)
#' kinship <- matrix(c(kinship_self, 0, 0, kinship_self), nrow=2)
#' # expected output of inbr (extracts inbreeding coefficients)
#' inbr_exp <- c(inbr, inbr)
#' # actual output from this function
#' inbr_obs <- inbr(kinship)
#' # verify that they match (up to machine precision)
#' stopifnot( all( abs(inbr_obs - inbr_exp) < .Machine$double.eps ) )
#' 
#' #########
#' # Construct toy data
#' X <- matrix(c(0,1,2,1,0,1,1,0,2), nrow=3, byrow=TRUE) # genotype matrix
#' subpops <- c(1,1,2) # subpopulation assignments for individuals
#' 
#' # NOTE: for BED-formatted input, use BEDMatrix!
#' # "file" is path to BED file (excluding .bed extension)
#' ## library(BEDMatrix)
#' ## X <- BEDMatrix(file) # load genotype matrix object
#'
#' # estimate the kinship matrix from the genotypes "X"!
#' kinship <- popkin(X, subpops) # calculate kinship from X and optional subpop labels
#'
#' # extract inbreeding coefficients from Kinship
#' inbr <- inbr(kinship)
#' 
#' @export
inbr <- function(kinship) {
    # die if this is missing
    if (missing(kinship))
        stop('`kinship` matrix is required!')
    
    # validate inputs
    validate_kinship(kinship)
    
    # returns vector of inbreeding coefficients!
    2 * diag(kinship) - 1
}

