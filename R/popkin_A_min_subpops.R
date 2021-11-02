#' Estimate the minimum expected value of a matrix `A` using subpopulations
#'
#' This function averages the values of a square matrix `A` between every subpopulation pair and returns the minimum of these averages.
#' If no subpopulation partition is provided, the function returns the minimum value of `A` excluding the diagonal, to agree when the code treats each individual as a subpopulation.
#' The return value can be used to adjust an `A` matrix to yield the kinship matrix.
#' 
#' @param A A symmetric `n`-by-`n` matrix with values between every individual pair, including self comparisons.
#' @param subpops A length-`n` vector of subpopulation assignments for each individual.
#' If missing, every individual is treated as a different subpopulation.
#'
#' @return The minimum of the average between-subpopulation `A` values, which estimates the minimum expected value of `A`
#'
#' @examples
#' # Construct toy data
#' X <- matrix(c(0,1,2,1,0,1,1,0,2), nrow=3, byrow=TRUE) # genotype matrix
#' subpops <- c(1,1,2) # subpopulation assignments for individuals
#'
#' # NOTE: for BED-formatted input, use BEDMatrix!
#' # "file" is path to BED file (excluding .bed extension)
#' ## library(BEDMatrix)
#' ## X <- BEDMatrix(file) # load genotype matrix object
#'
#' # calculate A from genotypes
#' A <- popkin_A(X)$A
#' 
#' # the recommended form using appropriate subpopulation labels
#' A_min_est <- popkin_A_min_subpops( A, subpops )
#'
#' # this recovers the popkin estimate
#' kinship <- 1 - A / A_min_est
#' stopifnot( kinship == popkin( X, subpops ) )
#' 
#' # a simple default for exploratory analysis, equals min( A ) for correctly-calculated A
#' A_min_est <- popkin_A_min_subpops( A )
#' stopifnot( A_min_est == min( A ) )
#'
#' @seealso
#' [popkin_A()] to generate the `A` matrix usually inputted into this function (`popkin_A_min_subpops`).
#' [popkin()] is the wrapper function around both of these.
#'
#' [avg_kinship_subpops()] for the full matrix of mean kinship values between subpopulations.
#'
#' @export
popkin_A_min_subpops <- function(A, subpops = NULL) {
    if ( missing( A ) )
        stop( '`A` matrix is required!' )

    # check the dimensions, etc
    # NOTE: A must be symmetric, otherwise below loop doesn't work
    validate_kinship(A, name = 'A')
    
    # handle a trivial case
    if ( is.null( subpops ) )
        return( min( A[ lower.tri( A ) ] ) )

    # this does the bulk of the work
    mean_subpops <- avg_kinship_subpops( A, subpops )
    
    # this is the minimum value to return!
    return( min( mean_subpops ) )
}

