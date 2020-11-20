#' Estimate the minimum expected value of a matrix `A` using subpopulations
#'
#' This function averages the values of a square matrix `A` between every subpopulation pair and returns the minimum of these averages.
#' The return value can be used to adjust an `A` matrix to yield the kinship matrix.
#' 
#' If no subpopulation partition is provided, the function returns the minimum value of `A`.
#' This default choice may be appropriate in some settings, but is susceptible to bias when there are few loci and many pairs of individuals with zero kinship (taking the most extreme estimate is clearly worse than averaging these values).
#' This default is provided for convenience, to explore the data when a correct choice of subpopulations is not clear, but is not recommended as a final approach.
#' 
#' @param A A symmetric `n`-by-`n` matrix with values between every individual pair, including self comparisons.
#' @param subpops A length-`n` vector of subpopulation assignments for each individual.
#' If missing, every individual is effectively treated as a different population.
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
#' # a simple default for exploratory analysis, equals min( A )
#' A_min_est <- popkin_A_min_subpops( A )
#' stopifnot( A_min_est == min( A ) )
#'
#' @seealso
#' The `\link[popkin_A]` to generate the A matrix normally inputted into this function (`popkin_A_min_subpops`), and `\link[popkin]` is the wrapper function around both of these.
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
        return( min( A ) )
    
    # otherwise check these dimensions too
    if( length( subpops ) != nrow( A ) )
        stop( 'Number of individuals in `subpops` (', length(subpops), ') and `A` (', nrow(A), ') disagree!' )
    
    # get unique subpopulations to navigate
    subpops_unique <- sort( unique( subpops ) )
    K <- length( subpops_unique )
    # we must have at least two subpopulations or the code below fails!
    if ( K < 2 )
        stop( 'Cannot estimate A_min with less than two subpopulations (K = ', K, ')' )
    
    # for simplicity, let's store these values in a submatrix (don't strictly need the whole thing)
    # initialize the matrix with the maximum value, so unfilled data won't interfere with min-finding
    mean_subpops <- matrix(
        max( A ),
        nrow = K,
        ncol = K
    )
    
    # now compare pairs: ignore self (is never minimum in positive-definite matrices) and only do each pair once!
    # NOTE: i and j index subpopulations here, not individuals like in our paper
    for ( i in 1 : ( K - 1 ) ) {
        # indexes of individuals from subpopulation i
        indexes_i <- subpops == subpops_unique[ i ]
        for ( j in ( i + 1 ) : K ) {
            # indexes of individuals from subpopulation j
            indexes_j <- subpops == subpops_unique[ j ]
            # compute desired mean, store one way (not both ways)
            mean_subpops[ i, j ] <- mean( A[ indexes_i, indexes_j] )
        }
    }
    
    # this is the minimum value to return!
    return( min( mean_subpops ) )
}

