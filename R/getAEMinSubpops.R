#' Estimate the minimum expected value of A using subpopulations
#'
#' @param A The \eqn{n \times n} A matrix from \code{getA}.
#' @param subpops The length-\eqn{n} vector of subpopulation assignments for each individual.
#'
#' @return The estimate of the minimum expected value of A (scalar).
#'
#' @examples
#' \dontrun{AEMinHat <- getAEminSubpops(A, subpops)}
#' 
#' @export
getAEminSubpops <- function(A, subpops) {
    ## assuming A is computed previously, estimate its minimum value and return it
    ## here the minimum value is estimated by averaging submatrices of A between subpopulations.
    ## the minimum of these averages is returned
    ## input: subpops is a vector of subpopulation assignments that matches the dimensions of A

    ## check the dimensions
    n <- nrow(A)
    stopifnot(ncol(A) == n)
    stopifnot(length(subpops) == n)

    ## get unique subpopulations to navigate
    uniqueSubpops <- sort(unique(subpops))
    K <- length(uniqueSubpops)
    if (K < 2) stop('Error: cannot estimate A_Emin with less than two subpopulations (K=',K,')') # we must have at least two subpopulations or the code below fails!
    ## just because, let's store these values in a submatrix (don't strictly need the whole thing, but meh)
    ## NOTE: since A values as computed above are always negative, 0 is never the min! (keeping missing comparisons at zero is ok)
    ABarSubpops <- matrix(0, nrow=K, ncol=K)
    ## now compare pairs: ignore self (is never minimum) and only do each pair once!
    ## NOTE: i and j index subpopulations here, not individuals like in our paper
    for (i in 1:(K-1)) {
        is <- which(subpops == uniqueSubpops[i]) # indexes of individuals from subpopulation i
        for (j in (i+1):K) {
            js <- which(subpops == uniqueSubpops[j]) # indexes of individuals from subpopulation j
            Aij <- mean(A[is,js]) # compute the desired mean
            ABarSubpops[i,j] <- Aij # stored one way (not both ways)
        }
    }
    AEMinHat <- min(ABarSubpops) # this is the minimum value to return!
    ## sanity checks
    stopifnot(AEMinHat < 0) # if it equals zero it's really bad (this goes in a denominator). A zero or larger value is a sure sign something went wrong...
    return(AEMinHat) # cool, now return!
}

