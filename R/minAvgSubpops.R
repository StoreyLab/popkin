## Estimate the minimum expected value of a matrix \eqn{Y} using subpopulations
##
## This function averages the values of a matrix \eqn{Y} between every subpopulation pair and returns the minimum of these averages.
## The return value can be used to adjust an \eqn{A} matrix to yield the kinship matrix, or rescale a given kinship matrix (see below).
## 
## When \eqn{Y=A} is the \eqn{A} matrix, the output is a stable estimate of \eqn{\hat{A}_{Emin}}.
## When \eqn{Y=\Phi} is a kinship matrix, the output is a stable estimate of the minimum kinship \eqn{\phi_{min}}.
## If the kinship values between the most distant pair of subpopulations is truly zero, then the estimated \eqn{\hat{A}_{Emin}} can be used to adjust \eqn{A} into consistent kinship estimates using \code{\link{getKinshipFromA}}, or the estimated \eqn{\phi_{min}} can be used to rescale the pre-existing kinship matrix \eqn{\Phi} used as input using \code{\link{rescalePopkin}}.
##
## If no subpopulation partition is provided, the function returns the minimum value of \eqn{Y}.
## This default choice may be appropriate in some settings, but is susceptible to bias when there are few SNPs and many pairs of individuals with zero kinship (taking the most extreme estimate is clearly worse than averaging these values would be).
## This default is provided for convenience, to explore the data when a correct choice of subpopulations is not clear, but is not recommended as a final approach.
## 
## @param Y A symmetric \eqn{n \times n} matrix with values between every individual pair, including self comparisons.
## @param subpops A length-\eqn{n} vector of subpopulation assignments for each individual.  If missing, every individual is effectively treated as a different population.
##
## @return The scalar estimate of the minimum expected value of \eqn{Y}.
##
## @examples
## \dontrun{
## ## the recommended form using appropriate subpopulation labels
## AEMinHat <- minAvgSubpops(A, subpops)
## ## a simple default for exploratory analysis
## AEMinHat <- minAvgSubpops(A) # == min(A)
## }
minAvgSubpops <- function(Y, subpops=NULL) {
    if (is.null(subpops)) return(min(Y))
    
    ## check the dimensions
    n <- nrow(Y)
    stopifnot(ncol(Y) == n)
    stopifnot(length(subpops) == n)

    ## get unique subpopulations to navigate
    uniqueSubpops <- sort(unique(subpops))
    K <- length(uniqueSubpops)
    if (K < 2) stop('Error: cannot estimate Y_min with less than two subpopulations (K=',K,')') # we must have at least two subpopulations or the code below fails!
    ## just because, let's store these values in a submatrix (don't strictly need the whole thing, but meh)
    ## initialize the matrix with the maximum value, so this choice doesn't interfere with min-finding
    YBarSubpops <- matrix(max(Y), nrow=K, ncol=K)
    ## now compare pairs: ignore self (is never minimum) and only do each pair once!
    ## NOTE: i and j index subpopulations here, not individuals like in our paper
    for (i in 1:(K-1)) {
        is <- which(subpops == uniqueSubpops[i]) # indexes of individuals from subpopulation i
        for (j in (i+1):K) {
            js <- which(subpops == uniqueSubpops[j]) # indexes of individuals from subpopulation j
            Yij <- mean(Y[is,js]) # compute the desired mean
            YBarSubpops[i,j] <- Yij # stored one way (not both ways)
        }
    }
    YMin <- min(YBarSubpops) # this is the minimum value to return!
    return(YMin) # cool, now return!
}

