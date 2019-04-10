# Estimate the minimum expected value of a matrix \eqn{Y} using subpopulations
#
# This function averages the values of a matrix \eqn{Y} between every subpopulation pair and returns the minimum of these averages.
# The return value can be used to adjust an \eqn{A} matrix to yield the kinship matrix, or rescale a given kinship matrix (see below).
# 
# When \eqn{Y=A} is the \eqn{A} matrix, the output is a stable estimate of \eqn{\hat{A}_{Emin}}.
# When \eqn{Y=\Phi} is a kinship matrix, the output is a stable estimate of the minimum kinship \eqn{\phi_{min}}.
#
# If no subpopulation partition is provided, the function returns the minimum value of \eqn{Y}.
# This default choice may be appropriate in some settings, but is susceptible to bias when there are few SNPs and many pairs of individuals with zero kinship (taking the most extreme estimate is clearly worse than averaging these values would be).
# This default is provided for convenience, to explore the data when a correct choice of subpopulations is not clear, but is not recommended as a final approach.
# 
# @param Y A symmetric \eqn{n \times n}{n-by-n} matrix with values between every individual pair, including self comparisons.
# @param subpops A length-\eqn{n} vector of subpopulation assignments for each individual.  If missing, every individual is effectively treated as a different population.
#
# @return The scalar estimate of the minimum expected value of \eqn{Y}.
#
# @examples
# # Construct toy data
# X <- matrix(c(0,1,2,1,0,1,1,0,2), nrow=3, byrow=TRUE) # genotype matrix
# subpops <- c(1,1,2) # subpopulation assignments for individuals
#
# # NOTE: for BED-formatted input, use BEDMatrix!
# # "file" is path to BED file (excluding .bed extension)
# ## library(BEDMatrix)
# ## X <- BEDMatrix(file) # load genotype matrix object
#
# A <- get_A(X) # calculate A from genotypes
# 
# # the recommended form using appropriate subpopulation labels
# AEMinHat <- min_mean_subpops(A, subpops)
# # a simple default for exploratory analysis
# AEMinHat <- min_mean_subpops(A) # == min(A)
# 
min_mean_subpops <- function(Y, subpops = NULL) {
    # handle a trivial case
    if (is.null(subpops))
        return(min(Y))
    
    # check the dimensions
    n <- nrow(Y)
    stopifnot(ncol(Y) == n)
    stopifnot(length(subpops) == n)

    # get unique subpopulations to navigate
    subpops_unique <- sort(unique(subpops))
    K <- length(subpops_unique)
    # we must have at least two subpopulations or the code below fails!
    if (K < 2)
        stop('Error: cannot estimate Y_min with less than two subpopulations (K = ', K, ')')
    
    # just because, let's store these values in a submatrix (don't strictly need the whole thing, but meh)
    # initialize the matrix with the maximum value, so this choice doesn't interfere with min-finding
    mean_subpops <- matrix(max(Y), nrow = K, ncol = K)
    # now compare pairs: ignore self (is never minimum) and only do each pair once!
    # NOTE: i and j index subpopulations here, not individuals like in our paper
    for (i in 1:(K-1)) {
        is <- which(subpops == subpops_unique[i]) # indexes of individuals from subpopulation i
        for (j in (i+1):K) {
            js <- which(subpops == subpops_unique[j]) # indexes of individuals from subpopulation j
            # compute desired mean, store one way (not both ways)
            mean_subpops[i, j] <- mean( Y[is, js] )
        }
    }

    # this is the minimum value to return!
    min(mean_subpops)
}

