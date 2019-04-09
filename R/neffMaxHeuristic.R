## Heuristic solution of nEffMax under non-negative weights

neffMaxHeuristic <- function(Phi, w) {
    # sanity checks
    if (missing(Phi)) stop('kinship matrix is missing!')
    n <- nrow(Phi)
    if (ncol(Phi) != n) stop('kinship matrix is not square!')
    # in this case we require the weights from the full solution!
    if (missing(w)) stop('weights from full solution are missing!')
    # usually the input weights have at least one negative value, die if this isn't so!
    if (min(w) >= 0) stop('input weights were supposed to have at least one negative value but were all non-negative!')
    
    # isPos is a vector that indicates which weights are positive (non-zero)
    isPos <- rep.int(TRUE, n) # initialize this to include everybody at first
    i <- which.min(w)
    while (w[i] < 0) { # loop until we don't have negative cases
        isPos[i] <- FALSE # set only the minimum weight to false
        w[i] <- 0 # set this value to zero, it will stay zero
        PhiInv <- solve(Phi[isPos,isPos]) # consider submatrix only
        w[isPos] <- rowSums( PhiInv ) # unnormalized weights (this is ok)
        i <- which.min(w)
    }
    # nEff is the sum of values
    # NOTE: should hold even when some weights are zero (as long as we're summing the elements of the inverse of the submatrix given above)
    nEff <- sum(PhiInv)
    
    return( list(neff=nEff, w=w/sum(w)) ) # weights were unnormalized until now!
}
