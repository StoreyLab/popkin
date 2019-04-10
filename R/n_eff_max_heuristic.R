# Heuristic solution of n_eff_max under non-negative weights

n_eff_max_heuristic <- function(kinship, weights) {
    # die if these are missing
    if (missing(kinship))
        stop('`kinship` matrix is required!')
    if (missing(weights))
        stop('weights from full solution are missing!')
    
    # usually the input weights have at least one negative value, die if this isn't so!
    if (min(weights) >= 0)
        stop('input weights were supposed to have at least one negative value but were all non-negative!')
    
    # run additional validations
    validate_kinship(kinship)

    # is_positive is a vector that indicates which weights are positive (non-zero)
    n <- nrow(kinship)
    is_positive <- rep.int(TRUE, n) # initialize this to include everybody at first
    i <- which.min(weights)
    while (weights[i] < 0) { # loop until we don't have negative cases
        is_positive[i] <- FALSE # set only the minimum weight to false
        weights[i] <- 0 # set this value to zero, it will stay zero
        inverse_kinship <- solve(kinship[is_positive, is_positive]) # consider submatrix only
        weights[is_positive] <- rowSums( inverse_kinship ) # unnormalized weights (this is ok)
        i <- which.min(weights)
    }
    # n_eff is the sum of values
    # NOTE: should hold even when some weights are zero (as long as we're summing the elements of the inverse of the submatrix given above)
    n_eff <- sum( inverse_kinship )
    # normalize weights
    weights <- weights / sum(weights)
    # return
    return(
        list(n_eff = n_eff, weights = weights)
    )
}
