# Gradient descent implementation of n_eff_max solver under non-negative weights

n_eff_max_gradient <- function(kinship, weights = NULL, tol = 1e-10, verbose = FALSE) {
    # die if this is missing
    if (missing(kinship))
        stop('`kinship` matrix is required!')
    
    # run additional validations
    validate_kinship(kinship)

    # default: use uniform weights as initial guess (best place to start from spatially)
    n <- nrow(kinship)
    if (is.null(weights)) {
        weights <- rep.int(1/n, n)
    } else {
        if (any(weights<0))
            stop('initial weights must be non-negative!')
        weights <- weights / sum(weights) # renormalize for good measure
    }
    
    # underlying variables are square root of w
    v <- sqrt(weights)
    mean_kin <- drop( weights %*% kinship %*% weights ) # for initializations
    # initialize max sol values for now
    weights_max <- weights # initialize weights that give current max n_eff 
    n_eff_max <- 1 / mean_kin # initialize this as the max value!
    
    step <- rep.int(1, n) # a big vector, to initialize loop below
    
    # stop when step becomes too small
    while(sum(step^2) > tol) {
        # will need vector of mean kinship values per row
        mean_kin_j <- drop( weights %*% kinship )
        # compute gradient
        D <- v * (mean_kin_j - mean_kin) # v components, with lambda=2*mean_kin approx
        # compute pseudo-optimal scaling factor (from a linear approx of direct optimization, plus assumption that previous weights already summed to one, which we ensure is the case)
        alpha <- (
            mean_kin^2 - drop(weights %*% mean_kin_j^2)
        ) / (
            drop(weights %*% mean_kin_j^3)
            - 7 * mean_kin * drop(weights %*% mean_kin_j^2)
            + 4 * mean_kin^3
            + 2 * drop( (weights * mean_kin_j) %*% kinship %*% (weights * mean_kin_j) )
        )
        if (verbose && alpha > 0)
            message('Alpha reversed gradient!')
        # get new step
        step <- alpha * D
        # compute new v after this step
        v <- v + step
        # new weights
        weights <- v^2
        # compute final things of interest
        weights <- weights / sum(weights) # normalize for good measure
        v <- sqrt(weights) # rebuild consistent sol
        mean_kin <- drop( weights %*% kinship %*% weights )
        n_eff <- 1 / mean_kin
        # compare to current max
        if (n_eff > n_eff_max) {
            # reset solution
            weights_max <- weights
            n_eff_max <- n_eff
            if (verbose)
                message('n_eff_max: ', n_eff_max)
        }
    }

    return(
        list(n_eff = n_eff_max, weights = weights_max)
    )
}
