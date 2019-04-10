# Newton's method implementation of n_eff_max solver under non-negative weights

n_eff_max_newton <- function(kinship, weights = NULL, tol = 1e-10, verbose = FALSE) {
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
    # underlying variables are square root of weights
    v <- sqrt(weights)
    mean_kin <- drop( weights %*% kinship %*% weights ) # for initializations
    # initialize max sol values for now
    weights_max <- weights # initialize weights that give current max n_eff 
    n_eff_max <- 1 / mean_kin # initialize this as the max value!
    # not sure what to set lambda to, but in the final solution it's twice the mean kinship
    lambda <- 2 * mean_kin # this is a scalar!

    step <- rep.int(1, n) # a big vector, to initialize loop below

    # stop when step becomes too small
    while(sum(step^2) > tol) {
        # will need vector of mean kinship values per row
        mean_kin_j <- drop( weights %*% kinship )
        stopifnot( length(mean_kin_j) == n)
        # compute Gradient/4
        G <- v * (mean_kin_j - lambda/2) # v components
        G <- c(G, (1-sum(weights))/4) # add lambda component in the end
        # compute Hessian/4
        H <- 2 * t(kinship * v) * v # this weird element-wise product is needed
        # diagonal has additional terms
        diag(H) <- diag(H) + mean_kin_j
        # add border terms
        H <- rbind(H, -v/2) # lower border
        H <- cbind(H, c(-v/2, 0)) # right border
        # get new step
        step <- solve(H, G)
        # compute new v after this step
        v <- v - step[1:n]
        # new weights
        weights <- v^2
        # new lambda
        lambda <- lambda - step[n+1]
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
