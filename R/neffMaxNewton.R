## Newton's method implementation of nEffMax solver under non-negative weights

neffMaxNewton <- function(Phi, w=NULL, tol=1e-10, verbose=FALSE) {
    # sanity checks
    if (missing(Phi)) stop('kinship matrix is missing!')
    n <- nrow(Phi)
    if (ncol(Phi) != n) stop('kinship matrix is not square!')
    
    # default: use uniform weights as initial guess (best place to start from spatially)
    if (is.null(w)) {
        w <- rep.int(1/n, n)
    } else {
        if (any(w<0)) stop('initial weights must be non-negative!')
        w <- w/sum(w) # renormalize for good measure
    }
    # underlying variables are square root of w
    v <- sqrt(w)
    phiBar <- drop(w %*% Phi %*% w) # for initializations
    # initialize max sol values for now
    wMax <- w # initialize weights that give current max nEff 
    nEffMax <- 1/phiBar # initialize this as the max value!
    # not sure what to set lambda to, but in the final solution it's twice the mean kinship
    lambda <- 2 * phiBar # this is a scalar!

    step <- rep.int(1, n) # a big vector, to initialize loop below

    # stop when step becomes too small
    while(sum(step^2) > tol) {
        # will need vector of mean kinship values per row
        phiBarJ <- drop(w %*% Phi)
        stopifnot( length(phiBarJ) == n)
        # compute Gradient/4
        G <- v * (phiBarJ - lambda/2) # v components
        G <- c(G, (1-sum(w))/4) # add lambda component in the end
        # compute Hessian/4
        H <- 2 * t(Phi * v) * v # this weird element-wise product is needed
        # diagonal has additional terms
        diag(H) <- diag(H) + phiBarJ
        # add border terms
        H <- rbind(H, -v/2) # lower border
        H <- cbind(H, c(-v/2, 0)) # right border
        # get new step
        step <- solve(H, G)
        # compute new v after this step
        v <- v - step[1:n]
        # new weights
        w <- v^2
        # new lambda
        lambda <- lambda - step[n+1]
        # compute final things of interest
        w <- w/sum(w) # normalize for good measure
        v <- sqrt(w) # rebuild consistent sol
        phiBar <- drop(w %*% Phi %*% w)
        nEff <- 1/phiBar
        # compare to current max
        if (nEff > nEffMax) {
            # reset solution
            wMax <- w
            nEffMax <- nEff
            if (verbose) message('nEffMax: ', nEffMax)
        }
    }

    return( list(neff=nEffMax, w=wMax) )
}
