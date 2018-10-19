## Gradient descent implementation of nEffMax solver under non-negative weights

neffMaxGradient <- function(Phi, w=NULL, tol=1e-10, verbose=FALSE) {
    # sanity checks
    if (missing(Phi)) stop('Fatal: kinship matrix is missing!')
    n <- nrow(Phi)
    if (ncol(Phi) != n) stop('Fatal: kinship matrix is not square!')
    
    # default: use uniform weights as initial guess (best place to start from spatially)
    if (is.null(w)) {
        w <- rep.int(1/n, n)
    } else {
        if (any(w<0)) stop('Fatal: initial weights must be non-negative!')
        w <- w/sum(w) # renormalize for good measure
    }
    # underlying variables are square root of w
    v <- sqrt(w)
    phiBar <- drop(w %*% Phi %*% w) # for initializations
    # initialize max sol values for now
    wMax <- w # initialize weights that give current max nEff 
    nEffMax <- 1/phiBar # initialize this as the max value!
    
    step <- rep.int(1, n) # a big vector, to initialize loop below
    
    # stop when step becomes too small
    while(sum(step^2) > tol) {
        # will need vector of mean kinship values per row
        phiBarJ <- drop(w %*% Phi)
        # compute gradient
        D <- v * (phiBarJ - phiBar) # v components, with lambda=2*phiBar approx
        # compute pseudo-optimal scaling factor (from a linear approx of direct optimization, plus assumption that previous weights already summed to one, which we ensure is the case)
        alpha <- (
            phiBar^2 - drop(w %*% phiBarJ^2)
        ) / (
            drop(w %*% phiBarJ^3)
            - 7 * phiBar * drop(w %*% phiBarJ^2)
            + 4 * phiBar^3
            + 2 * drop( (w * phiBarJ) %*% Phi %*% (w * phiBarJ) )
        )
        if (verbose && alpha > 0) message('Alpha reversed gradient!')
        # get new step
        step <- alpha * D
        # compute new v after this step
        v <- v + step
        # new weights
        w <- v^2
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
