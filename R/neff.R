#' Calculates the effective sample size of the data
#'
#' The effective sample size (\eqn{n_{eff}}{n_eff}) is the equivalent number of independent haplotypes that gives the same variance as that observed under the given population.
#' The variance in question is for the weighted sample mean ancestral allele frequency estimator.
#' It follows that \eqn{n_{eff}}{n_eff} equals the inverse of the weighted mean kinship.
#' If \code{max=TRUE}, a calculation is performed that implicitly uses optimal weights which maximize \eqn{n_{eff}}{n_eff}: here \eqn{n_{eff}}{n_eff} equals the sum of the elements of the inverse kinship matrix.
#' However, if \code{nonneg=TRUE} and if the above solution has negative weights (common), optimal non-negative weights are found instead (there are three algorithms available, see \code{algo}).
#' If \code{max=FALSE}, then the input weights are used in this calculation, and if weights are \code{NULL}, uniform weights are used.
#' 
#' The maximum \eqn{n_{eff}}{n_eff} possible is \eqn{2n}, where \eqn{n} is the number of individuals; this value is attained only when all haplotypes are independent (a completely unstructured population in Hardy-Weinberg equilibrium).
#' The minimum \eqn{n_{eff}}{n_eff} possible is 1, which is attained in an extremely structured population with \eqn{F_{ST}}{FST} of 1, where every individual has exactly the same haplotype at every locus (no heterozygotes).
#' Moreover, for \eqn{K} extremely-differentiated subpopulations (\eqn{F_{ST}}{FST}=1 per subpopulation) \eqn{n_{eff}}{n_eff} equals \eqn{K}.
#' In this way, \eqn{n_{eff}}{n_eff} is smaller than the ideal value of \eqn{2n} depending on the amount of kinship (covariance) in the data.
#'
#' @param Phi An \eqn{n \times n}{n-by-n} kinship matrix.
#' @param max If \code{TRUE}, returns the maximum \eqn{n_{eff}}{n_eff} value among those computed using all possible vectors of weights that sum to one (and which are additionally non-negative if \code{nonneg=TRUE}).  If \code{FALSE}, \eqn{n_{eff}}{n_eff} is computed using the specific weight vector provided.
#' @param w Weights for individuals (optional). If \code{NULL}, uniform weights are used.  This parameter is ignored if \code{max=TRUE}.
#' @param retW If \code{TRUE} (default), the (input or optimal) weights are returned in addition to \eqn{n_{eff}}{n_eff} (see return value below).  This is highly recommended for troubleshooting when \code{max=TRUE}, as optimal weights may be zero or negative.
#' @param nonneg If \code{TRUE} (default) and \code{max=TRUE}, non-negative weights that maximize \eqn{n_{eff}}{n_eff} are found.  See \code{algo}.  This has no effect if \code{max=FALSE}.
#' @param algo Algorithm for finding optimal non-negative weights (applicable only if \code{nonneg=TRUE} and \code{max=TRUE} and the weights found by matrix inversion are non-negative).
#' If "Gradient" (default), an optimized gradient descent algorithm is used (fastest; recommended).
#' If "Newton", the exact multivariate Newton's Method is used (slowest since \eqn{(n+1) \times (n+1)}{(n+1)-by-(n+1)} Hessian matrix needs to be inverted at every iteration; use if possible to confirm that "Gradient" gives the best answer).
#' If "Heuristic", if the optimal solution by the inverse matrix method contains negative weights, the most negative weight in an iteration is forced to be zero in all subsequent iterations and the rest of the weights are solved for using the inverse matrix method, repeating until all resulting weights are non-negative (also slow, since inversion of large matrices is required; least likely to find optimal solution).
#' @param tol Tolerance parameter for "Gradient" and "Newton" algorithms. The algorithms converge when the norm of the step vector is smaller than this tolerance value. 
#'
#' @return If \code{retW=TRUE}, a list containing \eqn{n_{eff}}{n_eff} and the (input [if \code{max=FALSE}] or optimal [if \code{max=TRUE}]) weights are returned.  Otherwise only \eqn{n_{eff}}{n_eff} is returned.
#'
#' @examples
#' ## Get nEff from a genotype matrix
#' 
#' ## Construct toy data
#' X <- matrix(c(0,1,2,1,0,1,1,0,2), nrow=3, byrow=TRUE) # genotype matrix
#' subpops <- c(1,1,2) # subpopulation assignments for individuals
#' 
#' ## NOTE: for BED-formatted input, use BEDMatrix!
#' ## "file" is path to BED file (excluding .bed extension)
#' # library(BEDMatrix)
#' # X <- BEDMatrix(file) # load genotype matrix object
#'
#' ## estimate the kinship matrix "Phi" from the genotypes "X"!
#' Phi <- popkin(X, subpops) # calculate kinship from X and optional subpop labels
#' w <- weightsSubpops(subpops) # can weigh individuals so subpopulations are balanced
#'
#' # use kinship matrix to calculate nEff
#' # default mode returns maximum nEff possible across all non-negative weights that sum to one
#' # also returns the weights that were optimal
#' obj <- neff(Phi)
#' nEffMax <- obj$neff
#' wMax <- obj$w
#'
#' # version that uses weights provided
#' obj <- neff(Phi, max=FALSE, w=w)
#' nEffW <- obj$neff
#' w <- obj$w # returns input weights renormalized for good measure
#' 
#' # no (or NULL) weights implies uniform weights
#' obj <- neff(Phi, max=FALSE)
#' nEffU <- obj$neff
#' w <- obj$w # uniform weights
#'
#' # get nEff only, disregard weights used
#' nEffMax <- neff(Phi, retW=FALSE)
#'
#' @export
neff <- function (Phi, max=TRUE, w=NULL, retW=TRUE, nonneg=TRUE, algo=c('Gradient', 'Newton', 'Heuristic'), tol=1e-10) {
    # make sure given choice matches one of the only three choices
    algo <- match.arg(algo)
    # sanity checks
    if (missing(Phi)) 
        stop('Fatal: the input kinship matrix is missing!')
    if (class(Phi) != 'matrix')
        stop('Fatal: the input kinship is not an R matrix!')
    n <- nrow(Phi)
    if (ncol(Phi) != n)
        stop('Fatal: the input kinship matrix is not a square matrix!')
    
    # main behavior depends on the value of the weights
    # the NULL default is to use optimal weights, yielding the maximum value of neff possible
    if (max) {
        # invert matrix
        PhiInv <- solve(Phi)
        # compute optimal weights for this problem (unnormalized, they get normalized just before return!)
        # this is not needed to get nEff but it's a good troubleshooting tool
        w <- rowSums(PhiInv)
        if (min(w) < 0 && nonneg) {
            # NOTE: all these algorithms return complete solutions (neff,w), with fully normalized weights
            if (algo == 'Gradient') {
                obj <- neffMaxGradient(Phi, tol=tol)
            } else if (algo == 'Newton') {
                obj <- neffMaxNewton(Phi, tol=tol)
            } else if (algo == 'Heuristic') {
                obj <- neffMaxHeuristic(Phi, w)
            } else stop('Fatal: algorithm "', algo, '" not implemented!')
        } else {
            # complete calculation of standard case!
            # nEff is the sum of values
            nEff <- sum(PhiInv)
            # complete main object to return
            obj <- list(neff=nEff, w=w/sum(w)) # normalize weights here!
        }
    } else {
        # compute mean kinship
        phiBar <- meanKin(Phi, w=w)
        # desired value is inverse of this
        nEff <- 1/phiBar
        # complete main object to return
        obj <- list(neff=nEff, w=w/sum(w)) # normalize weights here!
    }
    # return!
    if (retW) {
        return(obj)
    } else {
        return(obj$neff)
    }
}

# only used internally by neff
# NOTES:
# - allows zero or negative weights
# - renormalizes weights internally for good measure
meanKin <- function(Phi, w=NULL) {
    if (is.null(w)) {
        # this means use weights but use default uniform weights
        return ( mean(Phi) ) # this is the value to return
    } else {
        # weights should be defined, let's make sure they're fine...
        n <- nrow(Phi) # Phi was validated earlier inside neff, skip validations here
        if (n != length(w))
            stop('Fatal: number of individuals in weights and kinship matrix differ!')
        w <- w/sum(w) # force normalization here, don't check if it was already ok
        phiBar <- drop( w %*% Phi %*% w )
        return (phiBar)
    }
}
