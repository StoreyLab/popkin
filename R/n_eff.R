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
#' Occasionally, depending on the quality of the input kinship matrix, the estimated \eqn{n_{eff}}{n_eff} may be outside the theoretical \[1, 2n\] range, in which case the return value is set to the closest boundary value.
#' The quality of the results depends on the success of matrix inversion (which for numerical reasons may incorrectly contain negative eigenvalues, for example) or of the gradient optimization.
#'
#' @param kinship An \eqn{n \times n}{n-by-n} kinship matrix.
#' @param max If \code{TRUE}, returns the maximum \eqn{n_{eff}}{n_eff} value among those computed using all possible vectors of weights that sum to one (and which are additionally non-negative if \code{nonneg=TRUE}).  If \code{FALSE}, \eqn{n_{eff}}{n_eff} is computed using the specific weight vector provided.
#' @param weights Weights for individuals (optional). If \code{NULL}, uniform weights are used.  This parameter is ignored if \code{max=TRUE}.
#' @param nonneg If \code{TRUE} (default) and \code{max=TRUE}, non-negative weights that maximize \eqn{n_{eff}}{n_eff} are found.  See \code{algo}.  This has no effect if \code{max=FALSE}.
#' @param algo Algorithm for finding optimal non-negative weights (applicable only if \code{nonneg=TRUE} and \code{max=TRUE} and the weights found by matrix inversion are non-negative).
#' May be abbreviated.
#' If "gradient" (default), an optimized gradient descent algorithm is used (fastest; recommended).
#' If "newton", the exact multivariate newton's Method is used (slowest since \eqn{(n+1) \times (n+1)}{(n+1)-by-(n+1)} Hessian matrix needs to be inverted at every iteration; use if possible to confirm that "gradient" gives the best answer).
#' If "heuristic", if the optimal solution by the inverse matrix method contains negative weights, the most negative weight in an iteration is forced to be zero in all subsequent iterations and the rest of the weights are solved for using the inverse matrix method, repeating until all resulting weights are non-negative (also slow, since inversion of large matrices is required; least likely to find optimal solution).
#' @param tol Tolerance parameter for "gradient" and "newton" algorithms. The algorithms converge when the norm of the step vector is smaller than this tolerance value. 
#'
#' @return A list containing \code{n_eff} and \code{weights} (optimal weights if \code{max = TRUE}, input weights otherwise).
#'
#' @examples
#' # Get n_eff from a genotype matrix
#' 
#' # Construct toy data
#' X <- matrix(c(0,1,2,1,0,1,1,0,2), nrow=3, byrow=TRUE) # genotype matrix
#' subpops <- c(1,1,2) # subpopulation assignments for individuals
#' 
#' # NOTE: for BED-formatted input, use BEDMatrix!
#' # "file" is path to BED file (excluding .bed extension)
#' ## library(BEDMatrix)
#' ## X <- BEDMatrix(file) # load genotype matrix object
#'
#' # estimate the kinship matrix "kinship" from the genotypes "X"!
#' kinship <- popkin(X, subpops) # calculate kinship from X and optional subpop labels
#' weights <- weights_subpops(subpops) # can weigh individuals so subpopulations are balanced
#'
#' # use kinship matrix to calculate n_eff
#' # default mode returns maximum n_eff possible across all non-negative weights that sum to one
#' # also returns the weights that were optimal
#' obj <- n_eff(kinship)
#' n_eff_max <- obj$n_eff
#' w_max <- obj$weights
#'
#' # version that uses weights provided
#' obj <- n_eff(kinship, max = FALSE, weights = weights)
#' n_eff_w <- obj$n_eff
#' w <- obj$weights # returns input weights renormalized for good measure
#' 
#' # no (or NULL) weights implies uniform weights
#' obj <- n_eff(kinship, max = FALSE)
#' n_eff_u <- obj$n_eff
#' w <- obj$weights # uniform weights
#'
#' @export
n_eff <- function (kinship, max = TRUE, weights = NULL, nonneg = TRUE, algo = c('gradient', 'newton', 'heuristic'), tol = 1e-10) {
    # make sure given choice matches one of the only three choices
    algo <- match.arg(algo)
    
    # die if this is missing
    if (missing(kinship))
        stop('`kinship` matrix is required!')
    
    # run additional validations
    validate_kinship(kinship)

    # data dimensions (used for a validation below)
    n_ind <- ncol( kinship )

    # main behavior depends on the value of the weights
    # the NULL default is to use optimal weights, yielding the maximum value of n_eff possible
    if (max) {
        # invert matrix
        inverse_kinship <- solve(kinship)
        # compute optimal weights for this problem (unnormalized, they get normalized just before return!)
        # this is not needed to get n_eff but it's a good troubleshooting tool
        weights <- rowSums( inverse_kinship )
        if (min(weights) < 0 && nonneg) {
            # NOTE: all these algorithms return complete solutions (n_eff, weights), with fully normalized weights
            if (algo == 'gradient') {
                obj <- n_eff_max_gradient(kinship, tol = tol)
            } else if (algo == 'newton') {
                obj <- n_eff_max_newton(kinship, tol = tol)
            } else if (algo == 'heuristic') {
                obj <- n_eff_max_heuristic(kinship, weights)
            } else
                stop('algorithm not implemented: ', algo)
        } else {
            # complete calculation of standard case!
            # n_eff is the sum of values
            n_eff <- sum( inverse_kinship )
            # normalize weights
            weights <- weights / sum(weights)
            # complete main object to return
            obj <- list(n_eff = n_eff, weights = weights)
        }
    } else {
        # compute mean kinship
        mean_kin <- mean_kinship(kinship, weights)
        # desired value is inverse of this
        n_eff <- 1 / mean_kin
        # normalize weights
        weights <- weights / sum(weights)
        # complete main object to return
        obj <- list(n_eff = n_eff, weights = weights)
    }
    
    # check and fix range of n_eff, if needed
    if ( obj$n_eff < 1 ) {
        obj$n_eff <- 1
    } else if ( obj$n_eff > 2 * n_ind ) {
        obj$n_eff <- 2 * n_ind
    }

    # return!
    return(obj)
}

# stick deprecated function name here

#' @title Calculates the effective sample size of the data
#' @description Calculates the effective sample size of the data
#'
#' @param kinship A kinship matrix
#' @param max Return maximum \eqn{n_{eff}}{n_eff}; otherwise compute \eqn{n_{eff}}{n_eff} using provided weights.
#' @param w Weights for individuals.
#' @param retW Return weights along with \eqn{n_{eff}}{n_eff}.
#' @param nonneg Optimize constrained to non-negative weights.
#' @param algo Algorithm for finding optimal non-negative weights: gradient, newton, or heuristic.
#' @param tol Tolerance parameter for "gradient" and "newton" algorithms.
#'
#' @return A list containing \eqn{n_{eff}}{n_eff} and the weights.
#'
#' @name neff-deprecated
#' @usage neff(kinship, max = TRUE, w = NULL, retW = TRUE, nonneg = TRUE,
#' algo = c('gradient', 'newton', 'heuristic'), tol = 1e-10)
#' @seealso \code{\link{popkin-deprecated}}
#' @keywords internal
NULL

#' @rdname popkin-deprecated
#' @section \code{neff}:
#' For \code{neff}, use \code{\link{n_eff}}.
#'
#' @export
neff <- function(kinship, max = TRUE, w = NULL, retW = TRUE, nonneg = TRUE, algo = c('gradient', 'newton', 'heuristic'), tol = 1e-10) {
    # mark as deprecated
    .Deprecated('n_eff')
    # return as usual, to not break things just yet
    obj <- n_eff(kinship, max = max, weights = w, nonneg = nonneg, algo = algo, tol = tol)
    # obsolete behavior (new function has retW=TRUE hardcoded)
    if (retW) {
        return (obj)
    } else {
        warning('`retW = FALSE` option is deprecated (no longer available in n_eff)')
        return (obj$n_eff)
    }
}
