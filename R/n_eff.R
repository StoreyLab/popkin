#' Calculates the effective sample size of the data
#'
#' The effective sample size `n_eff` is the equivalent number of independent haplotypes that gives the same variance as that observed under the given population.
#' The variance in question is for the weighted sample mean ancestral allele frequency estimator.
#' It follows that `n_eff` equals the inverse of the weighted mean kinship.
#' If `max = TRUE`, a calculation is performed that implicitly uses optimal weights which maximize `n_eff`, which equals the sum of the elements of the inverse kinship matrix.
#' However, if `nonneg = TRUE` and if the above solution has negative weights (common), optimal non-negative weights are found instead (there are three algorithms available, see `algo`).
#' If `max = FALSE`, then the input weights are used in this calculation, and if weights are `NULL`, uniform weights are used.
#' 
#' The maximum `n_eff` possible is `2 * n`, where `n` is the number of individuals; this value is attained only when all haplotypes are independent (a completely unstructured population in Hardy-Weinberg equilibrium).
#' The minimum `n_eff` possible is 1, which is attained in an extremely structured population with FST of 1, where every individual has exactly the same haplotype at every locus (no heterozygotes).
#' Moreover, for `K` extremely-differentiated subpopulations (FST = 1 per subpopulation) `n_eff = K`.
#' In this way, `n_eff` is smaller than the ideal value of `2 * n` depending on the amount of kinship (covariance) in the data.
#' 
#' Occasionally, depending on the quality of the input kinship matrix, the estimated `n_eff` may be outside the theoretical `\[1, 2*n\]` range, in which case the return value is set to the closest boundary value.
#' The quality of the results depends on the success of matrix inversion (which for numerical reasons may incorrectly contain negative eigenvalues, for example) or of the gradient optimization.
#'
#' @param kinship An `n`-by-`n` kinship matrix.
#' @param max If `TRUE`, returns the maximum `n_eff` value among those computed using all possible vectors of weights that sum to one (and which are additionally non-negative if `nonneg = TRUE`).
#' If `FALSE`, `n_eff` is computed using the specific weight vector provided.
#' @param weights Weights for individuals (optional).
#' If `NULL`, uniform weights are used.
#' This parameter is ignored if `max = TRUE`.
#' @param nonneg If `TRUE` (default) and `max = TRUE`, non-negative weights that maximize `n_eff` are found.
#' See `algo`.
#' This has no effect if `max = FALSE`.
#' @param algo Algorithm for finding optimal non-negative weights (applicable only if `nonneg = TRUE` and `max = TRUE` and the weights found by matrix inversion are non-negative).
#' May be abbreviated.
#' If "gradient" (default), an optimized gradient descent algorithm is used (fastest; recommended).
#' If "newton", the exact multivariate newton's Method is used (slowest since `(n+1)`-by-`(n+1)` Hessian matrix needs to be inverted at every iteration; use if possible to confirm that "gradient" gives the best answer).
#' If "heuristic", if the optimal solution by the inverse matrix method contains negative weights, the most negative weight in an iteration is forced to be zero in all subsequent iterations and the rest of the weights are solved for using the inverse matrix method, repeating until all resulting weights are non-negative (also slow, since inversion of large matrices is required; least likely to find optimal solution).
#' @param tol Tolerance parameter for "gradient" and "newton" algorithms.
#' The algorithms converge when the norm of the step vector is smaller than this tolerance value. 
#'
#' @return A list containing `n_eff` and `weights` (optimal weights if `max = TRUE`, input weights otherwise).
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
