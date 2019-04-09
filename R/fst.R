#' Calculate FST from a population-level kinship matrix or vector of inbreeding coefficients
#'
#' This function simply returns the weighted mean inbreeding coefficient.
#' If weights are \code{NULL} (default), the regular mean is calculated.
#' If a kinship matrix is provided, then the inbreeding coefficients are extracted from its diagonal using \code{\link{inbr}} (requires the diagonal to contains self-kinship values (\eqn{\phi_{jj}^T = \frac{1}{2}(1+f_j^T)}{\phi_jj^T = (1+f_j^T)/2}) as \code{\link{popkin}} returns, and not inbreeding coefficients (\eqn{f_j^T}) as \code{\link{inbr_diag}} returns).
#'
#' The returned weighted mean inbreeding coefficient equals the generalized \eqn{F_{ST}}{FST} if all individuals are "locally outbred" (i.e. if the self-relatedness of every individual stems entirely from the population structure rather than due partly to having unusually closely related parents, such as first or second cousins).
#' Note most individuals in population-scale human data are locally outbred.
#' If there are locally-inbred individuals, the returned value will overestimate \eqn{F_{ST}}{FST}.
#' 
#' @param x The vector of inbreeding coefficients, or the kinship matrix if \code{x} is a matrix.
#' @param w Weights for individuals (optional, defaults to uniform weights)
#'
#' @return \eqn{F_{ST}}{FST}
#'
#' @examples
#' # Get FST from a genotype matrix
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
#' # estimate the kinship matrix "Phi" from the genotypes "X"!
#' Phi <- popkin(X, subpops) # calculate kinship from X and optional subpop labels
#' w <- weights_subpops(subpops) # can weigh individuals so subpopulations are balanced
#' Fst <- fst(Phi, w) # use kinship matrix and weights to calculate fst
#' 
#' Fst <- fst(Phi) # no (or NULL) weights implies uniform weights
#'
#' inbr <- inbr(Phi) # if you extracted inbr for some other analysis...
#' Fst <- fst(inbr, w) # ...use this inbreeding vector as input too!
#'
#' @export
fst <- function(x, w=NULL) {
    # validate inputs
    if (missing(x))
        stop('you must provide a kinship matrix or vector of inbreeding coefficients!')
    
    # if input is a matrix, let's assume it is the kinship matrix, so extract the inbreeding coefficients first
    # let inbr validate kinship matrix
    if (is.matrix(x))
        x <- inbr(x)
    
    # now x is a vector of inbreeding coefficients
    if (is.null(w)) {
        return( mean(x) ) # no weights implies uniform weights
    } else {
        if (length(w) != length(x)) {
            # sanity check: make sure these have matching dimensions!
            # note that x is a vector at this point (if it was a matrix earlier, the vector of inbreeding coefficients has been extracted), so this will always work when the data makes sense
            stop('the number of individuals differs between inbreeding vector (', length(x), ') and weight vector (', length(w), ')')
        } else {
            return( drop( x %*% w ) ) # weighted mean
        }
    }
}

