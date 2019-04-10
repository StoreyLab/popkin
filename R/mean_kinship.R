#' Calculate the weighted mean kinship
#'
#' This function computes a particular weighted mean kinship that arises in the context of kinship and \eqn{F_{ST}}{FST} estimators and in the definition of the effective sample size.
#' This function allows for weights to be zero or even negative, but they are internally normalized to sum to one.
#'
#' @param kinship The kinship matrix
#' @param weights Weights for individuals (optional).
#' If \code{NULL} (default), uniform weights are used.
#'
#' @return The weighted mean kinship matrix, equivalent to \code{ drop( weights \%*\% kinship \%*\% weights ) } for normalized weights (which sum to one).
#'
#' @examples
#' # construct a dummy kinship matrix
#' kinship <- matrix(c(0.5, 0, 0, 0.6), nrow=2)
#' # this is the ordinary mean
#' mean_kinship(kinship)
#' # weighted mean with twice as much weight on the second individual
#' # (weights are internally normalized to sum to one)
#' weights <- c(1, 2)
#' mean_kinship(kinship, weights)
#' 
#' @export
mean_kinship <- function(kinship, weights = NULL) {
    # die if this is missing
    if (missing(kinship))
        stop('`kinship` matrix is required!')
    
    # run additional validations
    validate_kinship(kinship)

    if (is.null(weights)) {
        # this means use weights but use default uniform weights
        return ( mean(kinship) ) # this is the value to return
    } else {
        # check dimensions
        n <- nrow(kinship)
        n2 <- length(weights)
        if (n != n2)
            stop('number of individuals in `kinship` and `weights` differ: ', n , ' != ', n2)

        # actual computations
        weights <- weights / sum(weights) # force normalization here, don't check if it was already ok
        mean_kin <- drop( weights %*% kinship %*% weights )
        return (mean_kin)
    }
}
