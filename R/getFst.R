#' Extract \eqn{F_{ST}} from a kinship matrix or a vector of inbreeding coefficients
#'
#' @param x The vector of inbreeding coefficients, or the kinship matrix (if class(x) is matrix).
#' @param w Weights for individuals to use in calculating \eqn{F_{ST}} (optional, defaults to uniform weights)
#'
#' @return \eqn{F_{ST}}, which is the weighted mean inbreeding coefficient.
#'
#' @examples
#' \dontrun{
#' Fst <- getFst(inbr, w) # use inbreeding vector as input
#' Fst <- getFst(Phi, w) # use kinship matrix as input
#' Fst <- getFst(Phi) # no weights implies uniform weights
#' }
#'
#' @export
getFst <- function(x, w) {
    ## if input is a matrix, let's assume it is the kinship matrix, so extract the inbreeding coefficients first
    if (class(x) == 'matrix') {
        x <- getInbr(x)
    }
    ## now x is a vector of inbreeding coefficients
    if (missing(w)) {
        return( mean(x) ) # no weights implies uniform weights
    } else {
        return( drop( x %*% w ) ) # weighted mean
    }
}

