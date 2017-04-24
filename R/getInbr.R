#' Extract inbreeding coefficients from a kinship matrix
#'
#' @param Phi The kinship matrix, with self-kinship values along the diagonal
#'
#' @return The length-\eqn{n} vector of inbreeding coefficients for each individual.
#'
#' @examples
#' \dontrun{inbr <- getInbr(Phi)}
#'
#' @export
getInbr <- function(Phi) {
    2 * diag(Phi) - 1  # returns vector of inbreeding coefficients!
}

