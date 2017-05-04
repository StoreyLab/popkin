#' Extract inbreeding coefficients from a kinship matrix
#'
#' The kinship matrix \eqn{\Phi^T} contains inbreeding coefficients \eqn{f_j^T} along the diagonal, present as \eqn{\phi_{jj}^T = \frac{1}{2}(1+f_j^T)}{\phi_jj^T = (1+f_j^T)/2}.  This function extracts the vector of \eqn{f_j^T} values from the input \eqn{\Phi^T}.
#' 
#' @param Phi The \eqn{n \times n}{n-by-n} kinship matrix \eqn{\Phi^T}.
#'
#' @return The length-\eqn{n} vector of inbreeding coefficients \eqn{f_j^T} for each individual \eqn{j}.
#'
#' @examples
#' \dontrun{inbr <- inbr(Phi)}
#'
#' @export
inbr <- function(Phi) {
    2 * diag(Phi) - 1  # returns vector of inbreeding coefficients!
}

