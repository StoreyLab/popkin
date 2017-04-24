#' Replace kinship diagonal with inbreeding coefficients
#'
#' The usual kinship matrix contains self-kinship values \eqn{(1+f_j)/2} where \eqn{f_j} are inbreeding coefficients.
#' This function returns a modified kinship matrix with the diagonal replaced with the \eqn{f_j} values.
#' This form produces more aesthetically pleasing visualizations, but is not appropriate for modeling (e.g. in GWAS or heritability).
#'
#' @param Phi The kinship matrix with self-kinship values along the diagonal
#'
#' @return The modified kinship matrix, with inbreeding coefficients along the diagonal
#'
#' @examples
#' \dontrun{PhiMod <- getKinshipWithInbrDiag(Phi)}
#'
#' @export
getKinshipWithInbrDiag <- function(Phi) {
    ## returns same kinship matrix but with inbreeding along diagonal instead of self-kinship
    ## these are always better for plots
    diag(Phi) <- getInbr(Phi) # this is the only transformation needed
    Phi # return edited matrix
}

