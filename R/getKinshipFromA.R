## Estimate kinship from \eqn{A} and an estimate of the minimum expected A
##
## Simply rescales the matrix \eqn{A} using the estimate of its minimum expected value \eqn{A_{Emin}} using the formula
## \deqn{\Phi = 1 - \frac{A}{A_{Emin}}.}
##
## @param A The \eqn{n \times n} A matrix from \code{getA}.
## @param AEMin The estimate of the minimum expected value of A.
##
## @return The estimated \eqn{n \times n} kinship matrix, with a minimum expected value of zero.
##
## @examples
## \dontrun{Phi <- getKinshipFromA(A, AEmin)}
getKinshipFromA <- function(A, AEMin) {
    ## implements the simple transformation of A and AEMinHat into PhiHat
    Phi <- 1 - A/AEMin # return this matrix!
}
