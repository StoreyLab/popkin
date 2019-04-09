## Estimate kinship from \eqn{A} and an estimate of the minimum expected A
##
## Simply rescales the matrix \eqn{A} using the estimate of its minimum expected value \eqn{A_{Emin}} using the formula
## \deqn{\Phi = 1 - \frac{A}{A_{Emin}}.}
##
## @param A The \eqn{n \times n}{n-by-n} A matrix from \code{getA}.
## @param AEMin The estimate of the minimum expected value of A.
##
## @return The estimated \eqn{n \times n}{n-by-n} kinship matrix, with a minimum expected value of zero.
##
## @examples
## ## Construct toy data
## X <- matrix(c(0,1,2,1,0,1,1,0,2), nrow=3, byrow=TRUE) # genotype matrix
## subpops <- c(1,1,2) # subpopulation assignments for individuals
##
## ## NOTE: for BED-formatted input, use BEDMatrix!
## ## "file" is path to BED file (excluding .bed extension)
## # library(BEDMatrix)
## # X <- BEDMatrix(file) # load genotype matrix object
##
## A <- getA(X) # calculate A from genotypes
## AEMin <- min_mean_subpops(A, subpops)
## Phi <- getKinshipFromA(A, AEmin)
getKinshipFromA <- function(A, AEMin) {
    ## implements the simple transformation of A and AEMinHat into PhiHat
    Phi <- 1 - A/AEMin # return this matrix!
}
