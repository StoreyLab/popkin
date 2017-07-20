#' @useDynLib popkin
#' @importFrom Rcpp sourceCpp
NULL

## Compute A matrix from genotypes
##
## Given the biallelic genotypes of \eqn{n} individuals, this function returns the \eqn{n}-by-\eqn{n} matrix \eqn{A} that satisfies
## \deqn{E[A] = \alpha(\Phi - 1),}
## where \eqn{\Phi} is the kinship matrix and \eqn{\alpha} is a nuisance scaling factor (determined by the unknown ancestral allele frequencies of each locus).
## Thus a \eqn{\Phi} estimate can be recovered from \eqn{A} after a separate step that estimates \eqn{\alpha = -\min E[A]}{M = -min E[A]} (see \code{\link{minAvgSubpops}} for one example).
##
## The matrix X (or the vectors returned by the function X) must have values only in c(0,1,2,NA), encoded to count the number of reference alleles at the locus, or NA for missing data.
##
## @param X Genotype matrix, BEDMatrix object, or a function X(mc) that returns the genotype matrix of all individuals at mc successive loci, and NULL when no loci are left.
## @param n Number of individuals (required only when X is a function, ignored otherwise)
## @param lociOnCols If true, X has loci on columns and individuals on rows; if false, loci are on rows and individuals on columns. Has no effect if X is a function.  If X is a BEDMatrix object, lociOnCols=TRUE is set automatically.
## @param memLim Memory limit in GB used to calculate the "chunk size" (numbers of SNPs). Note memory usage is somewhat underestimated and is not controlled strictly.  Default is 2GB, except in linux it is the free memory in the system times 0.7.
##
## @return The A matrix.
##
## @examples
## \dontrun{
## ## This example assumes input is in BED format and is loaded using BEDMatrix
## ## "file" is path to BED file (excluding .bed extension)
## library(BEDMatrix)
## X <- BEDMatrix(file) # load genotype matrix object
## A <- getA(X) # calculate A from genotypes
## }
getA <- function(X, n=NA, memLim=NA, lociOnCols=FALSE) {
    ## determine some behaviors depending on data type
    ## first validate class and set key booleans
    isFn <- FALSE
    if (class(X) == 'function') {
        isFn <- TRUE
        if (is.na(n)) stop('Fatal: missing number of individuals "n", which is required when X is a function.')
    } else if (class(X) == 'BEDMatrix') { # same as general matrix but transposed
        lociOnCols <- TRUE # this is always imposed for this particular format!
    } else if (class(X) != 'matrix') {
        stop('Fatal: X has unsupported class: ', class(X))
    } 
    
    ## extract dimensions from data (not possible for function version)
    if (!isFn) {
        if (lociOnCols) {
            if (!is.na(n) && n != nrow(X)) 
                warning('User set number of samples that does not match X dimensions (will go with latter): ', n, ' != ', nrow(X))
            n <- nrow(X)
            m <- ncol(X)
        } else {
            if (!is.na(n) && n != ncol(X)) 
                warning('User set number of samples that does not match X dimensions (will go with latter): ', n, ' != ', ncol(X))
            n <- ncol(X)
            m <- nrow(X)
        }
    } 
    
    ## initialize desired matrix
    A <- matrix(0, nrow=n, ncol=n)
    M <- matrix(0, nrow=n, ncol=n) # normalization now varies per individual pair (this tracks NAs, so subtract from overall m below)

    ## infer the number of SNPs to break data into, since we're limited by memory
    mc <- getMemLimM(m, n, memLim)

    ## navigate chunks
    mcis <- seq.int(1, m, mc)
    for (mci in mcis) {
        is <- mci:min(mci+mc-1,m) ## range of SNPs to extract in this chunk
        if (isFn) {
            Xi <- X( length(is) ) # get next SNPs
            if (is.null(Xi)) break # stop when SNPs run out (only happens for functions X, not matrices)
        } else if (lociOnCols) {
            Xi <- t(X[,is, drop=FALSE]) # transpose for our usual setup
        } else  {
            Xi <- X[is, , drop=FALSE]
        }

        ## before passing along to my RcppEigen code, I need to make sure the genotypes are treated by R as integers or RcppEigen dies on me
        ## I'm not sure if this always works though...
        if (storage.mode(Xi) != 'integer') storage.mode(Xi) <- 'integer'
        
        ## solve chunk using very efficient RcppEigen code!
        ## it is both runtime and memory efficient!
        obj <- getMAInt(Xi)
        ## increment each part
        A <- A + obj$SA # increment sum of A_{ijk}'s
        M <- M + obj$M # increment M's too
    }

    ## return final estimate!
    A/M - 1
}
