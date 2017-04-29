#' Compute A matrix from genotypes
#'
#' Given the biallelic genotypes of \eqn{n} individuals, this function returns the \eqn{n}-by-\eqn{n} matrix \eqn{A} that satisfies
#' \deqn{E[A] = M(\Phi - 1),}
#' where \eqn{\Phi} is the kinship matrix and \eqn{M} is a nuisance scaling factor (determined by the unknown ancestral allele frequencies of each locus).
#' Thus a \eqn{\Phi} estimate can be recovered from \eqn{A} after a separate step that estimates \eqn{M = -\min E[A]}{M = -min E[A]} (see \code{\link{getAEminSubpops}} for one example).
#'
#' The matrix X (or the vectors returned by the function X) must have values only in c(0,1,2,NA), encoded to count the number of reference alleles at the locus, or NA for missing data.
#'
#' @param X Genotype matrix, BEDMatrix object, or a function X(mc) that returns the genotype matrix of all individuals at mc successive loci, and NULL when no loci are left.
#' @param n Number of individuals (required only when X is a function, ignored otherwise)
#' @param m Number of loci (optional, may truncate input when X is a function; ignored when X is a matrix or BEDMatrix object)
#' @param lociOnCols If true, X has loci on columns and individuals on rows; if false, loci are on rows and individuals on columns. Has no effect if X is a function.  If X is a BEDMatrix object, lociOnCols=TRUE is set automatically.
#' @param memLim Memory limit in GB used to calculate the "chunk size" (numbers of SNPs). NOTE: memory will not be controlled strictly, just approximately.  Default is 2GB, except in linux it is the free memory in the system times 0.7.
#'
#' @return The A matrix.
#'
#' @examples
#' \dontrun{
#' ## This example assumes input is in BED format and is loaded using BEDMatrix
#' ## "file" is path to BED file (excluding .bed extension)
#' library(BEDMatrix)
#' X <- BEDMatrix(file) # load genotype matrix object
#' A <- getA(X) # calculate A from genotypes
#' }
#' 
#' @export
getA <- function(X, n=NA, m=NA, memLim=NA, lociOnCols=FALSE) {
    ## determine some behaviors depending on data type
    ## first validate class and set key booleans
    ## also immediately call high-memory version if requested
    classX <- class(X)
    isFn <- FALSE
    if (classX == 'function') {
        isFn <- TRUE
        if (is.na(n)) stop('Fatal: missing number of individuals "n", which is required when X is a function.')
    } else if (classX == 'BEDMatrix') { # same as general matrix but transposed
        lociOnCols <- TRUE # this is always imposed for this particular format!
    } else if (classX != 'matrix') {
        stop('Fatal: X has unsupported class: ', classX)
    } 
    
    ## extract dimensions from data (not possible for function version)
    if (!isFn) {
        if (lociOnCols) {
            if (!is.na(n) && n != nrow(X)) 
                warning('User set number of samples that does not match X dimensions (will go with latter): ', n, ' != ', nrow(X))
            if (!is.na(m) && m != ncol(X)) 
                warning('User set number of loci that does not match X dimensions (will go with latter): ', m, ' != ', ncol(X))
            n <- nrow(X)
            m <- ncol(X)
        } else {
            if (!is.na(n) && n != ncol(X)) 
                warning('User set number of samples that does not match X dimensions (will go with latter): ', n, ' != ', ncol(X))
            if (!is.na(m) && m != nrow(X)) 
                warning('User set number of loci that does not match X dimensions (will go with latter): ', m, ' != ', nrow(X))
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
        message('chunk ', match(mci,mcis), '/', length(mcis)) # DEBUGGING
        is <- mci:min(mci+mc-1,m) ## range of SNPs to extract in this chunk
        if (isFn) {
            Xi <- X( length(is) ) # get next SNPs
            if (is.null(Xi)) break # stop when SNPs run out (only happens for functions X, not matrices)
        } else if (lociOnCols) {
            Xi <- t(X[,is]) # transpose for our usual setup
        } else  {
            Xi <- X[is,]
        }

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

#' @useDynLib popkin
#' @importFrom Rcpp sourceCpp
NULL


## http://dirk.eddelbuettel.com/code/rcpp/Rcpp-package.pdf
## http://dirk.eddelbuettel.com/code/rcpp/Rcpp-FAQ.pdf
## https://cran.r-project.org/web/packages/RcppEigen/vignettes/RcppEigen-Introduction.pdf
## https://stackoverflow.com/questions/13956292/rcpp-inline-creating-and-calling-additional-functions
## https://stackoverflow.com/questions/27490659/rcppeigen-going-from-inline-to-a-cpp-function-in-a-package-and-map

## library(RcppEigen)
## library(inline)

## crossprodSelfIntCpp <- '
## using Eigen::Map;
## using Eigen::MatrixXi;
## using Eigen::Lower;
## // Map the integer matrix AA from R
## const Map<MatrixXi> A(as<Map<MatrixXi>  >(AA));
## // evaluate and return crossprod(A)
## const int           n(A.cols());
## MatrixXi          AtA(MatrixXi(n, n).setZero().selfadjointView<Lower>().rankUpdate(A.adjoint()));
## // wrap and return!
## return wrap(AtA);
## '
## crossprodSelfInt <- cxxfunction(signature(AA = "matrix"), crossprodSelfIntCpp, "RcppEigen")
## SA2 <- crossprodSelfInt(X) # verified this returns ints!

getMemLimM <- function(m, n, M) {
    ## NOTE: in estimating the chunk size, we don't know ahead of time if there are NAs or not!
    ## so these calculations assume M in getA_himem is an n*n matrix
    ## if there aren't any NAs, M is a scalar and we end up underestimating memory usage (better than the other way around)
    
    ## try to get total memory from the system if M wasn't specified, so it works reasonably
    if (is.na(M)) {
        M <- getMemLim() # infer a reasonable default from system!
    } else {
        ## assuming M is in GB, let's convert to bytes
        M <- M*1024*1024*1024
    }

    ## OLD...
    ## estimating total memory usage in bytes, for getA_hiMem
    ## mem(X) = m*n*4+40 # genotypes are stored as ints!
    ## mem(FHat) = m*8+40 # double vector
    ## fixed cases double X at least temporarily?
    ## X <- X-1 may also double X (but later, not at the same time)
    ## mem(M) = n*n*4+40 # as int (RcppEigen), double otherwise
    ## mem(is.na(X)) same as X
    ## mem(!is.na(X)) same as X
    ## mem(SA) = n*n*4+40 # as int (RcppEigen), double otherwise
    
    ## total mem = M # version for int M,SA (RcppEigen), no more FHat
    ## = 2*(m*n*4+40) + 2*(n*n*4+40)
    ## = 2*m*n*4 + 2*n*n*4 + 4*40
    ## = 8*m*n + 8*n*n + 160
    ## = 8*(m*n + n*n + 20)
    ## given fixed n, solve for m:
    ## M = 8*(m*n + n*n + 20)
    ## m = (M/8 - 20 - n*n)/n

    ## apply formula, get maximum m (number of SNPs) given n and the memory requested
    mc <- (M/8 - 20 - n*n)/n
    
    if (m < mc) {
        mc <- m # use the smaller one
    } else {
        ## should "redistribute" based on number of chunks, to lower memory even more per iteration
        mc <- floor( m/ceiling(m/mc) ) # this lowers mc even more, balances load better
    }
    
    Mact <- 8*(mc*n + n*n + 20)
    message('Choice of mc should limit mem to ', round( Mact/(1024*1024*1024), 2 ), ' GB')

    return(mc) # return desired value
}

