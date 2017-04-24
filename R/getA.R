#' Compute A matrix from genotypes
#'
#' Given the biallelic genotypes of \eqn{n} individuals, this function returns the \eqn{n}-by-\eqn{n} matrix \eqn{A} that satisfies
#' \deqn{E[A] = M(\Phi - 1),}
#' where \eqn{\Phi} is the kinship matrix and \eqn{M} is a nuisance scaling factor (determined by the unknown ancestral allele frequencies of each locus).
#' Thus a \eqn{\Phi} estimate can be recovered from \eqn{A} after a separate step that estimates \eqn{M = -\min E[A]}{M = -min E[A]} (see \code{\link{getAEminSubpops}} for one example).
#'
#' The matrix X (or the vectors returned by the function X) must have values only in c(0,1,2,NA), encoded to count the number of reference alleles at the locus, or NA for missing data.
#'
#' @param X Genotype matrix, BEDMatrix object, or a function that returns the genotypes of all individuals at successive loci each time it is called, and NULL when no loci are left.
#' @param n Number of individuals (required only when X is a function, ignored otherwise)
#' @param m Number of loci (optional, may truncate input when X is a function; ignored when X is a matrix or BEDMatrix object)
#' @param lociOnCols If true, X has loci on columns and individuals on rows; if false, loci are on rows and individuals on columns. Has no effect if X is a function.  If X is a BEDMatrix object, lociOnCols=TRUE is set automatically.
#' @param lowMem If true, code that runs through each SNP is used, which uses very low memory but is slower. If false, code that uses matrix algebra is used, which uses much more memory but is also faster.  Set to TRUE automaticaly when X is a function.
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
getA <- function(X, n=NA, m=NA, lociOnCols=FALSE, lowMem=FALSE) {
    ## SUPER low-mem version that processes input as it is read
    ## X(n) is a genetic SNP-reading function that also dies if row isn't length n
    ### NOTE: low-memory version of getA (handles missing data for large X but is slower than getA)
    ## since this is a big problem, let's go full low-mem!

    ## determine some behaviors depending on data type
    ## first validate class and set key booleans
    ## also immediately call high-memory version if requested
    classX <- class(X)
    isFn <- FALSE
    if (classX == 'function') {
        isFn <- TRUE
        if (is.na(n)) stop('Fatal: missing number of individuals "n", which is required when X is a function.')
    } else if (classX == 'BEDMatrix') { # same as general matrix but transposed
        ## this strange notation X[,] forces the BEDMatrix object to load into a regular matrix (fully in memory) that can then be transposed
        ## NOTE: calling directly t(X) on a BEDMatrix object X does nothing (no errors either, weird), and such an X dies with colMeans and potentially other matrix operations (didn't try crossprod())
        if (!lowMem) return( getA_hiMem(t(X[,])) )
        lociOnCols <- TRUE # this is always imposed for this particular format!
    } else if (classX == 'matrix') {
        if (!lowMem) return( getA_hiMem(X) )
    } else stop('Fatal: X has unsupported class: ', classX)
    
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
    i <- 1 # SNP counter (always advances)
    m2 <- 0 # actual number of processed SNPs (which were not fixed, so may be smaller than m and final i-1)
    a <- 0 # scalar adjustments to A...
    ## constants
    n2 <- n*n # handy shortcut
    
    ## navigate SNPs!
    while(TRUE) { # start an infinite loop (don't usually know when SNPs will run out)
        if (!is.na(m) && i > m) break # stop if we've reached SNP limit (only when m is known, either from X or set manually)

        if (isFn) {
            xi <- X() # get next SNP
            if (is.null(xi)) break # stop when SNPs run out (only happens for functions X, not matrices)
        } else if (lociOnCols)  {
            xi <- X[,i] # get i^th COLUMN
        } else  {
            xi <- X[i,] # get i^th row (assumes our usual setup)
        }
        i <- i+1 # increment SNP index (must do here for matrix/BEDMatrix cases, not at end, so we don't get stuck in an infinite loop when a fixed SNP is skipped or something like it)

        ## skip fixed SNPs
        pi <- mean(xi, na.rm=TRUE)/2 # standard mean times half
        if (0 >= pi || pi >= 1) next # always skip fixed SNPs (contribute zero to A)

        ## center, handle missingness as needed
        xi <- xi - 1 # center before missingness hack
        isn <- is.na(xi) # TRUE/FALSE status of NA positions
        if (any(isn)) { # missingness happens rarely, so the following can be skipped often
            M[,isn] <- M[,isn] + 1 # this will reduce counts in each column containing these individuals
            M[isn,!isn] <- M[isn,!isn] + 1 # ditto rows (but don't repeat columns we already did, so we don't double count)
            xi[isn] <- 0 # before applying cross product, to prevent NA errors, just set those values to zero and it works out!
        }

        ## update A running sums (main part of the calculation/estimator)
        ## goal here is to identify blocks of cases, figure out which strategy minimizes number of edits
        ## count homozygotes (heterozygotes xi == 0 contribute to zero in our formula)
        isp <- xi == 1
        ism <- xi == -1
        ## count the cases
        np <- sum(isp)
        nm <- sum(ism)
        ## now count the size of the matrix edits
        nnp <- np*np + nm*nm # number of 1's
        nnz <- n2 - (np+nm)^2
        ## this test decides which kind of update will be faster (which will require the smallest number of matrix edits)
        if (nnz >= nnp) {
            ## updates assuming the zero case is the most common case...
            A[isp,isp] <- A[isp,isp] + 1 # same blocks are +1
            A[ism,ism] <- A[ism,ism] + 1
            A[isp,ism] <- A[isp,ism] - 1 # diff blocks are -1
            A[ism,isp] <- A[ism,isp] - 1
        } else {
            ## in this case 1 is most common value, keep track of it in a separate scalar (faster updates)
            a <- a+1
            ## update values that differ...
            A[isp,ism] <- A[isp,ism] - 2 # diff blocks are -1 (adjusting for the overall +1 makes this a -2)
            A[ism,isp] <- A[ism,isp] - 2
            isz <- !isp & !ism # xi == 0 # need these indexes only in this case
            A[,isz] <- A[,isz] - 1 # columns, old zero blocks are -1 here (this and next)
            A[isz,!isz] <- A[isz,!isz] - 1 # rows, but not repeating columns we already did
        }
        m2 <- m2+1 # increment number of non-fixed SNPs
    }

    ## normalize estimate and return!
    A <- (A+a)/(m2-M) - 1 # turn into average, subtract -1 now!
}

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
#' @param lowMem If true, code that runs through each SNP is used, which uses very low memory but is slower. If false, code that uses matrix algebra is used, which uses much more memory but is also faster.  Set to TRUE automaticaly when X is a function.
#' @param mc For lowMem=TRUE, number of SNPs of each "chunk" to use high memory code on.
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
getA2 <- function(X, n=NA, m=NA, mc=100000, lociOnCols=FALSE, lowMem=FALSE) {
    ## TODO: specify chunks indirectly as a memory constraint (need to map dimensions to memory usage!)

    ## determine some behaviors depending on data type
    ## first validate class and set key booleans
    ## also immediately call high-memory version if requested
    classX <- class(X)
    isFn <- FALSE
    if (classX == 'function') {
        isFn <- TRUE
        if (is.na(n)) stop('Fatal: missing number of individuals "n", which is required when X is a function.')
    } else if (classX == 'BEDMatrix') { # same as general matrix but transposed
        ## this strange notation X[,] forces the BEDMatrix object to load into a regular matrix (fully in memory) that can then be transposed
        ## NOTE: calling directly t(X) on a BEDMatrix object X does nothing (no errors either, weird), and such an X dies with colMeans and potentially other matrix operations (didn't try crossprod())
        if (!lowMem) return( getA_hiMem(t(X[,])) )
        lociOnCols <- TRUE # this is always imposed for this particular format!
    } else if (classX == 'matrix') {
        if (!lowMem) return( getA_hiMem(X) )
    } else stop('Fatal: X has unsupported class: ', classX)
    
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

    ## navigate chunks
    mcis <- seq.int(1, m, mc)
    ## for (mci=1; mci < m; mci = mci+mc) {
    for (mci in mcis) {
        is <- mci:min(mci+mc-1,m) ## range of SNPs to extract in this chunk
        if (isFn) {
            Xi <- X( length(is) ) # get next SNPs
            if (is.null(Xi)) break # stop when SNPs run out (only happens for functions X, not matrices)
        } else if (lociOnCols) {
            Xi <- t(X[,is]) # transpose for our usual setup
        } else  {
            Xi <- X[is,]
        }

        ## solve chunk using high-mem code (turned into low mem only because we break problem up into chunks)
        obj <- getA_hiMem(Xi, returnParts=TRUE)
        ## increment each part
        A <- A + obj$SA # increment sum of A_{ijk}'s
        M <- M + obj$M # increment M's too
    }

    ## return final estimate!
    A/M - 1
}

## internal high-memory function that is called by main getA() as needed
getA_hiMem <- function(X, returnParts=FALSE) { # new version that handles missingness correctly
    ## less normalizations, no weights or ancestral AFs needed! this is the coolest and simplest of my estimators so far!
    ## this version handles missingness but stays efficient when there is no missingness
    ## also tests for and removes fixed SNPs
    
    ## handles missing values by setting relevant bits to zero this way:
    ## https://stackoverflow.com/questions/16535084/matrix-multiplication-with-scattered-na-values

    ## estimating total memory usage in bytes...
    ## GB = B*1024*1024*1024 
    ## mem(X) = m*n*4+40 # genotypes are stored as ints!
    ## mem(FHat) = m*8+40 # prediction, didn't test...
    ## fixed cases double X at least temporarily?
    ## X <- X-1 may also double X (but later, not at the same time)
    ## mem(M) = n*n*4+40 # assuming these are encoded as ints too
    ## mem(is.na(X)) ? I noticed saving this was somehow worse, I don't know if I should count it or now
    ## mem(SA) = n*n*4+40 # also assuming int, not sure...
    ## total mem
    ## = 2*(m*n*4+40) + m*8+40 + 2*(n*n*4+40)
    ## = 2*m*n*4 + m*8 + 2*n*n*4 + 5*40

    ## Part 1: remove fixed SNPs
    ## only count polymorphic SNPs
    FHat <- rowMeans(X, na.rm=TRUE)/2 # weights don't matter here, we only want to quickly identify fixed vs non-fixed
    if (any(FHat == 0) || any(FHat == 1)) { # test for problem cases
        X <- X[ 0<FHat & FHat<1 ,] # clean up X if needed!!!
    }

    ## Part 2, compute matrix of interest
    X <- X - 1 # center this way
    if (any(is.na(X))) {
        M <- crossprod( !is.na(X) ) # normalization factor is now a matrix, varies per pair, this should compute correct values
        ## NOTE: if X is missing "pc" proportion at random, then the mean expected value of M is = nrow(X) * (1-pc)^2
        X[is.na(X)] <- 0 # before applying cross product, to prevent NA errors, just set those values to zero and it works out!
    } else {
        M <- nrow(X) # when nothing is missing, the normalization is a scalar and it's simply the number of SNPs (rows of X)
    }
    SA <- crossprod(X) # compute sum_i A_{ijk} terms ("Sum A")
    if (returnParts) {
        return( list(SA=SA, M=M) ) # return these two matrices separately
    } else {
        return( SA/M - 1 ) # return final A estimate of interest!!!
    }
}

