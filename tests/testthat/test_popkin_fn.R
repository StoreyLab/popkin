context('popkin_fn')

## loads Rdata matrices to test
load('Xs.RData')

## construct artificial function that reads m loci at the time from matrix
## input is so tiny only one chunk will be read probably, but let's have a nice general solution just in case...
ml <- 1 # keeps track of last locus not yet read from (will be edited as global)
Xf <- function(m) {
    ## set up range of loci to access
    if (ml > nrow(X)) return(NULL) # return NULL if we've read everything (tells popkin to stop)
    is <- ml:min( (ml+m-1), nrow(X)) # this gives correct max length of "m"
    ## if (max(is) > nrow(X)) { # reduce range if it's out of bounds of matrix
    ##     is <- is[is <= nrow(X)] # reduce here
    ##     if (length(is) == 0) return(NULL) # if nothing was left, 
    ## }
    ml <<- max(is)+1 # update global with next value to read...
    X[is,] # return this subset!
}
n <- ncol(X) # need to pass it separately for function inputs

# NOTE: function version doesn't preserve names, so let's erase names from the matrices to compare to make sure they agree in that sense
# this is the proper way to do this (setting each of colnames/rownames to NULL is not enough!)
dimnames(A) <- NULL
dimnames(Phi0) <- NULL
dimnames(Phi) <- NULL

## only repeat tests where genotypes X is input!

test_that("function returns precomputed values: get_A", {
    expect_equal(get_A(Xf, n=n), A)
    ml <<- 1 # reset function after we're done!
})

## higher-level tests now!

test_that("function returns precomputed values: popkin", {
    expect_equal(popkin(Xf, n=n), Phi0)
    ml <<- 1 # reset function after we're done!
    expect_equal(popkin(Xf, subpops0), Phi0) # NOTE: n==length(subpops0) is inferred
    ml <<- 1 # reset function after we're done!
    expect_equal(popkin(Xf, subpops), Phi)   # NOTE: n==length(subpops) is inferred
    ml <<- 1 # reset function after we're done!
})

test_that("popkin individual names do not exist in this case", {
    # looks dumb to be testing the precomputed value, but earlier test proved that this agrees with fresh outputs so meh
    expect_null(colnames(Phi))
    expect_null(rownames(Phi))
})
