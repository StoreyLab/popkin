context('popkin_Rdata')

## message('getwd: ', getwd())

## ## introduce an error on purpose
## test_that("testing is run at all!", {
##     expect_equal(1, 2)
## })

## loads Rdata matrices to test
load('Xs.RData')

## start with lower-level/internal tests, more informative that higher-level function errors

test_that("validate_kinship works", {
    # validate positive examples
    expect_silent( validate_kinship( Phi ) )
    expect_silent( validate_kinship( Phi0 ) )
    expect_silent( validate_kinship( A ) ) # not real kinship but satisfies requirements

    # negative examples
    # dies if input is missing
    expect_error( validate_kinship() )
    
    # and if input is not a matrix
    expect_error( validate_kinship( 1:5 ) )
    
    # and for non-numeric matrices
    char_mat <- matrix(c('a', 'b', 'c', 'd'), nrow=2)
    expect_error( validate_kinship( char_mat ) )
    
    # and non-square matrices
    non_kinship <- matrix(1:2, nrow=2)
    expect_error( validate_kinship( non_kinship ) )
})


test_that("get_mem_lim returns positive numbers", {
    mem <- get_mem_lim()
    expect_equal(class(mem), 'numeric')
    expect_true(mem > 0)
})

test_that("get_mem_lim_m returns positive numbers", {
    mc <- get_mem_lim_m(n=1000, mem=2, m=100000) # chunk size, setting memory manually, set number of SNPs too
    expect_equal(length(mc), 1)
    expect_equal(class(mc), 'numeric')
    expect_true(mc > 0)
    mc <- get_mem_lim_m(n=1000, mem=2) # chunk size, setting memory manually, omit m
    expect_equal(length(mc), 1)
    expect_equal(class(mc), 'numeric')
    expect_true(mc > 0)
    mc <- get_mem_lim_m(n=1000) # chunk size, inferring free memory from system, omit m
    expect_equal(length(mc), 1)
    expect_equal(class(mc), 'numeric')
    expect_true(mc > 0)
})

test_that("function returns precomputed values: weights_subpops", {
    expect_equal(weights_subpops(subpops0), w0)
    expect_equal(weights_subpops(subpops), w)
    ## make sure dimensions match
    expect_equal(length(w0), nrow(Phi0))
    expect_equal(length(w), nrow(Phi))
    ## test the basic qualities of weights
    expect_equal(sum(w0), 1)
    expect_equal(sum(w), 1)
    expect_true(all(w0 > 0))
    expect_true(all(w > 0))
    expect_true(all(w0 < 1))
    expect_true(all(w < 1))
})

test_that("function returns precomputed values: getA", {
    expect_equal(getA(X), A)
    expect_equal(getA(X+0), A) # turns numeric
    expect_equal(getA(2L-X), A)
    expect_equal(getA(2-X), A) # numeric version again
    expect_equal(nrow(A), ncol(A))
})

test_that("function returns precomputed values: min_mean_subpops", {
    expect_equal(min_mean_subpops(A), min(A))
    expect_equal(min_mean_subpops(A), Amin0)
    expect_equal(min_mean_subpops(A, subpops0), Amin0)
    expect_equal(min_mean_subpops(A, subpops), Amin)
})

## higher-level tests now!

test_that("function returns precomputed values: popkin", {
    expect_equal(popkin(X), Phi0)
    expect_equal(popkin(X, subpops0), Phi0)
    expect_equal(popkin(X, subpops), Phi)
    expect_equal(popkin(X+0, subpops), Phi)
    expect_equal(popkin(2L-X, subpops), Phi)
    expect_equal(popkin(2-X, subpops), Phi)
})

test_that("popkin preserves names of individuals", {
    # in this case X is the regular tall matrix (individuals along columns)
    expect_equal(colnames(X), colnames(Phi))
    expect_equal(colnames(X), rownames(Phi))
})

test_that("function returns precomputed values: rescale_popkin", {
    expect_equal(rescale_popkin(Phi0, min_kinship = phiMin0), Phi)
    expect_equal(rescale_popkin(Phi0, subpops), Phi)
    expect_equal(rescale_popkin(Phi, subpops0), Phi0)
    expect_equal(rescale_popkin(Phi), Phi0)
})

test_that("function returns precomputed values: fst", {
    expect_equal(fst(Phi), fst)
    expect_equal(fst(Phi, w0), fst)
    expect_equal(fst(Phi, w), fstW)
    ## type of return value
    expect_equal(length(fstW), 1)
    expect_equal(length(fst), 1)
    ## Fst inequalities
    expect_true(fstW >= 0)
    expect_true(fst >= 0)
    expect_true(fstW <= 1)
    expect_true(fst <= 1)
})

test_that("function returns precomputed values: inbr", {
    expect_equal(inbr(Phi), inbr)
})

test_that("function returns precomputed values: pwfst", {
    expect_equal(pwfst(Phi), pwF)
    expect_equal(pwfst(Phi0), pwF)
    expect_equivalent(diag(pwF), rep.int(0, nrow(pwF))) # test that diagonal is zero ("equivalent" ignores label mismatches)
    expect_true(max(pwF) <= 1)
    ## note estimates may be slightly negative though
})

test_that("inbr_diag works", {
    # dies when kinship matrix is missing
    expect_error( inbr_diag() )
    # make sure precomputed values match
    expect_equal( inbr_diag(Phi), PhiInbr )
    # test list version (just duplicates things)
    expect_equal( inbr_diag(list(Phi, Phi)), list(PhiInbr, PhiInbr) )
})

test_that("neff works", {

    # a small toy case!
    F <- 0.4 # for this case to produce negative weights, need F > 1/3
    K <- (1+F)/2 # the self kinship values
    Phi3 <- matrix(c(K,F,0,F,K,F,0,F,K),nrow=3)
    w3 <- c(0.4, 0.2, 0.4) # dummy weights for a test

    # default case is sum of elements of inverse matrix!
    expect_equal(neff(Phi3, retW=FALSE, nonneg=FALSE), sum(solve(Phi3)))
    # non-max case returns inverse of mean kinship (with uniform weights)
    expect_equal(neff(Phi3, retW=FALSE, max=FALSE), 1/mean(Phi3))
    # non-max case returns inverse of weighted mean kinship (with provided non-uniform weights)
    expect_equal(neff(Phi3, retW=FALSE, max=FALSE, w=w3), 1/drop(w3 %*% Phi3 %*% w3))
    # test that gradient descent is the default
    expect_equal(neff(Phi3, retW=FALSE), neff(Phi3, retW=FALSE, algo='G'))
    # basic inequalities
    expect_true(neff(Phi3, retW=FALSE) >= 1) # min possible value
    expect_true(neff(Phi3, retW=FALSE, algo='G') >= 1) # min possible value
    expect_true(neff(Phi3, retW=FALSE, algo='N') >= 1) # min possible value
    expect_true(neff(Phi3, retW=FALSE, algo='H') >= 1) # min possible value
    expect_true(neff(Phi3, retW=FALSE, nonneg=FALSE) >= 1) # min possible value
    expect_true(neff(Phi3, retW=FALSE, max=FALSE) >= 1) # min possible value
    expect_true(neff(Phi3, retW=FALSE, max=FALSE, w=w3) >= 1) # min possible value
    expect_true(neff(Phi3, retW=FALSE) <= 2*nrow(Phi3)) # max possible value
    expect_true(neff(Phi3, retW=FALSE, algo='G') <= 2*nrow(Phi3)) # max possible value
    expect_true(neff(Phi3, retW=FALSE, algo='N') <= 2*nrow(Phi3)) # max possible value
    expect_true(neff(Phi3, retW=FALSE, algo='H') <= 2*nrow(Phi3)) # max possible value
    expect_true(neff(Phi3, retW=FALSE, nonneg=FALSE) <= 2*nrow(Phi3)) # max possible value
    expect_true(neff(Phi3, retW=FALSE, nonneg=FALSE) >= neff(Phi3, retW=FALSE)) # the max nEff should really be larger than a non-max case (non-optimal weights)
    expect_true(neff(Phi3, retW=FALSE, nonneg=FALSE) >= neff(Phi3, retW=FALSE, max=FALSE)) # the max nEff should really be larger than a non-max case (non-optimal weights)
    expect_true(neff(Phi3, retW=FALSE, nonneg=FALSE) >= neff(Phi3, retW=FALSE, max=FALSE, w=w3)) # the max nEff should really be larger than a non-max case (non-optimal weights)
    expect_true(neff(Phi3, retW=FALSE) >= neff(Phi3, retW=FALSE, max=FALSE)) # hope the numeric max with non-negative weights gives a larger value than a uniform weights estimate (actually not gauranteed to perform better)
    expect_true(neff(Phi3, retW=FALSE, algo='G') >= neff(Phi3, retW=FALSE, max=FALSE))
    expect_true(neff(Phi3, retW=FALSE, algo='N') >= neff(Phi3, retW=FALSE, max=FALSE))
    expect_true(neff(Phi3, retW=FALSE, algo='H') >= neff(Phi3, retW=FALSE, max=FALSE))

    # construct some other toy data tests
    # check that we get 2*n on unstructured kinship matrix
    expect_equal(neff(diag(1/2, 10, 10), retW=FALSE), 20)
    # check that we get K on extreme Fst=1 independent subpopulations
    # NOTE: FAILS because it's not invertible... what should we do here???
    #    expect_equal(neff(matrix(c(1, 1, 0, 1, 1, 0, 0, 0, 1), nrow=3), retW=FALSE), 2)
    # maybe construct Fst<1 case and compare to theoretical expectation there

    # normally weights are returned too, but their values are less constrained (by construction they sum to one, that's it!)
    # do test both max versions since nEff is not directly constructed from the weights (test that it's what it should be)
    # test outputs in that setting now
    obj <- neff(Phi3, algo='G') # this tests Gradient version
    expect_equal(class(obj), 'list') # return class is list
    expect_equal(length(obj), 2) # only have two elements
    # roughly retest that first element is an nEff
    expect_true(obj$neff >= 1) # min possible value
    expect_true(obj$neff <= 2*nrow(Phi3)) # max possible value
    # roughly test weights
    expect_equal(sum(obj$w), 1) # verify that weights sum to 1
    expect_equal(obj$neff, 1/meanKin(Phi3, obj$w)) # verify that neff has the value it should have given the weights
    
    obj <- neff(Phi3, algo='N') # this tests Newton version
    expect_equal(class(obj), 'list') # return class is list
    expect_equal(length(obj), 2) # only have two elements
    # roughly retest that first element is an nEff
    expect_true(obj$neff >= 1) # min possible value
    expect_true(obj$neff <= 2*nrow(Phi3)) # max possible value
    # roughly test weights
    expect_equal(sum(obj$w), 1) # verify that weights sum to 1
    expect_equal(obj$neff, 1/meanKin(Phi3, obj$w)) # verify that neff has the value it should have given the weights
    
    obj <- neff(Phi3, algo='H') # this tests Heuristic version
    expect_equal(class(obj), 'list') # return class is list
    expect_equal(length(obj), 2) # only have two elements
    # roughly retest that first element is an nEff
    expect_true(obj$neff >= 1) # min possible value
    expect_true(obj$neff <= 2*nrow(Phi3)) # max possible value
    # roughly test weights
    expect_equal(sum(obj$w), 1) # verify that weights sum to 1
    expect_equal(obj$neff, 1/meanKin(Phi3, obj$w)) # verify that neff has the value it should have given the weights
    
    obj <- neff(Phi3, nonneg=FALSE) # this tests optimal version (with possibly negative weights)
    expect_equal(class(obj), 'list') # return class is list
    expect_equal(length(obj), 2) # only have two elements
    # roughly retest that first element is an nEff
    expect_true(obj$neff >= 1) # min possible value
    expect_true(obj$neff <= 2*nrow(Phi3)) # max possible value
    # roughly test weights
    expect_equal(sum(obj$w), 1) # verify that weights sum to 1
    expect_equal(obj$neff, 1/meanKin(Phi3, obj$w)) # verify that neff has the value it should have given the weights
    expect_true(min(obj$w) < 0) # this example must have negative weights, or it's useless!
})

