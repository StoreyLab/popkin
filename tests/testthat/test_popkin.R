context('popkin_Rdata')

## message('getwd: ', getwd())

## ## introduce an error on purpose
## test_that("testing is run at all!", {
##     expect_equal(1, 2)
## })

## loads Rdata matrices to test
load('Xs.RData')

## start with lower-level/internal tests, more informative that higher-level function errors

test_that("getMemLim returns positive numbers", {
    mem <- getMemLim()
    expect_equal(class(mem), 'numeric')
    expect_true(mem > 0)
})

test_that("getMemLimM returns positive numbers", {
    mc <- getMemLimM(n=1000, mem=2, m=100000) # chunk size, setting memory manually, set number of SNPs too
    expect_equal(length(mc), 1)
    expect_equal(class(mc), 'numeric')
    expect_true(mc > 0)
    mc <- getMemLimM(n=1000, mem=2) # chunk size, setting memory manually, omit m
    expect_equal(length(mc), 1)
    expect_equal(class(mc), 'numeric')
    expect_true(mc > 0)
    mc <- getMemLimM(n=1000) # chunk size, inferring free memory from system, omit m
    expect_equal(length(mc), 1)
    expect_equal(class(mc), 'numeric')
    expect_true(mc > 0)
})

test_that("function returns precomputed values: weightsSubpops", {
    expect_equal(weightsSubpops(subpops0), w0)
    expect_equal(weightsSubpops(subpops), w)
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

test_that("function returns precomputed values: minAvgSubpops", {
    expect_equal(minAvgSubpops(A), min(A))
    expect_equal(minAvgSubpops(A), Amin0)
    expect_equal(minAvgSubpops(A, subpops0), Amin0)
    expect_equal(minAvgSubpops(A, subpops), Amin)
})

test_that("function returns precomputed values: getKinshipFromA", {
    expect_equal(getKinshipFromA(A, Amin0), Phi0)
    expect_equal(getKinshipFromA(A, Amin), Phi)
    expect_equal(nrow(Phi0), ncol(Phi0))
    expect_equal(nrow(Phi), ncol(Phi))
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

test_that("function returns precomputed values: rescalePopkin", {
    expect_equal(rescalePopkin(Phi0, phiMin=phiMin0), Phi)
    expect_equal(rescalePopkin(Phi0, subpops), Phi)
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

test_that("function returns precomputed values: inbrDiag", {
    expect_equal(inbrDiag(Phi), PhiInbr)
})

