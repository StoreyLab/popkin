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
    expect_equal(class(mc), 'numeric')
    expect_true(mc > 0)
    mc <- getMemLimM(n=1000, mem=2) # chunk size, setting memory manually, omit m
    expect_equal(class(mc), 'numeric')
    expect_true(mc > 0)
    mc <- getMemLimM(n=1000) # chunk size, inferring free memory from system, omit m
    expect_equal(class(mc), 'numeric')
    expect_true(mc > 0)
})

test_that("function returns precomputed values: getA", {
    expect_equal(getA(X), A)
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
})

## higher-level tests now!

test_that("function returns precomputed values: popkin", {
    expect_equal(popkin(X), Phi0)
    expect_equal(popkin(X, subpops0), Phi0)
    expect_equal(popkin(X, subpops), Phi)
})

test_that("function returns precomputed values: rescalePopkin", {
    expect_equal(rescalePopkin(Phi0, phiMin=phiMin0), Phi)
    expect_equal(rescalePopkin(Phi0, subpops), Phi)
})

test_that("function returns precomputed values: fst", {
    expect_equal(fst(Phi), fst)
})

test_that("function returns precomputed values: inbr", {
    expect_equal(inbr(Phi), inbr)
})

test_that("function returns precomputed values: inbrDiag", {
    expect_equal(inbrDiag(Phi), PhiInbr)
})

