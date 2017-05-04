context('popkin_BEDMatrix')

## message('getwd: ', getwd())

if (suppressMessages(suppressWarnings(require(BEDMatrix)))) {

    ## loads Rdata matrices to test
    load('Xs.RData')
    rm(X) # pretend we don't have X!

    ## load the BED X instead
    X <- suppressMessages(suppressWarnings(BEDMatrix('Xs')))

    ## only repeat tests where genotypes X is input!
    
    test_that("function returns precomputed values: getA", {
        expect_equal(getA(X), A)
    })

    ## higher-level tests now!

    test_that("function returns precomputed values: popkin", {
        expect_equal(popkin(X), Phi0)
        expect_equal(popkin(X, subpops0), Phi0)
        expect_equal(popkin(X, subpops), Phi)
    })

}
