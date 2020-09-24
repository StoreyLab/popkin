context('popkin_BEDMatrix')

# message('getwd: ', getwd())

if (suppressMessages(suppressWarnings(require(BEDMatrix)))) {

    # loads Rdata matrices to test
    load('Xs.RData')
    rm(X) # pretend we don't have X!

    # load the BED X instead
    X <- suppressMessages(suppressWarnings(BEDMatrix('Xs')))

    # only repeat tests where genotypes X is input!
    
    test_that("function returns precomputed values: get_A", {
        expect_silent( obj <- get_A( X ) )
        expect_equal( obj$A, A )
        expect_equal( obj$M, M )
    })

    # higher-level tests now!

    test_that("function returns precomputed values: popkin", {
        expect_equal(popkin(X), Phi0)
        expect_equal(popkin(X, subpops0), Phi0)
        expect_equal(popkin(X, subpops), Phi)
    })

    test_that("popkin preserves names of individuals", {
        # in this case X is the wide matrix (individuals along rows)
        expect_equal(rownames(X), colnames(Phi))
        expect_equal(rownames(X), rownames(Phi))
    })

}
