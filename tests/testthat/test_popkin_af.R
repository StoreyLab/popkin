#################
### popkin_af ###
#################

# for the following tests, simulate a matrix with subpopulation allele frequencies
n <- 10
m <- 100 # P's is smaller than X's m
# group subpopulations in pairs, just for this test
subpops <- ceiling( (1 : n) / 2 )
# unstructured allele frequencies, meh
P <- matrix(
    runif( n * m ),
    nrow = m,
    ncol = n
)
# add labels for tests
colnames( P ) <- paste0( 'S', 1 : n )
rownames( P ) <- paste0( 'L', 1 : m )
# true M when there's no missingness
M_true <- matrix( m, n, n )
colnames( M_true ) <- colnames( P )
rownames( M_true ) <- colnames( P )

# make a copy of P with random missingness
miss <- 0.1 # missingness rate
P_miss <- P
P_miss[ sample( n * m, n * m * miss ) ] <- NA
# NOTE: if any rows are entirely fixed or missing, popkin_af should still work!
# this is true M for missingness case
M_miss_true <- crossprod( !is.na( P_miss ) )

test_that( "true M are ok for popkin_af", {
    # M matrix didn't have equivalent in past runs, test in more detail here
    expect_true( is.matrix( M_true ) )
    expect_equal( nrow( M_true ), n )
    expect_equal( ncol( M_true ), n )
    expect_true( all( M_true == m ) )
    
    # repeat with missingness
    expect_true( is.matrix( M_miss_true ) )
    expect_equal( nrow( M_miss_true ), n )
    expect_equal( ncol( M_miss_true ), n )
    # with missingness we just have bounds
    expect_true( all( M_miss_true <= m ) )
    expect_true( all( M_miss_true >= 0 ) )
})

test_that( "popkin_af_basic works", {
    # expect error if P is missing
    expect_error( popkin_af_basic( ) )
    # cause error if subpopulations and P don't match in dimensions
    expect_error( popkin_af_basic( P, subpops = subpops[ -1 ] ) )

    # now proper run
    # NOTE: saving `coancestry` as global so next block can compare against it
    expect_silent(
        coancestry <<- popkin_af_basic( P )
    )
    # check dimensions and data ranges
    expect_true( is.numeric( coancestry ) )
    expect_true( is.matrix( coancestry ) )
    expect_equal( nrow( coancestry ), n )
    expect_equal( ncol( coancestry ), n )
    expect_true( !anyNA( coancestry ) )
    # because of the way Amin is picked here (subpops = NULL), the minimum is exactly zero
    expect_true( all( coancestry >= 0 ) )
    expect_true( all( coancestry <= 1 ) )
    # test label transfer
    expect_equal( colnames( coancestry ), colnames( P ) )
    expect_equal( rownames( coancestry ), colnames( P ) )

    # test want_M version
    expect_silent(
        obj <- popkin_af_basic( P, want_M = TRUE )
    )
    # require equality here with previous run
    expect_equal( obj$coancestry, coancestry )
    # this is not a matrix in this case, but just to prepare for version with missingness later
    expect_equal( obj$M, m )
    
    # run transposed version
    coancestry2 <- popkin_af_basic( t( P ), loci_on_cols = TRUE )
    # require equality here with previous run
    expect_equal( coancestry2, coancestry )
    
    # run with subpops
    expect_silent(
        coancestry <- popkin_af_basic( P, subpops = subpops )
    )
    # check dimensions and data ranges
    expect_true( is.numeric( coancestry ) )
    expect_true( is.matrix( coancestry ) )
    expect_equal( nrow( coancestry ), n )
    expect_equal( ncol( coancestry ), n )
    expect_true( !anyNA( coancestry ) )
    # because of the way Amin is picked here (paired subpops), the minimum is slightly below zero
    expect_true( all( coancestry >= -1 ) )
    expect_true( all( coancestry <= 1 ) )
    # test label transfer
    expect_equal( colnames( coancestry ), colnames( P ) )
    expect_equal( rownames( coancestry ), colnames( P ) )
})

test_that( "popkin_af_basic_na works", {
    # expect error if P is missing
    expect_error( popkin_af_basic_na( ) )
    # cause error if subpopulations and P don't match in dimensions
    expect_error( popkin_af_basic_na( P, subpops = subpops[ -1 ] ) )

    # now proper run
    # here we use P without missingness, expect to match previous version
    expect_silent(
        coancestry2 <- popkin_af_basic_na( P )
    )
    expect_equal( coancestry2, coancestry )

    # and want_M version
    expect_silent(
        obj <- popkin_af_basic_na( P, want_M = TRUE )
    )
    expect_equal( obj$coancestry, coancestry )
    expect_equal( obj$M, M_true )

    # now run with P with missingness!
    # export this one too
    expect_silent(
        coancestry_miss <<- popkin_af_basic_na( P_miss )
    )
    # check dimensions and data ranges
    expect_true( is.numeric( coancestry_miss ) )
    expect_true( is.matrix( coancestry_miss ) )
    expect_equal( nrow( coancestry_miss ), n )
    expect_equal( ncol( coancestry_miss ), n )
    expect_true( !anyNA( coancestry_miss ) )
    # because of the way Amin is picked here (subpops = NULL), the minimum is exactly zero
    expect_true( all( coancestry_miss >= 0 ) )
    expect_true( all( coancestry_miss <= 1 ) )
    # test label transfer
    expect_equal( colnames( coancestry_miss ), colnames( P ) )
    expect_equal( rownames( coancestry_miss ), colnames( P ) )

    # want_M version with missigness
    expect_silent(
        obj_miss <- popkin_af_basic_na( P_miss, want_M = TRUE )
    )
    # compare to proper calculation from data with missingness
    expect_equal( obj_miss$coancestry, coancestry_miss )
    expect_equal( obj_miss$M, M_miss_true )
    
    # and transposed version (stick with want_M to test more things)
    expect_silent(
        obj_miss2 <- popkin_af_basic_na( t( P_miss ), want_M = TRUE, loci_on_cols = TRUE )
    )
    # test whole object at once
    expect_equal( obj_miss2, obj_miss )

    # run with subpops
    expect_silent(
        obj_miss_subpops <- popkin_af_basic_na( P_miss, want_M = TRUE, subpops = subpops )
    )
    # extract elements
    coancestry_miss_subpops <- obj_miss_subpops$coancestry
    # check dimensions and data ranges
    expect_true( is.numeric( coancestry_miss_subpops ) )
    expect_true( is.matrix( coancestry_miss_subpops ) )
    expect_equal( nrow( coancestry_miss_subpops ), n )
    expect_equal( ncol( coancestry_miss_subpops ), n )
    expect_true( !anyNA( coancestry_miss_subpops ) )
    # because of the way Amin is picked here (paired subpops), the minimum is slightly below zero
    expect_true( all( coancestry_miss_subpops >= -1 ) )
    expect_true( all( coancestry_miss_subpops <= 1 ) )
    # test label transfer
    expect_equal( colnames( coancestry_miss_subpops ), colnames( P ) )
    expect_equal( rownames( coancestry_miss_subpops ), colnames( P ) )
    # M is the same here whether subpops are used or not
    expect_equal( obj_miss_subpops$M, obj_miss$M )
})

# this tests final, exported version (the other ones are for checks only)
test_that( "popkin_af works", {
    # expect error if P is missing
    expect_error( popkin_af( ) )
    # cause error if subpopulations and P don't match in dimensions
    expect_error( popkin_af( P, subpops = subpops[ -1 ] ) )

    # now proper run
    # here we use P without missingness, expect to match previous version
    expect_silent(
        coancestry2 <- popkin_af( P )
    )
    expect_equal( coancestry2, coancestry )
    # repeat with want_M
    expect_silent(
        obj <- popkin_af( P, want_M = TRUE )
    )
    expect_equal( obj$coancestry, coancestry )
    expect_equal( obj$M, M_true )

    # repeat with P with missingness
    expect_silent(
        coancestry2 <- popkin_af( P_miss )
    )
    expect_equal( coancestry2, coancestry_miss )
    # repeat with want_M
    expect_silent(
        obj <- popkin_af( P_miss, want_M = TRUE )
    )
    expect_equal( obj$coancestry, coancestry_miss )
    expect_equal( obj$M, M_miss_true )

    # test transposed version
    expect_silent(
        obj2 <- popkin_af( t( P_miss ), want_M = TRUE, loci_on_cols = TRUE )
    )
    expect_equal( obj2, obj )
})

