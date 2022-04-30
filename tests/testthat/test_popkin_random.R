# the original popkin tests are all for recovering precomputed values
# here we test random genotype data

# first round: no missingness
n <- 100
m <- 1000
X <- matrix( rbinom( n*m, 2, 0.5 ), m, n )

test_that( "popkin_A works on random data without missingness", {
    # expected values
    A_exp <- crossprod( X-1 ) / m - 1
    M_exp <- matrix( m, n, n )
    
    # test function now!
    expect_silent( obj <- popkin_A( X ) )
    expect_equal( obj$A, A_exp )
    expect_equal( obj$M, M_exp )
})

test_that( "popkin_A MOR works on random data without missingness", {
    # expected values
    # need MAFs
    p <- rowMeans( X ) / 2
    # exclude fixed loci explicitly
    # (this is a list of loci to keep)
    indexes <- 0 < p & p < 1
    p <- p[ indexes ]
    Xnf <- X[ indexes, ] # nf = non-fixed
    mnf <- nrow( Xnf )
    # normalize data now
    # to be extra different from implementation, explicitly normalize and sum A over loci
    A_exp <- matrix( 0, n, n )
    M_exp <- matrix( mnf, n, n )
    for ( i in 1 : mnf ) {
        xi <- Xnf[ i, ]
        pi <- p[ i ]
        A_exp <- A_exp + ( tcrossprod( xi - 1 ) - 1 ) / ( pi * ( 1 - pi ) )
    }
    # final normalization
    A_exp <- A_exp / mnf
    
    # test function now!
    expect_silent( obj <- popkin_A( X, mean_of_ratios = TRUE ) )
    expect_equal( obj$A, A_exp )
    expect_equal( obj$M, M_exp )
})

# now add missingness
p_miss <- 0.1
X[ sample( n*m, n*m*p_miss ) ] <- NA

test_that( "popkin_A works on random data with missingness", {
    # expected values
    M_exp <- crossprod( !is.na( X ) )
    # modify copy of X, not original!
    X_tmp <- X - 1
    X_tmp[ is.na( X_tmp ) ] <- 0
    A_exp <- crossprod( X_tmp ) / M_exp - 1

    # test function now!
    expect_silent( obj <- popkin_A( X ) )
    expect_equal( obj$A, A_exp )
    expect_equal( obj$M, M_exp )
})

test_that( "popkin_A MOR works on random data with missingness", {
    # expected values
    # need MAFs
    p <- rowMeans( X, na.rm = TRUE ) / 2
    # exclude fixed loci explicitly
    # (this is a list of loci to keep)
    indexes <- 0 < p & p < 1
    p <- p[ indexes ]
    Xnf <- X[ indexes, ] # nf = non-fixed
    mnf <- nrow( Xnf )
    # use non-fixed loci only!
    M_exp <- crossprod( !is.na( Xnf ) )
    # normalize data now
    # to be extra different from implementation, explicitly normalize and sum A over loci
    A_exp <- matrix( 0, n, n )
    for ( i in 1 : mnf ) {
        pi <- p[ i ]
        xi <- Xnf[ i, ] - 1
        Mi <- tcrossprod( !is.na( xi ) ) # subtract ones only for loci counted!
        xi[ is.na( xi ) ] <- 0
        A_exp <- A_exp + ( tcrossprod( xi ) - Mi ) / ( pi * ( 1 - pi ) )
    }
    # final normalization
    A_exp <- A_exp / M_exp
    
    # test function now!
    expect_silent( obj <- popkin_A( X, mean_of_ratios = TRUE ) )
    expect_equal( obj$A, A_exp )
    expect_equal( obj$M, M_exp )
})
