# used in some examples
library(ape)

# loads Rdata matrices to test
load('Xs.RData')

test_that("phylo_max_edge works", {
    # plot examples with trees
    tree <- rtree( 3 )

    expect_silent(
        max_length <- phylo_max_edge( tree )
    )
    expect_true( is.numeric( max_length ) )
    expect_equal( length( max_length ), 1 )
    expect_true( max_length >= 0 )
    expect_true( max_length <= sum( tree$edge.length ) )
})

test_that("plot_phylo works", {
    # just in case some standalone defaults don't work as they should
    
    # set up a temporary path to write to
    fo <- tempfile('test-plot-phylo', fileext = '.pdf')
    # random sample tree
    tree <- rtree( 3 )
    
    # this should work
    pdf( fo, width = 14 )
    par(oma = c(0, 1.5, 0, 3))
    par(mar = c(0, 0, 2, 0) + 0.2)
    expect_silent( plot_phylo( tree ) )
    invisible( dev.off() )
    invisible( file.remove(fo) )
})

test_that("plot_popkin works", {
    # set up a temporary path to write to
    fo <- tempfile('test-plot-popkin', fileext = '.pdf')
    
    # singleton works
    pdf( fo )
    par(oma = c(0, 1.5, 0, 3))
    par(mar = c(0, 0, 2, 0) + 0.2)
    expect_silent( plot_popkin( inbr_diag(Phi), labs = subpops ) )
    invisible( dev.off() )
    invisible( file.remove(fo) )
    
    # list works
    pdf( fo, width = 14 ) # make wider
    par(oma = c(0, 1.5, 0, 3))
    par(mar = c(0, 0, 2, 0) + 0.2)
    expect_silent( plot_popkin( inbr_diag( list(Phi, Phi0) ) ) )
    invisible( dev.off() )
    invisible( file.remove(fo) )

    # list with NULL works
    pdf( fo, width = 14 )
    par(oma = c(0, 1.5, 0, 3))
    par(mar = c(0, 0, 2, 0) + 0.2)
    expect_silent( plot_popkin( inbr_diag( list(Phi, NULL) ) ) )
    invisible( dev.off() )
    invisible( file.remove(fo) )
    
    # test mismatch cases with NULLs and titles
    
    # this should work in old version
    pdf( fo, width = 14 )
    par(oma = c(0, 1.5, 0, 3))
    par(mar = c(0, 0, 2, 0) + 0.2)
    expect_silent( plot_popkin( inbr_diag( list(Phi, NULL) ), titles = c('a', 'b'), null_panel_data = TRUE ) )
    invisible( dev.off() )
    invisible( file.remove(fo) )
    
    # this should fail in old version
    pdf( fo, width = 14 )
    par(oma = c(0, 1.5, 0, 3))
    par(mar = c(0, 0, 2, 0) + 0.2)
    expect_error( plot_popkin( inbr_diag( list(Phi, NULL) ), titles = 'a', null_panel_data = TRUE ) )
    invisible( dev.off() )
    invisible( file.remove(fo) )
    
    # this should fail in new version
    pdf( fo, width = 14 )
    par(oma = c(0, 1.5, 0, 3))
    par(mar = c(0, 0, 2, 0) + 0.2)
    expect_error( plot_popkin( inbr_diag( list(Phi, NULL) ), titles = c('a', 'b'), null_panel_data = FALSE ) )
    invisible( dev.off() )
    invisible( file.remove(fo) )
    
    # this should work in new version
    pdf( fo, width = 14 )
    par(oma = c(0, 1.5, 0, 3))
    par(mar = c(0, 0, 2, 0) + 0.2)
    expect_silent( plot_popkin( inbr_diag( list(Phi, NULL) ), titles = 'a', null_panel_data = FALSE ) )
    invisible( dev.off() )
    invisible( file.remove(fo) )

    # plot examples with trees
    tree <- rtree( 3 )

    # this should work
    pdf( fo, width = 14 )
    par(oma = c(0, 1.5, 0, 3))
    par(mar = c(0, 0, 2, 0) + 0.2)
    expect_silent( plot_popkin( list( tree, inbr_diag( Phi ), NULL) ) )
    invisible( dev.off() )
    invisible( file.remove(fo) )
    
})

test_that("plot_admix works", {
    # set up a temporary path to write to
    fo <- tempfile('test-plot-admix', fileext = '.pdf')

    # same example from vignette
    # some subpopulation sizes
    n1 <- 10
    n2 <- n1
    n3 <- 20
    # construct unadmixed individuals
    Q_pop1 <- matrix( c( 1, 0 ), ncol = 2, nrow = n1, byrow = TRUE )
    Q_pop2 <- Q_pop1[ , 2:1 ] # reverse columns
    # random admixture proportions
    Q_admix <- rev( sort( runif( n3 ) ) )
    Q_admix <- cbind( Q_admix, 1 - Q_admix )
    # combine data for all populations
    # add some informative row/col names too
    rownames( Q_pop1 ) <- paste0( 'S1-ind', 1 : n1 )
    rownames( Q_pop2 ) <- paste0( 'S2-ind', 1 : n2 )
    rownames( Q_admix ) <- paste0( 'S3-ind', 1 : n3 )
    Q <- rbind( Q_pop1, Q_pop2, Q_admix )
    colnames( Q ) <- c('A1', 'A2')
    # create subpopulation labels for later
    labs1 <- c( rep.int( 'S1', n1 ), rep.int( 'S2', n2 ), rep.int( 'S3', n3 ) )
    labs2 <- c( rep.int( 'Unadmixed', n1+n2 ), rep.int( 'Admixed', n3 ) )

    # no labels/names
    pdf( fo )
    par(mar = c(1, 3, 1, 0.2) + 0.2)
    expect_silent(
        plot_admix( Q )
    )
    invisible( dev.off() )
    invisible( file.remove(fo) )

    # names
    pdf( fo )
    par(mar = c(6, 3, 1, 0.2) + 0.2)
    expect_silent( 
        plot_admix( Q, names = TRUE, xlab_line = 5 )
    )
    invisible( dev.off() )
    invisible( file.remove(fo) )

    # labels version
    pdf( fo )
    par(mar = c(3, 3, 1, 0.2) + 0.2)
    expect_silent( 
        plot_admix( Q, labs = labs1, xlab_line = 2 )
    )
    invisible( dev.off() )
    invisible( file.remove(fo) )

    # labels version, legend omitted
    pdf( fo )
    par(mar = c(3, 3, 1, 0.2) + 0.2)
    expect_silent( 
        plot_admix( Q, labs = labs1, xlab_line = 2, leg_omit = TRUE )
    )
    invisible( dev.off() )
    invisible( file.remove(fo) )

})



### these don't create graphics, but are plotting-related


test_that( 'admix_order_cols works', {
    # write a toy Q where the desired order is evident
    Q <- matrix(
        c(
            0.1, 0.8, 0.1,
            0.1, 0.7, 0.2,
            0.0, 0.4, 0.6,
            0.0, 0.3, 0.7,
            0.9, 0.0, 0.1
        ),
        nrow = 5,
        ncol = 3,
        byrow = TRUE
    )
    order_exp <- c(2, 3, 1)
    expect_silent( 
        order_obs <- admix_order_cols( Q )
    )
    expect_equal( order_obs, order_exp )
    
})

test_that( 'make_unique works', {
    # simple function just feed it some quick expectations
    x <- 1:10
    # was already unique
    expect_equal( make_unique( x ), x )
    # simple repeats
    expect_equal( make_unique( c('A', 'A', 'B') ), c('A1', 'A2', 'B') )
    expect_equal( make_unique( c('A', 'B', 'B') ), c('A', 'B1', 'B2') )
    # more complex repeats
    expect_equal( make_unique( c('A', 'Z', 'B', 'Z', 'A') ), c('A1', 'Z1', 'B', 'Z2', 'A2') )
})

test_that( 'admix_mean_subpops_single and admix_label_cols work', {
    # same toy data, but now add labels for individuals
    Q <- matrix(
        c(
            0.1, 0.8, 0.1,
            0.1, 0.7, 0.2,
            0.0, 0.4, 0.6,
            0.0, 0.3, 0.7,
            0.9, 0.0, 0.1
        ),
        nrow = 5,
        ncol = 3,
        byrow = TRUE
    )
    labs <- c('X', 'X', 'Y', 'Y', 'Z')
    
    # repeat "single" test for each column
    for ( k in 1 : ncol(Q) ) {
        q <- Q[, k]
        # run code
        admix_mean_obs <- admix_mean_subpops_single( q, labs )
        # explicit calculation hardcoded to the above labels
        admix_mean_exp <- c( mean( q[1:2] ), mean( q[3:4] ), q[5] )
        names( admix_mean_exp ) <- c('X', 'Y', 'Z')
        expect_equal( admix_mean_obs, admix_mean_exp )
    }

    # now test the overall procedure of assigning labels to columns
    col_labs_obs <- admix_label_cols( Q, labs )
    col_labs_exp <- c('Z', 'X', 'Y')
    expect_equal( col_labs_obs, col_labs_exp )

    # a case where there's repeats that we make unique
    # (there's 3 ancestries but 2 labels, so obviously there's repeats)
    labs <- c('A', 'A', 'A', 'A', 'B')
    col_labs_obs <- admix_label_cols( Q, labs )
    col_labs_exp <- c('B', 'A1', 'A2')
    expect_equal( col_labs_obs, col_labs_exp )
})

