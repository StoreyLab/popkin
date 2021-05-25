context('popkin_Rdata')

# used in some examples
library(ape)

# loads Rdata matrices to test
load('Xs.RData')

# start with lower-level/internal tests, more informative that higher-level function errors

test_that( "M is correct", {
    # the M used in other tests was pre-calculated by popkin_A
    # here we calculate it from X the way we do it elsewhere (including popkin_af), make sure it agrees

    # popkin removes loci that are fixed (not in M count)
    x_hat <- rowMeans( X, na.rm = TRUE )
    # indexes to keep
    indexes <- !is.na( x_hat ) & x_hat > 0 & x_hat < 2
    # filter loci
    X2 <- X[ indexes, ]

    # now calculate M
    M_direct <- crossprod( !is.na( X2 ) )
    # in this case M_direct gets labels from X, let's remove those (M doesn't have them)
    dimnames( M_direct ) <- NULL
    # matching requires labels to agree
    expect_equal( M, M_direct )
})

test_that("validate_kinship works", {
    # try default "error" versions alongside logical versions
    
    # validate positive examples
    expect_silent( validate_kinship( Phi ) )
    expect_silent( validate_kinship( Phi0 ) )
    expect_silent( validate_kinship( A, name = 'A' ) ) # not real kinship but satisfies requirements
    expect_true( validate_kinship( Phi, logical = TRUE ) )
    expect_true( validate_kinship( Phi0, logical = TRUE ) )
    expect_true( validate_kinship( A, name = 'A', logical = TRUE ) ) # not real kinship but satisfies requirements

    # negative examples
    # dies if input is missing
    expect_error( validate_kinship() )
    expect_error( validate_kinship( logical = TRUE ) ) # still dies here
    
    # NULL values should fail (important for `plot_popkin`)
    expect_error( validate_kinship( NULL ) )
    expect_false( validate_kinship( NULL, logical = TRUE ) )
    
    # and if input is not a matrix
    expect_error( validate_kinship( 1:5 ) )
    expect_false( validate_kinship( 1:5, logical = TRUE ) )
    
    # and for non-numeric matrices
    char_mat <- matrix(c('a', 'b', 'c', 'd'), nrow=2)
    expect_error( validate_kinship( char_mat ) )
    expect_false( validate_kinship( char_mat, logical = TRUE ) )
    
    # and non-square matrices
    non_kinship <- matrix(1:2, nrow=2)
    expect_error( validate_kinship( non_kinship ) )
    expect_false( validate_kinship( non_kinship, logical = TRUE ) )
    
    # and non-symmetric matrices
    # (also test `sym = FALSE` option that lets this case pass)
    non_kinship <- matrix(1:4, nrow=2)
    expect_error( validate_kinship( non_kinship ) )
    expect_silent( validate_kinship( non_kinship, sym = FALSE ) )
    expect_false( validate_kinship( non_kinship, logical = TRUE ) )
    expect_true( validate_kinship( non_kinship, sym = FALSE, logical = TRUE ) )
})


test_that("get_mem_lim returns positive numbers", {
    mem <- get_mem_lim()
    expect_equal(class(mem), 'numeric')
    expect_true(mem > 0)
})

test_that("solve_m_mem_lim works", {
    # mem is in GB here, inferred from system if missing
    n <- 1000
    m <- 100000
    # create failures on purpose
    expect_error( solve_m_mem_lim() ) # mandatory values missing
    expect_error( solve_m_mem_lim( mem = 1, n = n ) ) # counts of m matrices are zero
    expect_error( solve_m_mem_lim( mem = 1/GB, n = n, vec_m = 1, mat_n_n = 1 ) ) # memory is too low for the n given
    # now reasonable success cases
    solve_m_mem_lim_TESTER <- function(
                                       mat_m_n = 0,
                                       mat_n_n = 0,
                                       vec_m = 0,
                                       vec_n = 0,
                                       mem = NA
                                       ) {
        data <- solve_m_mem_lim(
            mem = mem,
            n = n,
            m = m,
            mat_m_n = mat_m_n,
            mat_n_n = mat_n_n,
            vec_m = vec_m,
            vec_n = vec_n
        )
        expect_equal(class(data), 'list')
        expect_equal(length(data), 3)
        expect_equal(names(data), c('m_chunk', 'mem_chunk', 'mem_lim'))
        # data sets bounds too!
        expect_true( data$m_chunk > 0 )
        expect_true( data$m_chunk <= m )
        expect_true( data$mem_chunk > 0 )
        # don't test this if no memory was specified
        if ( !is.na(mem) )
            expect_true( data$mem_chunk/GB < mem ) # strict ineq because of factor = 0.7
    }
    
    # repeat with various settings
    # we can't have both mat_m_n and vec_m be zero, but all other cases are tested (12 cases)
    # here we set memory manually to 1GB
    solve_m_mem_lim_TESTER(mem = 1, mat_m_n = 1)
    solve_m_mem_lim_TESTER(mem = 1, mat_m_n = 1, mat_n_n = 1)
    solve_m_mem_lim_TESTER(mem = 1, mat_m_n = 1, vec_n = 1)
    solve_m_mem_lim_TESTER(mem = 1, mat_m_n = 1, vec_n = 1, mat_n_n = 1)
    solve_m_mem_lim_TESTER(mem = 1, vec_m = 1)
    solve_m_mem_lim_TESTER(mem = 1, vec_m = 1, mat_n_n = 1)
    solve_m_mem_lim_TESTER(mem = 1, vec_m = 1, vec_n = 1)
    solve_m_mem_lim_TESTER(mem = 1, vec_m = 1, vec_n = 1, mat_n_n = 1)
    solve_m_mem_lim_TESTER(mem = 1, mat_m_n = 1, vec_m = 1)
    solve_m_mem_lim_TESTER(mem = 1, mat_m_n = 1, vec_m = 1, mat_n_n = 1)
    solve_m_mem_lim_TESTER(mem = 1, mat_m_n = 1, vec_m = 1, vec_n = 1)
    solve_m_mem_lim_TESTER(mem = 1, mat_m_n = 1, vec_m = 1, vec_n = 1, mat_n_n = 1)
    # repeat inferring memory from system
    solve_m_mem_lim_TESTER(mat_m_n = 1)
    solve_m_mem_lim_TESTER(mat_m_n = 1, mat_n_n = 1)
    solve_m_mem_lim_TESTER(mat_m_n = 1, vec_n = 1)
    solve_m_mem_lim_TESTER(mat_m_n = 1, vec_n = 1, mat_n_n = 1)
    solve_m_mem_lim_TESTER(vec_m = 1)
    solve_m_mem_lim_TESTER(vec_m = 1, mat_n_n = 1)
    solve_m_mem_lim_TESTER(vec_m = 1, vec_n = 1)
    solve_m_mem_lim_TESTER(vec_m = 1, vec_n = 1, mat_n_n = 1)
    solve_m_mem_lim_TESTER(mat_m_n = 1, vec_m = 1)
    solve_m_mem_lim_TESTER(mat_m_n = 1, vec_m = 1, mat_n_n = 1)
    solve_m_mem_lim_TESTER(mat_m_n = 1, vec_m = 1, vec_n = 1)
    solve_m_mem_lim_TESTER(mat_m_n = 1, vec_m = 1, vec_n = 1, mat_n_n = 1)

    # test that we now successfully avoid a particularly bad integer overflow error
    # caused internally by an `n*n` term, so only n has to be huge to cause it
    # increase memory requested so that doesn't cause problems
    n <- 50000L # pass as integer!
    mem <- 32
    solve_m_mem_lim_TESTER(mem = mem, mat_m_n = 1)
    solve_m_mem_lim_TESTER(mem = mem, mat_m_n = 1, mat_n_n = 1)
    solve_m_mem_lim_TESTER(mem = mem, mat_m_n = 1, vec_n = 1)
    solve_m_mem_lim_TESTER(mem = mem, mat_m_n = 1, vec_n = 1, mat_n_n = 1)
    solve_m_mem_lim_TESTER(mem = mem, vec_m = 1)
    solve_m_mem_lim_TESTER(mem = mem, vec_m = 1, mat_n_n = 1)
    solve_m_mem_lim_TESTER(mem = mem, vec_m = 1, vec_n = 1)
    solve_m_mem_lim_TESTER(mem = mem, vec_m = 1, vec_n = 1, mat_n_n = 1)
    solve_m_mem_lim_TESTER(mem = mem, mat_m_n = 1, vec_m = 1)
    solve_m_mem_lim_TESTER(mem = mem, mat_m_n = 1, vec_m = 1, mat_n_n = 1)
    solve_m_mem_lim_TESTER(mem = mem, mat_m_n = 1, vec_m = 1, vec_n = 1)
    solve_m_mem_lim_TESTER(mem = mem, mat_m_n = 1, vec_m = 1, vec_n = 1, mat_n_n = 1)
})

test_that("function returns precomputed values: weights_subpops", {
    expect_equal(weights_subpops(subpops0), w0)
    expect_equal(weights_subpops(subpops), w)
    # make sure dimensions match
    expect_equal(length(w0), nrow(Phi0))
    expect_equal(length(w), nrow(Phi))
    # test the basic qualities of weights
    expect_equal(sum(w0), 1)
    expect_equal(sum(w), 1)
    expect_true(all(w0 > 0))
    expect_true(all(w > 0))
    expect_true(all(w0 < 1))
    expect_true(all(w < 1))
})

test_that("function returns precomputed values: popkin_A", {
    # standard test (X is integer)
    expect_silent( obj <- popkin_A( X ) )
    expect_equal( obj$A, A )
    expect_equal( obj$M, M )
    
    # turns X into doubles
    expect_silent( obj <- popkin_A( X + 0 ) )
    expect_equal( obj$A, A )
    expect_equal( obj$M, M )

    # reflect, keep X integer
    expect_silent( obj <- popkin_A( 2L - X ) )
    expect_equal( obj$A, A )
    expect_equal( obj$M, M )

    # reflect and turn X doubles too
    expect_silent( obj <- popkin_A( 2 - X ) )
    expect_equal( obj$A, A )
    expect_equal( obj$M, M )

    # make sure that non-default m_chunk_max version works
    # this tests the edge case `m_chunk_max = 1` in particular
    expect_silent( obj <- popkin_A( X, m_chunk_max = 1 ) )
    expect_equal( obj$A, A )
    expect_equal( obj$M, M )
    
    # these should be square matrices
    expect_equal(nrow(A), ncol(A))
    expect_equal(nrow(M), ncol(M))

    # M (pairwise sample sizes, so excluding NA pairs) must satisfy obvious range limits
    expect_true( min(M) >= 0 )
    expect_true( max(M) <= nrow(X) )
})

test_that("function returns precomputed values: popkin_A_min_subpops", {
    # trigger expected errors
    # missing A
    expect_error( popkin_A_min_subpops() )
    # subpops with wrong length
    expect_error( popkin_A_min_subpops( A, subpops[1:2] ) )
    # subpops with a single subpop really (but right length)
    n_ind <- length( subpops )
    subpops_one <- rep.int( 1, n_ind )
    expect_error( popkin_A_min_subpops( A, subpops_one ) )
    # but two subpopulations should be ok
    n2 <- round( n_ind / 2 )
    subpops_two <- c( rep.int( 1, n2 ), rep.int( 2, n_ind - n2 ) )
    expect_silent( popkin_A_min_subpops( A, subpops_two ) )
    
    # now positive examples
    expect_equal(popkin_A_min_subpops(A), min(A))
    expect_equal(popkin_A_min_subpops(A), Amin0)
    expect_equal(popkin_A_min_subpops(A, subpops0), Amin0)
    expect_equal(popkin_A_min_subpops(A, subpops), Amin)
})

# higher-level tests now!

test_that("function returns precomputed values: popkin", {
    # cause problems on purpose
    expect_error( popkin() )
    
    # good runs
    expect_equal(popkin(X), Phi0)
    expect_equal(popkin(X, subpops0), Phi0)
    expect_equal(popkin(X, subpops), Phi)
    expect_equal(popkin(X+0, subpops), Phi)
    expect_equal(popkin(2L-X, subpops), Phi)
    expect_equal(popkin(2-X, subpops), Phi)

    # test that non-default m_chunk_max options don't fail
    expect_equal(popkin(X, subpops, m_chunk_max = 1), Phi)
    
    # test want_M in these two cases only
    # though kinship is scaled differently in each version, M is the same for both
    expect_silent( obj <- popkin( X, want_M = TRUE ) )
    expect_equal( obj$kinship, Phi0 )
    expect_equal( obj$M, M )
    
    expect_silent( obj <- popkin( X, subpops, want_M = TRUE ) )
    expect_equal( obj$kinship, Phi )
    expect_equal( obj$M, M )
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
    # type of return value
    expect_equal(length(fstW), 1)
    expect_equal(length(fst), 1)
    # Fst inequalities
    expect_true(fstW >= 0)
    expect_true(fst >= 0)
    expect_true(fstW <= 1)
    expect_true(fst <= 1)
})

test_that("fst with local inbreeding works", {
    # construct a random fake local kinship matrix
    n_ind <- nrow(Phi)
    kinship_local <- matrix(
        runif( n_ind * n_ind ),
        nrow = n_ind,
        ncol = n_ind
    )
    # make pos-def by cross-multiplying with itself
    # normalize to keep in [0,1]
    kinship_local <- crossprod( kinship_local ) / n_ind
    # now create "total matrix" by combining effects
    # however, do this in coancestry mode (inbr diag)
    # - kinship_local has not been transformed yet so only Phi needs this
    kinship_total <- inbr_diag(Phi) + kinship_local - inbr_diag(Phi) * kinship_local
    # turn both local and total into proper kinship by altering diagonal
    diag( kinship_local ) <- ( 1 + diag( kinship_local ) ) / 2
    diag( kinship_total ) <- ( 1 + diag( kinship_total ) ) / 2
    # now estimate FST
    # correcting appropriately in all cases should result in the same FST as before...
    expect_equal(fst(kinship_total, x_local = kinship_local), fst)
    expect_equal(fst(kinship_total, x_local = kinship_local, w0), fst)
    expect_equal(fst(kinship_total, x_local = kinship_local, w), fstW)
    # test vector versions of both, and each separately
    expect_equal(fst(inbr(kinship_total), x_local = kinship_local), fst)
    expect_equal(fst(inbr(kinship_total), x_local = kinship_local, w0), fst)
    expect_equal(fst(inbr(kinship_total), x_local = kinship_local, w), fstW)
    expect_equal(fst(inbr(kinship_total), x_local = inbr(kinship_local)), fst)
    expect_equal(fst(inbr(kinship_total), x_local = inbr(kinship_local), w0), fst)
    expect_equal(fst(inbr(kinship_total), x_local = inbr(kinship_local), w), fstW)
    expect_equal(fst(kinship_total, x_local = inbr(kinship_local)), fst)
    expect_equal(fst(kinship_total, x_local = inbr(kinship_local), w0), fst)
    expect_equal(fst(kinship_total, x_local = inbr(kinship_local), w), fstW)
})

test_that("function returns precomputed values: inbr", {
    expect_equal(inbr(Phi), inbr)
})

test_that("function returns precomputed values: pwfst", {
    expect_equal(pwfst(Phi), pwF)
    expect_equal(pwfst(Phi0), pwF)
    expect_equivalent(diag(pwF), rep.int(0, nrow(pwF))) # test that diagonal is zero ("equivalent" ignores label mismatches)
    expect_true(max(pwF) <= 1)
    # note estimates may be slightly negative though
})

test_that("inbr_diag works", {
    # dies when kinship matrix is missing
    expect_error( inbr_diag() )
    # make sure precomputed values match
    expect_equal( inbr_diag(Phi), PhiInbr )
    # test list version (just duplicates things)
    expect_equal( inbr_diag(list(Phi, Phi)), list(PhiInbr, PhiInbr) )
})

test_that("mean_kinship works", {
    # dies when kinship matrix is missing
    expect_error(mean_kinship())
    # equals ordinary mean without weights
    expect_equal(mean_kinship(Phi), mean(Phi))
    expect_equal(mean_kinship(Phi0), mean(Phi0))
    # equals expected double weight formula under non-trivial weights
    expect_equal(mean_kinship(Phi, w0), mean(Phi)) # this is still uniform weights actually
    # draw random uniform for a very non-trivial test
    weights <- runif(ncol(Phi))
    # normalize for comparison
    weights <- weights / sum(weights)
    # actually compare formulas
    expect_equal(mean_kinship(Phi, weights), drop(weights %*% Phi %*% weights))
})

test_that("n_eff works", {

    # a small toy case!
    F <- 0.4 # for this case to produce negative weights, need F > 1/3
    K <- (1+F)/2 # the self kinship values
    Phi3 <- matrix(c(K,F,0,F,K,F,0,F,K),nrow=3)
    w3 <- c(0.4, 0.2, 0.4) # dummy weights for a test

    # default case is sum of elements of inverse matrix!
    expect_equal(n_eff(Phi3, nonneg=FALSE)$n_eff, sum(solve(Phi3)))
    # non-max case returns inverse of mean kinship (with uniform weights)
    expect_equal(n_eff(Phi3, max=FALSE)$n_eff, 1/mean(Phi3))
    # non-max case returns inverse of weighted mean kinship (with provided non-uniform weights)
    expect_equal(n_eff(Phi3, max=FALSE, w=w3)$n_eff, 1/drop(w3 %*% Phi3 %*% w3))
    # test that gradient descent is the default
    expect_equal(n_eff(Phi3)$n_eff, n_eff(Phi3, algo='g')$n_eff)
    # basic inequalities
    expect_true(n_eff(Phi3)$n_eff >= 1) # min possible value
    expect_true(n_eff(Phi3, algo='g')$n_eff >= 1) # min possible value
    expect_true(n_eff(Phi3, algo='n')$n_eff >= 1) # min possible value
    expect_true(n_eff(Phi3, algo='h')$n_eff >= 1) # min possible value
    expect_true(n_eff(Phi3, nonneg=FALSE)$n_eff >= 1) # min possible value
    expect_true(n_eff(Phi3, max=FALSE)$n_eff >= 1) # min possible value
    expect_true(n_eff(Phi3, max=FALSE, w=w3)$n_eff >= 1) # min possible value
    expect_true(n_eff(Phi3)$n_eff <= 2*nrow(Phi3)) # max possible value
    expect_true(n_eff(Phi3, algo='g')$n_eff <= 2*nrow(Phi3)) # max possible value
    expect_true(n_eff(Phi3, algo='n')$n_eff <= 2*nrow(Phi3)) # max possible value
    expect_true(n_eff(Phi3, algo='h')$n_eff <= 2*nrow(Phi3)) # max possible value
    expect_true(n_eff(Phi3, nonneg=FALSE)$n_eff <= 2*nrow(Phi3)) # max possible value
    expect_true(n_eff(Phi3, nonneg=FALSE)$n_eff >= n_eff(Phi3)$n_eff) # the max n_eff should really be larger than a non-max case (non-optimal weights)
    expect_true(n_eff(Phi3, nonneg=FALSE)$n_eff >= n_eff(Phi3, max=FALSE)$n_eff) # the max n_eff should really be larger than a non-max case (non-optimal weights)
    expect_true(n_eff(Phi3, nonneg=FALSE)$n_eff >= n_eff(Phi3, max=FALSE, w=w3)$n_eff) # the max n_eff should really be larger than a non-max case (non-optimal weights)
    expect_true(n_eff(Phi3)$n_eff >= n_eff(Phi3, max=FALSE)$n_eff) # hope the numeric max with non-negative weights gives a larger value than a uniform weights estimate (actually not gauranteed to perform better)
    expect_true(n_eff(Phi3, algo='g')$n_eff >= n_eff(Phi3, max=FALSE)$n_eff)
    expect_true(n_eff(Phi3, algo='n')$n_eff >= n_eff(Phi3, max=FALSE)$n_eff)
    expect_true(n_eff(Phi3, algo='h')$n_eff >= n_eff(Phi3, max=FALSE)$n_eff)

    # construct some other toy data tests
    # check that we get 2*n on unstructured kinship matrix
    expect_equal(n_eff(diag(1/2, 10, 10))$n_eff, 20)
    # check that we get K on extreme Fst=1 independent subpopulations
    # NOTE: FAILS because it's not invertible... what should we do here???
    #    expect_equal(n_eff(matrix(c(1, 1, 0, 1, 1, 0, 0, 0, 1), nrow=3)$n_eff), 2)
    # maybe construct Fst<1 case and compare to theoretical expectation there

    # normally weights are returned too, but their values are less constrained (by construction they sum to one, that's it!)
    # do test both max versions since n_eff is not directly constructed from the weights (test that it's what it should be)
    # test outputs in that setting now
    obj <- n_eff(Phi3, algo='g') # this tests Gradient version
    expect_equal(class(obj), 'list') # return class is list
    expect_equal(length(obj), 2) # only have two elements
    # roughly retest that first element is an n_eff
    expect_true(obj$n_eff >= 1) # min possible value
    expect_true(obj$n_eff <= 2*nrow(Phi3)) # max possible value
    # roughly test weights
    expect_equal(sum(obj$weights), 1) # verify that weights sum to 1
    expect_equal(obj$n_eff, 1/mean_kinship(Phi3, obj$weights)) # verify that n_eff has the value it should have given the weights
    
    obj <- n_eff(Phi3, algo='n') # this tests Newton version
    expect_equal(class(obj), 'list') # return class is list
    expect_equal(length(obj), 2) # only have two elements
    # roughly retest that first element is an n_eff
    expect_true(obj$n_eff >= 1) # min possible value
    expect_true(obj$n_eff <= 2*nrow(Phi3)) # max possible value
    # roughly test weights
    expect_equal(sum(obj$weights), 1) # verify that weights sum to 1
    expect_equal(obj$n_eff, 1/mean_kinship(Phi3, obj$weights)) # verify that n_eff has the value it should have given the weights
    
    obj <- n_eff(Phi3, algo='h') # this tests Heuristic version
    expect_equal(class(obj), 'list') # return class is list
    expect_equal(length(obj), 2) # only have two elements
    # roughly retest that first element is an n_eff
    expect_true(obj$n_eff >= 1) # min possible value
    expect_true(obj$n_eff <= 2*nrow(Phi3)) # max possible value
    # roughly test weights
    expect_equal(sum(obj$weights), 1) # verify that weights sum to 1
    expect_equal(obj$n_eff, 1/mean_kinship(Phi3, obj$weights)) # verify that n_eff has the value it should have given the weights
    
    obj <- n_eff(Phi3, nonneg=FALSE) # this tests optimal version (with possibly negative weights)
    expect_equal(class(obj), 'list') # return class is list
    expect_equal(length(obj), 2) # only have two elements
    # roughly retest that first element is an n_eff
    expect_true(obj$n_eff >= 1) # min possible value
    expect_true(obj$n_eff <= 2*nrow(Phi3)) # max possible value
    # roughly test weights
    expect_equal(sum(obj$weights), 1) # verify that weights sum to 1
    expect_equal(obj$n_eff, 1/mean_kinship(Phi3, obj$weights)) # verify that n_eff has the value it should have given the weights
    expect_true(min(obj$weights) < 0) # this example must have negative weights, or it's useless!
})

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

test_that( "true M are ok", {
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
