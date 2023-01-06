#' Estimate coancestry from an allele frequency matrix and subpopulation assignments
#'
#' Given the individual-specific allele frequencies of `n` individuals, this function returns the `n`-by-`n` coancestry matrix.
#' This function is the analog of [popkin()] for allele frequencies rather than genotypes, and as a consequence estimates coancestry instead of kinship.
#' These coancestry estimates are unbiased if the true allele frequencies are provided, but may be less accurate when the allele frequencies themselves are estimated.
#' This function is intended for cases where allele frequencies, but not individual genotypes, are available; otherwise it is best to use the individual genotypes and [popkin()].
#' An application of interest is the allele frequency matrices from admixture models, in which case the columns correspond to subpopulations rather than individuals, and `subpops = NULL` is an acceptable choice.
#'
#' @inheritParams popkin
#' @param P `m`-by-`n` matrix of individual-specific allele frequencies, which should have values between `[0, 1]` (range is not strictly required) or `NA` for missing data.
#' @param loci_on_cols If `TRUE`, `P` has loci on columns and individuals on rows; if `FALSE` (default), loci are on rows and individuals on columns.
#'
#' @return If `want_M = FALSE`, returns the estimated `n`-by-`n` coancestry matrix only.
#' If `P` has names for the individuals, they will be copied to the rows and columns of this coancestry matrix.
#' If `want_M = TRUE`, a named list is returned, containing:
#'
#' - `coancestry`: the estimated `n`-by-`n` coancestry matrix
#' - `M`: the `n`-by-`n` matrix of non-missing pair counts (see `want_M` option).
#'
#' @examples
#' # a matrix of random uniform allele frequencies
#' # (unstructured, unlike real data)
#' P <- matrix( runif( 9 ), nrow = 3 )
#'
#' coancestry <- popkin_af( P )
#'
#' @seealso
#' [popkin()] for kinship estimation from genotype matrices.
#'
#' @export
popkin_af <- function(
                      P,
                      subpops = NULL,
                      loci_on_cols = FALSE,
                      mem_factor = 0.7,
                      mem_lim = NA,
                      want_M = FALSE,
                      m_chunk_max = 1000
                      ) {
    # check for required data
    if ( missing( P ) )
        stop( '`P` is required!' )

    # get dimensions
    # number of loci and "individuals"
    # (in expected admixture-like data, each individual is actually a subpopulation, meh)
    if ( loci_on_cols ) {
        m_loci <- ncol( P )
        n_ind <- nrow( P )
        names_P <- rownames( P )
    } else {
        m_loci <- nrow( P )
        n_ind <- ncol( P )
        names_P <- colnames( P )
    }
    
    # test coherence between subpops and P
    if ( !is.null( subpops ) ) {
        if ( length( subpops ) != n_ind ) {
            stop('`subpops` length (', length( subpops ), ') disagreed with number of individuals in `P` (', n_ind, ')!')
        }
    }

    # initialize desired matrices
    A <- matrix( 0, nrow = n_ind, ncol = n_ind )
    M <- matrix( 0, nrow = n_ind, ncol = n_ind )
    
    # transfer names from P to A if present
    # this will carry over all the way to the final coancestry matrix!
    # (M need not have names at all)
    if ( !is.null( names_P ) ) {
        colnames( A ) <- names_P
        rownames( A ) <- names_P
    }
    
    # infer the number of SNPs to break data into, since we're limited by memory
    # given fixed n, solve for m:
    # get maximum m (number of SNPs) given n and the memory requested
    data <- solve_m_mem_lim(
        n = n_ind,
        m = m_loci,
        mat_m_n = 2, # P x2
        mat_n_n = 2, # A + M (1 + 1)
        mem = mem_lim,
        mem_factor = mem_factor
    )
    m_chunk <- data$m_chunk
    # cap value to a nice performing value (very good speed, minimal memory)
    if ( m_chunk > m_chunk_max )
        m_chunk <- m_chunk_max

    # navigate chunks
    i_chunk <- 1 # start of first chunk (needed for matrix inputs only; as opposed to function inputs)
    while (TRUE) { # start an infinite loop, break inside as needed

        # this means all SNPs have been covered!
        if (i_chunk > m_loci)
            break
        
        # range of SNPs to extract in this chunk
        indexes_loci_chunk <- i_chunk : min( i_chunk + m_chunk - 1, m_loci )
        
        if ( loci_on_cols ) {
            Pi <- t( P[, indexes_loci_chunk, drop = FALSE] ) # transpose for our usual setup
        } else  {
            Pi <- P[indexes_loci_chunk, , drop = FALSE]
        }
        i_chunk <- i_chunk + m_chunk # update starting point for next chunk! (overshoots at the end, that's ok)
        
        # construct M matrix, to normalize according to missingness
        # NOTES: 
        # - there are no "fixed" loci in the allele frequency world (exact zeroes or ones are practically impossible), so only missingness matters (compared to popkin, where we can easily draw data that is all zeroes or all 2s)
        # - completely missing loci don't contribute to M and also don't to A as is, so there's no need to remove them
        # this is the actual M matrix
        Mi <- crossprod( !is.na( Pi ) )
        
        # now we center Pi
        # original Pi is no longer needed, let's just overwrite it
        Pi <- Pi - 1/2
        # now zero out the NAs
        Pi[ is.na( Pi ) ] <- 0

        # now we can calculate A
        Ai <- crossprod( Pi )
        # increment each part
        A <- A + Ai # increment sum of A_{ijk}'s
        M <- M + Mi # increment M's too
    }
    
    # complete average and finish the calculation of A
    A <- A / M - 1/4
    
    # the denominator is a simple average, a scalar shared by all individuals
    A_min <- popkin_A_min_subpops( A, subpops )
    
    # this is the coancestry estimate
    coancestry <- 1 - A / A_min
    
    # figure out what to return
    if ( want_M ) {
        return(
            list(
                coancestry = coancestry,
                M = M
            )
        )
    } else {
        return( coancestry )
    }
}



# this is a primitive version of popkin for an allele frequency matrix
# intended for the P matrices from admixture and such, to estimate Psi matrices from them
# good things (vs official popkin for genotypes):
# - works on continuous data (official genotype version requires integers in 0, 1, 2)
# - as this is for allele frequencies covariance is *coancestry*, not kinship (no need for further diagonal transforms)
# bad things (this basic version only):
# - does not control memory usage
# - does not handle any missing data
popkin_af_basic <- function(
                            P,
                            subpops = NULL,
                            loci_on_cols = FALSE,
                            want_M = FALSE
                            ) {
    # check for required data
    if ( missing( P ) )
        stop( '`P` is required!' )

    # get dimensions
    # number of loci and "individuals"
    # (in expected admixture-like data, each individual is actually a subpopulation, meh)
    if ( loci_on_cols ) {
        m_loci <- ncol( P )
        n_ind <- nrow( P )
        names_P <- rownames( P )
    } else {
        m_loci <- nrow( P )
        n_ind <- ncol( P )
        names_P <- colnames( P )
    }
    
    # test coherence between subpops and P
    if ( !is.null( subpops ) ) {
        if ( length( subpops ) != n_ind ) {
            stop('`subpops` length (', length( subpops ), ') disagreed with number of individuals in `P` (', n_ind, ')!')
        }
    }
    
    # construct equivalent to A matrix
    if ( loci_on_cols ) {
        A <- tcrossprod( P - 1/2 ) / m_loci - 1/4
    } else {
        A <- crossprod( P - 1/2 ) / m_loci - 1/4
    }
    
    # transfer names from P to A if present
    # this will carry over all the way to the final coancestry matrix!
    # (M need not have names at all)
    if ( !is.null( names_P ) ) {
        colnames( A ) <- names_P
        rownames( A ) <- names_P
    }
    
    # the denominator is a simple average, a scalar shared by all individuals
    A_min <- popkin_A_min_subpops( A, subpops )
    
    # this is the coancestry estimate
    coancestry <- 1 - A / A_min
    
    # figure out what to return
    if ( want_M ) {
        return(
            list(
                coancestry = coancestry,
                M = m_loci # dummy version
            )
        )
    } else {
        return( coancestry )
    }
}

# extend basic case to handle missingness
# bad things remaining:
# - does not control memory usage
popkin_af_basic_na <- function(
                               P,
                               subpops = NULL,
                               loci_on_cols = FALSE,
                               want_M = FALSE
                               ) {
    # check for required data
    if ( missing( P ) )
        stop( '`P` is required!' )

    # get dimensions
    # number of loci and "individuals"
    # (in expected admixture-like data, each individual is actually a subpopulation, meh)
    if ( loci_on_cols ) {
        m_loci <- ncol( P )
        n_ind <- nrow( P )
        names_P <- rownames( P )
    } else {
        m_loci <- nrow( P )
        n_ind <- ncol( P )
        names_P <- colnames( P )
    }
    
    # test coherence between subpops and P
    if ( !is.null( subpops ) ) {
        if ( length( subpops ) != n_ind ) {
            stop('`subpops` length (', length( subpops ), ') disagreed with number of individuals in `P` (', n_ind, ')!')
        }
    }

    # construct M matrix, to normalize according to missingness
    # this is an intermediate
    M <- !is.na( P )
    # this is the actual M matrix
    M <- if ( loci_on_cols ) tcrossprod( M ) else crossprod( M )
    
    # now we center P
    # original P is no longer needed, let's just overwrite it
    P <- P - 1/2
    # now zero out the NAs
    P[ is.na( P ) ] <- 0

    # now we can calculate A
    A <- if ( loci_on_cols ) tcrossprod( P ) else crossprod( P )
    # complete average and finish the calculation of A
    A <- A / M - 1/4
    
    # transfer names from P to A if present
    # this will carry over all the way to the final coancestry matrix!
    # (M need not have names at all)
    if ( !is.null( names_P ) ) {
        colnames( A ) <- names_P
        rownames( A ) <- names_P
    }
    
    # the denominator is a simple average, a scalar shared by all individuals
    A_min <- popkin_A_min_subpops( A, subpops )
    
    # this is the coancestry estimate
    coancestry <- 1 - A / A_min
    
    # figure out what to return
    if ( want_M ) {
        return(
            list(
                coancestry = coancestry,
                M = M
            )
        )
    } else {
        return( coancestry )
    }
}


