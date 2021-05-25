#' Validate a kinship matrix
#'
#' Tests that the input is a valid kinship matrix (a numeric, square, and optionally symmetric R matrix).
#' Intended for matrices to plot and for other uses, including biased estimates, so there is flexibility as to what constitutes a valid kinship matrix.
#' Throws errors if the input is not as above.
#' Can instead return `TRUE`/`FALSE` if `logical = TRUE`.
#' 
#' True kinship matrices have values strictly between 0 and 1, and diagonal values strictly between 0.5 and 1.
#' However, estimated matrices may contain values slightly out of range.
#' For greater flexibility, this function does not check for out-of-range values.
#'
#' @param kinship The kinship matrix to validate.
#' @param sym If `TRUE` (default), the matrix is required to be symmetric.  Otherwise this particular test is skipped.
#' @param name Default "kinship".
#' Change to desired variable name for more informative error messages (i.e. "A" when used to validate the `A` matrix inside `popkin_A_min_subpops`).
#' @param logical If `FALSE` (default), function stops with an error message if the input is not a kinship matrix.
#' If `TRUE`, function instead returns `TRUE` if the input passed all tests (appears to be a valid kinship matrix) or `FALSE` otherwise.
#'
#' @return If `logical = FALSE` (default), nothing.
#' If `logical = TRUE`, returns `TRUE` if the input is a valid kinship matrix, `FALSE` otherwise.
#'
#' @examples
#' # this is a valid (positive) example
#' kinship <- matrix(c(0.5, 0, 0, 0.6), nrow=2)
#' # this will run without errors or warnings
#' validate_kinship(kinship)
#'
#' # negative examples
#' 
#' # dies if input is missing
#' try( validate_kinship() )
#' 
#' # and if input is not a matrix
#' try( validate_kinship( 1:5 ) )
#' 
#' # and for non-numeric matrices
#' char_mat <- matrix(c('a', 'b', 'c', 'd'), nrow=2)
#' try( validate_kinship( char_mat ) )
#' 
#' # and non-square matrices
#' non_kinship <- matrix(1:2, nrow=2)
#' try( validate_kinship( non_kinship ) )
#'
#' # and non-symmetric matrices
#' non_kinship <- matrix(1:4, nrow=2)
#' try( validate_kinship( non_kinship ) )
#' # but example passes if we drop symmetry requirement this way
#' validate_kinship( non_kinship, sym = FALSE )
#'
#' # instead of stopping, can get a logical value
#' # this returns FALSE
#' validate_kinship( non_kinship, logical = TRUE )
#'
#' @export
validate_kinship <- function(kinship, sym = TRUE, name = 'kinship', logical = FALSE) {
    # die if this is missing
    # this should die even if `logical = TRUE`!
    if ( missing( kinship ) )
        stop( '`', name, '` matrix is required!' )
    
    # make sure it is an ordinary matrix or equivalent
    if ( !is.matrix( kinship ) )
        if (logical) return (FALSE) else stop( '`', name, '` must be an R matrix!' )
    # make sure it is numeric
    if ( !is.numeric( kinship ) )
        if (logical) return (FALSE) else stop( '`', name, '` must be numeric!' )
    # check dimensions
    m <- nrow( kinship )
    n <- ncol( kinship )
    if (n != m)
        if (logical) return (FALSE) else stop( '`', name, '` must be a square matrix!  (nrow ', m, ' != ncol ', n, ')' )
    # test symmetry
    if ( sym && !isSymmetric( kinship ) )
        if (logical) return (FALSE) else stop( '`', name, '` must be a symmetric matrix!' )
    
    # if everything passed and we wanted a logical, return TRUE now
    # don't return anything otherwise
    if (logical)
        return (TRUE)
}
