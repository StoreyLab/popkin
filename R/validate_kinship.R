#' Validate a kinship matrix
#'
#' Tests that the input is a valid kinship matrix (a numeric, square, symmetric R matrix).
#' Throws errors if the input is not as above.
#' 
#' True kinship matrices have values strictly between 0 and 1, and diagonal values strictly between 0.5 and 1.
#' However, estimated matrices may contain values slightly out of range.
#' For greater flexibility, this function does not check for out-of-range values.
#'
#' @param kinship The kinship matrix to validate.
#' @param sym If `TRUE` (default), the matrix is required to be symmetric.  Othewise this particular test is skipped.
#'
#' @return Nothing
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
#' @export
validate_kinship <- function(kinship, sym = TRUE) {
    # die if this is missing
    if (missing(kinship))
        stop('`kinship` matrix is required!')
    # make sure it is an ordinary matrix or equivalent
    if (!is.matrix(kinship))
        stop('`kinship` must be an R matrix!')
    # make sure it is numeric
    if (!is.numeric(kinship))
        stop('`kinship` must be numeric!')
    # check dimensions
    m <- nrow(kinship)
    n <- ncol(kinship)
    if (n != m)
        stop('`kinship` must be a square matrix!  (nrow ', m, ' != ncol ', n, ')')
    # test symmetry
    if ( sym && !isSymmetric( kinship ) )
        stop( '`kinship` must be a symmetric matrix!' )
}
