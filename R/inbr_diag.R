#' Replace kinship diagonal with inbreeding coefficients
#'
#' The usual kinship matrix contains self-kinship values along their diagonal given by `kinship_jj = ( 1 + inbr_j ) / 2`, where `inbr_j` is the inbreeding coefficients of individual `j`.
#' This function returns a modified kinship matrix with diagonal values replaced with `inbr_j` values (off-diagonal `j != k` values stay the same).
#' The resulting matrix is better for visualization, but is often not appropriate for modeling (e.g. in mixed-effects models for association or heritability estimation).
#'
#' @param kinship A kinship matrix with self-kinship values along the diagonal.
#' Can pass multiple kinship matrices contained in a list.
#' If `NULL`, it is returned as-is.
#'
#' @return The modified kinship matrix, with inbreeding coefficients along the diagonal, preserving column and row names.
#' If the input was a list of kinship matrices, the output is the corresponding list of transformed matrices.
#' `NULL` inputs are preserved without causing errors.
#'
#' @examples
#' #########
#' # illustrate the main transformation on a 2x2 kinship matrix:
#' # same inbreeding values for both individuals
#' inbr <- 0.2
#' # corresponding self kinship (diagonal values) for both individuals
#' kinship_self <- (1 + inbr)/2
#' # kinship between the two individuals
#' kinship_between <- 0.1
#' # actual kinship matrix
#' kinship <- matrix(c(kinship_self, kinship_between, kinship_between, kinship_self), nrow=2)
#' # expected output of inbr_diag (replaces self kinship with inbreeding)
#' kinship_inbr_diag_exp <- matrix(c(inbr, kinship_between, kinship_between, inbr), nrow=2)
#' # actual output from this function
#' kinship_inbr_diag_obs <- inbr_diag(kinship)
#' # verify that they match (up to machine precision)
#' stopifnot( all( abs(kinship_inbr_diag_obs - kinship_inbr_diag_exp) < .Machine$double.eps ) )
#'
#' # for a list of matrices, returns list of transformed matrices:
#' inbr_diag( list(kinship, kinship) )
#' 
#' # a list with NULL values also works
#' inbr_diag( list(kinship, NULL, kinship) )
#' 
#' #########
#' # Construct toy data (to more closely resemble real data analysis)
#' X <- matrix(c(0,1,2,1,0,1,1,0,2), nrow=3, byrow=TRUE) # genotype matrix
#' subpops <- c(1,1,2) # subpopulation assignments for individuals
#' 
#' # NOTE: for BED-formatted input, use BEDMatrix!
#' # "file" is path to BED file (excluding .bed extension)
#' ## library(BEDMatrix)
#' ## X <- BEDMatrix(file) # load genotype matrix object
#'
#' # estimate the kinship matrix from the genotypes "X"!
#' kinship <- popkin(X, subpops) # calculate kinship from X and optional subpop labels
#'
#' # lastly, replace diagonal of kinship matrix with inbreeding coefficients
#' kinship_inbr_diag <- inbr_diag(kinship)
#'
#' @seealso
#' The inverse function is given by `\link[bnpsd]{coanc_to_kinship}`.
#' 
#' @export
inbr_diag <- function(kinship) {
    # die if this is missing
    if (missing(kinship))
        stop('`kinship` matrix is required!')
    
    # if input is a list, process each element (recursively calling self)
    # returns a list of the same length in this case
    if ( is.list(kinship) )
        return ( lapply( kinship, inbr_diag ) )

    # we now assume we're in the singleton case now
    # handle NULL case
    if ( is.null(kinship) )
        return(NULL)
    
    # additional validations
    validate_kinship(kinship)

    # ready to actually process!
    # returns same kinship matrix but with inbreeding along diagonal instead of self-kinship
    diag(kinship) <- inbr(kinship) # this is the only transformation needed

    # return edited matrix
    return(kinship)
}
