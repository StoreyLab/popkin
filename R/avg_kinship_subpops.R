#' Calculate a kinship matrix between subpopulations by averaging individual data
#'
#' This function calculates a kinship matrix between subpopulations, whose values are the average kinship values between all individual pairs where one individual is in the first subpopulation and the other individual is in the second subpopulation.
#' To estimate coanestry instead of kinship, which is recommended to get more interpretable diagonal values, the input kinship matrix should be transformed using [inbr_diag()].
#'
#' @param kinship A symmetric `n`-by-`n` kinship matrix.
#' @param subpops The length-`n` vector of subpopulation assignments for each individual.
#' @param subpop_order The optional order of subpopulations in the output matrix.
#' `subpop_order` must contain every unique subpopulation in `subpops`.
#' Any additional subpopulations in `subpop_order` (missing in `subpops`) are ignored.
#' By default, subpopulations are in the order of first appearance in `subpops`.
#'
#' @return The symmetric `K`-by-`K` kinship matrix between subpopulations, where `K` is the number of unique subpopulations in `subpops`, ordered as in `subpop_order`.
#'
#' @examples
#' # a toy kinship matrix with 5 individuals belonging to 2 subpopulations
#' kinship <- matrix(
#'     c(
#'         0.7, 0.4, 0.4, 0.1, 0.0,
#'         0.4, 0.7, 0.4, 0.2, 0.1,
#'         0.4, 0.4, 0.7, 0.2, 0.0,
#'         0.1, 0.2, 0.2, 0.6, 0.1,
#'         0.0, 0.1, 0.0, 0.1, 0.6
#'     ),
#'     nrow = 5,
#'     ncol = 5
#' )
#' subpops <- c(1, 1, 1, 2, 2)
#' 
#' # calculate mean kinship between (and within) subpopulations
#' # a 2x2 matrix
#' avg_kinship_subpops( kinship, subpops )
#' 
#' # calculate coancestry estimate instead (difference is diagonal)
#' avg_kinship_subpops( inbr_diag( kinship ), subpops )
#' 
#' @export
avg_kinship_subpops <- function(
                                kinship,
                                subpops,
                                subpop_order = unique( subpops )
                                ) {
    # check inputs
    if ( missing( kinship ) )
        stop( '`kinship` is required!' )
    if ( missing( subpops ) )
        stop( '`subpops` is required!' )

    # check the dimensions, etc
    # NOTE: `kinship` must be symmetric, otherwise below loop doesn't work
    validate_kinship( kinship )
    
    # otherwise check these dimensions too
    if ( length( subpops ) != nrow( kinship ) )
        stop( 'Number of individuals in `subpops` (', length( subpops ), ') and `kinship` (', nrow( kinship ), ') disagree!' )

    # check if user provided `subpop_order` (if it is not missing)
    if ( !missing( subpop_order ) ) {
        # get unique subset
        subpops_unique <- unique( subpops )
        # missing subpops, if any
        subpops_missing <- subpops_unique[ !( subpops_unique %in% subpop_order) ]
        
        # make sure subpop_order is complete
        if ( length( subpops_missing ) > 0 )
            stop( 'These subpopulations are missing in `subpop_order`:', toString( subpops_missing ) )
        
        # also clean up list: remove values in subpop_order that are missing in subpops
        # (that way we don't have to do it outside, we ensure agreement more generally)
        subpop_order <- subpop_order[ subpop_order %in% subpops ]
    }
    
    # initiate output matrix
    K <- length( subpop_order )
    # we must have at least two subpopulations or the code below fails!
    if ( K < 2 )
        stop( 'At least two subpopulations are required (K = ', K, ')!' )
    kinship_subpops <- matrix( 0, nrow = K, ncol = K )
    
    # start averaging
    for ( i in 1 : K ) {
        indexes_i <- subpops == subpop_order[ i ] # booleans for individuals that are of this population
        # note `i = j` case is included 
        for ( j in 1 : i ) {
            indexes_j <- subpops == subpop_order[ j ] # booleans for individuals that are of this population
            kinship_ij <- mean( kinship[ indexes_i, indexes_j ] ) # the mean value we want
            kinship_subpops[ i, j ] <- kinship_ij # store both ways
            kinship_subpops[ j, i ] <- kinship_ij # store both ways
        }
    }
    
    # store names of subpopulations on matrix
    colnames( kinship_subpops ) <- subpop_order
    rownames( kinship_subpops ) <- subpop_order
    
    return( kinship_subpops )
}
