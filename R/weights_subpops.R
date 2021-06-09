#' Get weights for individuals that balance subpopulations
#'
#' This function returns positive weights that sum to one for individuals using subpopulation labels, such that every subpopulation receives equal weight.
#' In particular, if there are `K` subpopulations, then the sum of weights for every individuals of a given subpopulation will equal `1 / K`.
#' The weight of every individual is thus inversely proportional to the number of individuals in its subpopulation.
#' If the optional sub-subpopulation labels are also provided, then each sub-subpopulation within a given subpopulation is also weighted equally.
#'
#' @param subpops The length-`n` vector of subpopulation assignments for each individual.
#' @param subsubpops The optional length-`n` vector of sub-subpopulation assignments for each individual.
#' Each sub-subpopulation must belong to a single subpopulation (a nested hierarchy) or an error is produced.
#'
#' @return The length-`n` vector of weights for each individual.
#'
#' @examples
#' # if every individual has a different subpopulation, weights are uniform:
#' subpops <- 1:10
#' weights <- weights_subpops( subpops )
#' stopifnot( all( weights == rep.int( 1/10, 10 ) ) )
#'
#' # subpopulations can be strings too
#' subpops <- c('a', 'b', 'c')
#' weights <- weights_subpops( subpops )
#' stopifnot( all( weights == rep.int( 1/3, 3 ) ) )
#' 
#' # if there are two subpopulations
#' # and the first has twice as many individuals as the second
#' # then the individuals in this first subpopulation weight half as much 
#' # as the ones in the second subpopulation
#' subpops <- c(1, 1, 2)
#' weights <- weights_subpops( subpops )
#' stopifnot( all( weights == c( 1/4, 1/4, 1/2 ) ) )
#'
#' # hierarchy example
#' subpops <- c(1, 1, 1, 2, 2)
#' subsubpops <- c('a', 'b', 'b', 'c', 'd')
#' weights <- weights_subpops( subpops, subsubpops )
#' stopifnot( all( weights == c( 1/4, 1/8, 1/8, 1/4, 1/4 ) ) )
#' 
#' @export
weights_subpops <- function(subpops, subsubpops = NULL) {
    # validate inputs
    if ( missing( subpops ) )
        stop('`subpops` is required!')
    if ( is.null( subpops ) )
        stop('`subpops` cannot be NULL!')
    if ( anyNA( subpops ) )
        stop('`subpops` cannot contain NAs!')
    
    if ( is.null( subsubpops ) ) {
        # count number of individuals in each subpopulation
        subpop_counts <- table( subpops )
        # count number of subpopulations
        K <- length( subpop_counts )
        # construct weights this way
        # unusual query `subpop_counts[ subpops ]` actually works with `subpop_counts` of class `table`
        weights <- 1 / ( K * subpop_counts[ subpops ] )
    } else {
        # need to have more complicated weights
        # validate second input
        if ( anyNA( subsubpops ) )
            stop('`subsubpops` cannot contain NAs!')

        # count number of individuals in each subsubpopulation
        subsubpop_counts <- table( subsubpops )

        # NOTE: while `K` and `subpop_counts` below look deceivingly like they did in the `subsubpops = NULL` case, here they are calculated differenty by counting sub-subpopulations rather than counting individuals!

        # now need to aggregate subsubpops into subpops
        # here rows are individuals, which repeat for all in the same subsubpops
        tab <- data.frame( subpops = subpops, subsubpops = subsubpops )
        # now each row is a subsubpop
        tab <- unique( tab )
        # get table at this level, number of subsubpopulations in each subpopulation
        subpop_counts <- table( tab$subpops )
        # count number of subpopulations
        K <- length( subpop_counts )
        
        # check that every subsubpop actually belongs to a single subpop
        # a discrepancy for this number reveals the problem
        if ( sum( subpop_counts ) != length( subsubpop_counts ) )
            stop( 'The number of unique sub-subpopulations (', length( subsubpop_counts ), ') does not match the total number of sub-subpopulations summed over subpopulations (', sum( subpop_counts ), '), which implies that there is at least one sub-subpopulation with membership in more than one subpopulation!  This is not allowed since procedure assumes sub-subpopulations are nested within subpopulations!' )

        # now we're good with weights!
        weights <- 1 / ( K * subpop_counts[ subpops ] * subsubpop_counts[ subsubpops ] )
    }
    # done, return!
    return( as.numeric( weights ) )
}
