#' Get weights for individuals that balance subpopulations
#'
#' This function returns positive weights that sum to one for individuals using subpopulation labels, such that every subpopulation receives equal weight.
#' In particular, if there are \eqn{K} subpopulations, then the sum of weights for every individuals of a given subpopulation will equal \eqn{\frac{1}{K}}{1/K}.
#' The weight of every individual is thus inversely proportional to the number of individuals in its subpopulation.
#'
#' @param subpops The length-\eqn{n} vector of subpopulation assignments for each individual.
#'
#' @return The length-\eqn{n} vector of weights for each individual.
#'
#' @examples
#' # if every individual has a different subpopulation, weights are uniform:
#' subpops <- 1:10
#' w <- weightsSubpops(subpops)
#' stopifnot(all(w == rep.int(1/10,10)))
#'
#' # subpopulations can be strings too
#' subpops <- c('a', 'b', 'c')
#' w <- weightsSubpops(subpops)
#' stopifnot(all(w == rep.int(1/3,3)))
#' 
#' # if there are two subpopulations
#' # and the first has twice as many individuals as the second
#' # then the individuals in this first subpopulation weight half as much 
#' # as the ones in the second subpopulation
#' subpops <- c(1, 1, 2)
#' w <- weightsSubpops(subpops)
#' stopifnot(all(w == c(1/4,1/4,1/2)))
#' 
#' @export
weightsSubpops <- function(subpops) {
    if (missing(subpops) || is.null(subpops)) stop('Fatal: subpopulation assignments are missing!')
    if (anyNA(subpops)) stop('Fatal: subpopulations vector contains NAs!')
    ## count number of individuals in each subpopulation
    subpop2c <- table(subpops)
    ## count number of subpopulations
    K <- length(subpop2c)
    ## construct weights, return!
    w <- 1/(K*subpop2c[match(subpops, names(subpop2c))])
}
