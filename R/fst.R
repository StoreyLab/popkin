#' Calculate FST from a population-level kinship matrix or vector of inbreeding coefficients
#'
#' This function simply returns the weighted mean inbreeding coefficient.
#' If weights are `NULL` (default), the regular mean is calculated.
#' If a kinship matrix is provided, then the inbreeding coefficients are extracted from its diagonal using `\link{inbr}` (requires the diagonal to contains self-kinship values as `\link{popkin}` returns, and not inbreeding coefficients as `\link{inbr_diag}` returns).
#' If there is local inbreeding and it can be estimated (from known pedigrees, for example), it can be subtracted from the total inbreeding coefficients, resulting in a vector of structural inbreeding that correctly averages into FST.
#'
#' The returned weighted mean inbreeding coefficient equals the generalized FST if all individuals are "locally outbred" (i.e. if the self-relatedness of every individual stems entirely from the population structure rather than due partly to having unusually closely related parents, such as first or second cousins).
#' Note most individuals in population-scale human data are locally outbred.
#' If there are locally-inbred individuals, but their local inbreeding cannot be estiamted, then the returned value will overestimate FST.
#' Good estimates of local inbreeding can be passed (parameter `x_local`), in which case the code will subtract their effect and FST will be more accurate.
#' 
#' @param x The vector of inbreeding coefficients, or the kinship matrix if `x` is a matrix.
#' @param weights Weights for individuals (optional, defaults to uniform weights)
#' @param x_local An optional vector of inbreeding coefficients, or a local kinship matrix if `x_local` is a matrix.
#'
#' @return FST
#'
#' @examples
#' # Get FST from a genotype matrix
#' 
#' # Construct toy data
#' X <- matrix(c(0,1,2,1,0,1,1,0,2), nrow = 3, byrow = TRUE) # genotype matrix
#' subpops <- c(1,1,2) # subpopulation assignments for individuals
#' 
#' # NOTE: for BED-formatted input, use BEDMatrix!
#' # "file" is path to BED file (excluding .bed extension)
#' ## library(BEDMatrix)
#' ## X <- BEDMatrix(file) # load genotype matrix object
#'
#' # estimate the kinship matrix "kinship" from the genotypes "X"!
#' kinship <- popkin(X, subpops) # calculate kinship from X and optional subpop labels
#' weights <- weights_subpops(subpops) # can weigh individuals so subpopulations are balanced
#' Fst <- fst(kinship, weights) # use kinship matrix and weights to calculate fst
#' 
#' Fst <- fst(kinship) # no (or NULL) weights implies uniform weights
#'
#' inbr <- inbr(kinship) # if you extracted inbr for some other analysis...
#' Fst <- fst(inbr, weights) # ...use this inbreeding vector as input too!
#'
#' @export
fst <- function(x, weights = NULL, x_local = NULL) {
    # validate inputs
    if (missing(x))
        stop('you must provide a kinship matrix or vector of inbreeding coefficients!')
    
    # if input is a matrix, let's assume it is the kinship matrix, so extract the inbreeding coefficients first
    # let inbr validate kinship matrix
    if (is.matrix(x))
        x <- inbr(x)
    
    # process and apply local data correction, if present
    if ( !is.null( x_local ) ) {
        # get inbreeding if needed
        if ( is.matrix( x_local ) )
            x_local <- inbr( x_local )
        # subtract local inbreeding from total inbreeding to get structural inbreeding vector only!
        x <- ( x - x_local ) / ( 1 - x_local )
    }
    
    # now x is a vector of inbreeding coefficients, adjusted if that was necessary
    # only choice left is how to weigh average
    if (is.null(weights)) {
        return( mean(x) ) # no weights implies uniform weights
    } else {
        if (length(weights) != length(x)) {
            # sanity check: make sure these have matching dimensions!
            # note that x is a vector at this point (if it was a matrix earlier, the vector of inbreeding coefficients has been extracted), so this will always work when the data makes sense
            stop('the number of individuals differs between inbreeding vector (', length(x), ') and weight vector (', length(weights), ')')
        } else {
            return( drop( x %*% weights ) ) # weighted mean
        }
    }
}

