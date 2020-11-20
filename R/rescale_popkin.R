#' Rescale kinship matrix to set a given kinship value to zero.
#'
#' If you already have a population kinship matrix, and you desire to estimate the kinship matrix in a subset of the individuals, you could do it the slow way (reestimating starting from the genotypes of the subset of individuals) or you can do it the fast way: first subset the kinship matrix to only contain the individuals of interest, then use this function to rescale this kinship matrix so that the minimum kinship is zero.
#' This rescaling is required when subsetting results in a more recent Most Recent Common Ancestor (MRCA) population compared to the original dataset (for example, if the original data had individuals from across the world but the subset only contains individuals from a single continent).
#' 
#' This function rescales the input `kinship` matrix so that the value `min_kinship` in the original kinship matrix becomes zero, using the formula
#' `kinship_rescaled = ( kinship - min_kinship ) / ( 1 - min_kinship )`.
#' This is equivalent to changing the ancestral population of the data.
#' If subpopulation labels `subpops` are provided (recommended), they are used to estimate `min_kinship` using the function `\link[popkin_A_min_subpops]`, which is the recommended way to set the MRCA population correctly.
#' If both `subpops` and `min_kinship` are provided, only `min_kinship` is used.
#' If both `subpops` and `min_kinship` are omitted, the function sets `min_kinship = min( kinship )`.
#'
#' @param kinship An `n`-by-`n` kinship matrix.
#' @param subpops The length-`n` vector of subpopulation assignments for each individual.
#' @param min_kinship A scalar kinship value to define the new zero kinship.
#'
#' @return The rescaled `n`-by-`n` kinship matrix, with the desired level of relatedness set to zero.
#'
#' @examples
#' # Construct toy data
#' X <- matrix(c(0,1,2,1,0,1,1,0,2), nrow=3, byrow=TRUE) # genotype matrix
#' subpops <- c(1,1,2) # subpopulation assignments for individuals
#' subpops2 <- 1:3 # alternate labels treat every individual as a different subpop
#' 
#' # NOTE: for BED-formatted input, use BEDMatrix!
#' # "file" is path to BED file (excluding .bed extension)
#' ## library(BEDMatrix)
#' ## X <- BEDMatrix(file) # load genotype matrix object
#'
#' # suppose we first estimate kinship without subpopulations, which will be more biased
#' kinship <- popkin(X) # calculate kinship from genotypes, WITHOUT subpops
#' # then we visualize this matrix, figure out a reasonable subpopulation partition
#'
#' # now we can adjust the kinship matrix!
#' kinship2 <- rescale_popkin(kinship, subpops)
#' # prev is faster but otherwise equivalent to re-estimating kinship from scratch with subpops:
#' # kinship2 <- popkin(X, subpops) 
#'
#' # can also manually set the level of relatedness min_kinship we want to be zero:
#' min_kinship <- min(kinship) # a naive choice for example
#' kinship2 <- rescale_popkin(kinship, min_kinship = min_kinship)
#'
#' # lastly, omiting both subpops and min_kinship sets the minimum value in kinship to zero
#' kinship3 <- rescale_popkin(kinship2)
#' # equivalent to both of:
#' # kinship3 <- popkin(X)
#' # kinship3 <- rescale_popkin(kinship2, min_kinship = min(kinship))
#'
#' @export
rescale_popkin <- function(kinship, subpops = NULL, min_kinship = NA) {
    # die if this is missing
    if (missing(kinship))
        stop('`kinship` matrix is required!')
    
    # run additional validations
    validate_kinship(kinship)

    if (is.na(min_kinship))
        min_kinship <- popkin_A_min_subpops(kinship, subpops)
    
    # finally, perform a simple IBD rescaling
    kinship <- (kinship - min_kinship)/(1 - min_kinship) # return this matrix!
}
