#' A package for estimating kinship and FST under arbitrary population structure
#'
#' The heart of this package is the `\link{popkin}` function, which estimates the kinship matrix of all individual pairs from their genotype matrix.
#' Inbreeding coefficients, the generalized `FST`, and the individual-level pairwise `FST` matrix are extracted from the kinship matrix using `\link{inbr}`, `\link{fst}`, and `\link{pwfst}`, respectively.
#' `\link{fst}` accepts weights for individuals to balance subpopulations obtained with `\link{weights_subpops}`.
#' Kinship matrices can be renormalized (to change the most recent common ancestor population or MRCA) using `\link{rescale_popkin}`.
#' Lastly, kinship and pairwise FST matrices can be visualized using `\link{plot_popkin}` (with the help of `\link{inbr_diag}` for kinship matrices only).
#' 
#' @examples
#' # estimate and visualize kinship and FST from a genotype matrix
#'
#' # Construct toy data
#' X <- matrix(c(0,1,2,1,0,1,1,0,2), nrow = 3, byrow = TRUE) # genotype matrix
#' subpops <- c(1,1,2) # subpopulation assignments for individuals
#' subpops2 <- 1:3 # alternate labels treat every individual as a different subpop
#' 
#' # NOTE: for BED-formatted input, use BEDMatrix!
#' # "file" is path to BED file (excluding .bed extension)
#' ## library(BEDMatrix)
#' ## X <- BEDMatrix(file) # load genotype matrix object
#'
#' # estimate the kinship matrix from the genotypes "X"!
#' # all downstream analysis require "kinship", none use "X" after this
#' kinship <- popkin(X, subpops) # calculate kinship from X and optional subpop labels
#'
#' # plot the kinship matrix, marking the subpopulations
#' # note inbr_diag replaces the diagonal of kinship with inbreeding coefficients
#' plot_popkin( inbr_diag(kinship), labs = subpops )
#'
#' # extract inbreeding coefficients from kinship
#' inbreeding <- inbr(kinship)
#' 
#' # estimate FST
#' weights <- weights_subpops(subpops) # weigh individuals so subpopulations are balanced
#' Fst <- fst(kinship, weights) # use kinship matrix and weights to calculate fst
#' Fst <- fst(inbreeding, weights) # estimate more directly from inbreeding vector (same result)
#'
#' # estimate and visualize the pairwise FST matrix
#' pairwise_fst <- pwfst(kinship) # estimated matrix
#' leg_title <- expression(paste('Pairwise ', F[ST])) # fancy legend label
#' # NOTE no need for inbr_diag() here!
#' plot_popkin(pairwise_fst, labs = subpops, leg_title = leg_title)
#'
#' # rescale the kinship matrix using different subpopulations (implicitly changes the MRCA)
#' kinship2 <- rescale_popkin(kinship, subpops2)
#'
#' @docType package
#' @name popkin-package
#' @aliases popkin-package
"_PACKAGE"

