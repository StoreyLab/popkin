% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/popkin-package.R
\docType{package}
\name{popkin-package}
\alias{popkin-package}
\title{A package for estimating kinship and FST under arbitrary population structure}
\description{
The heart of this package is the \code{\link[=popkin]{popkin()}} function, which estimates the kinship matrix of all individual pairs from their genotype matrix.
Inbreeding coefficients, the generalized \code{FST}, and the individual-level pairwise \code{FST} matrix are extracted from the kinship matrix using \code{\link[=inbr]{inbr()}}, \code{\link[=fst]{fst()}}, and \code{\link[=pwfst]{pwfst()}}, respectively.
\code{\link[=fst]{fst()}} accepts weights for individuals to balance subpopulations obtained with \code{\link[=weights_subpops]{weights_subpops()}}.
Kinship matrices can be renormalized (to change the most recent common ancestor population or MRCA) using \code{\link[=rescale_popkin]{rescale_popkin()}}.
Coancestry matrices can be estimated from allele frequency matrices using \code{\link[=popkin_af]{popkin_af()}}.
Lastly, kinship, coancestry, and pairwise FST matrices can be visualized using \code{\link[=plot_popkin]{plot_popkin()}} (with the help of \code{\link[=inbr_diag]{inbr_diag()}} for kinship matrices only).
}
\examples{
# estimate and visualize kinship and FST from a genotype matrix

# Construct toy data
X <- matrix(c(0,1,2,1,0,1,1,0,2), nrow = 3, byrow = TRUE) # genotype matrix
subpops <- c(1,1,2) # subpopulation assignments for individuals
subpops2 <- 1:3 # alternate labels treat every individual as a different subpop

# NOTE: for BED-formatted input, use BEDMatrix!
# "file" is path to BED file (excluding .bed extension)
## library(BEDMatrix)
## X <- BEDMatrix(file) # load genotype matrix object

# estimate the kinship matrix from the genotypes "X"!
# all downstream analysis require "kinship", none use "X" after this
kinship <- popkin(X, subpops) # calculate kinship from X and optional subpop labels

# plot the kinship matrix, marking the subpopulations
# note inbr_diag replaces the diagonal of kinship with inbreeding coefficients
plot_popkin( inbr_diag(kinship), labs = subpops )

# extract inbreeding coefficients from kinship
inbreeding <- inbr(kinship)

# estimate FST
weights <- weights_subpops(subpops) # weigh individuals so subpopulations are balanced
Fst <- fst(kinship, weights) # use kinship matrix and weights to calculate fst
Fst <- fst(inbreeding, weights) # estimate more directly from inbreeding vector (same result)

# estimate and visualize the pairwise FST matrix
pairwise_fst <- pwfst(kinship) # estimated matrix
leg_title <- expression(paste('Pairwise ', F[ST])) # fancy legend label
# NOTE no need for inbr_diag() here!
plot_popkin(pairwise_fst, labs = subpops, leg_title = leg_title)

# rescale the kinship matrix using different subpopulations (implicitly changes the MRCA)
kinship2 <- rescale_popkin(kinship, subpops2)

# toy allele frequency data
P <- matrix( runif( 9 ), nrow = 3 )

# estimate coancestry from allele frequency matrix
# (useful from P matrices from admixture models)
coancestry <- popkin_af( P )

}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/StoreyLab/popkin/}
  \item Report bugs at \url{https://github.com/StoreyLab/popkin/issues}
}

}
\author{
\strong{Maintainer}: Alejandro Ochoa \email{alejandro.ochoa@duke.edu} (\href{https://orcid.org/0000-0003-4928-3403}{ORCID})

Authors:
\itemize{
  \item John D. Storey \email{jstorey@princeton.edu} (\href{https://orcid.org/0000-0001-5992-402X}{ORCID})
}

}
