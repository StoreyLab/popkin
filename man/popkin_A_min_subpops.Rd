% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/popkin_A_min_subpops.R
\name{popkin_A_min_subpops}
\alias{popkin_A_min_subpops}
\title{Estimate the minimum expected value of a matrix \code{A} using subpopulations}
\usage{
popkin_A_min_subpops(A, subpops = NULL)
}
\arguments{
\item{A}{A symmetric \code{n}-by-\code{n} matrix with values between every individual pair, including self comparisons.}

\item{subpops}{A length-\code{n} vector of subpopulation assignments for each individual.
If missing, every individual is treated as a different subpopulation.}
}
\value{
The minimum of the average between-subpopulation \code{A} values, which estimates the minimum expected value of \code{A}
}
\description{
This function averages the values of a square matrix \code{A} between every subpopulation pair and returns the minimum of these averages.
If no subpopulation partition is provided, the function returns the minimum value of \code{A} excluding the diagonal, to agree when the code treats each individual as a subpopulation.
The return value can be used to adjust an \code{A} matrix to yield the kinship matrix.
}
\examples{
# Construct toy data
X <- matrix(c(0,1,2,1,0,1,1,0,2), nrow=3, byrow=TRUE) # genotype matrix
subpops <- c(1,1,2) # subpopulation assignments for individuals

# NOTE: for BED-formatted input, use BEDMatrix!
# "file" is path to BED file (excluding .bed extension)
## library(BEDMatrix)
## X <- BEDMatrix(file) # load genotype matrix object

# calculate A from genotypes
A <- popkin_A(X)$A

# the recommended form using appropriate subpopulation labels
A_min_est <- popkin_A_min_subpops( A, subpops )

# this recovers the popkin estimate
kinship <- 1 - A / A_min_est
stopifnot( kinship == popkin( X, subpops ) )

# a simple default for exploratory analysis, equals min( A ) for correctly-calculated A
A_min_est <- popkin_A_min_subpops( A )
stopifnot( A_min_est == min( A ) )

}
\seealso{
\code{\link[=popkin_A]{popkin_A()}} to generate the \code{A} matrix usually inputted into this function (\code{popkin_A_min_subpops}).
\code{\link[=popkin]{popkin()}} is the wrapper function around both of these.

\code{\link[=avg_kinship_subpops]{avg_kinship_subpops()}} for the full matrix of mean kinship values between subpopulations.
}
