Estimate Kinship and FST under Arbitrary Population Structure with popkin
===

The `popkin` ("population kinship") R package estimates the kinship matrix of individuals and FST from their biallelic genotypes.
Our estimation framework is the first to be practically unbiased under arbitrary population structures.

Installation
===

The stable version of the package is now on CRAN and can be installed using
```R
install.packages("popkin")
```

The current development version can be installed from the GitHub repository using `devtools`:
```R
install.packages("devtools") # if needed
library(devtools)
devtools::install_github('StoreyLab/popkin', build_opts=c())
```

You can see the package vignette, which has additional documentation, by typing this into your R session:
``` r
vignette('popkin')
```


Input data
===

The examples below assume the following R data variables are present for `n` individuals and `m` loci:
* The `m`-by-`n` genotype matrix `X`, containing only unphased biallelic variants encoded as 0,1,2 counting a given reference allele per locus.
* The length-`n` vector `subpops` that assigns each individual to a subpopulation.

The `subpops` vector is not required, but its use is recommended to improve estimation of the baseline kinship value treated as zero.

If your data is in BED format, `popkin` will process it efficiently using BEDMatrix.
If `file` is the path to the BED file (excluding .bed extension):
```R
library(BEDMatrix)
X <- BEDMatrix(file) # load genotype matrix object
```

Synopsis of commands
===

This is a quick overview of every `popkin` function, covering estimation and visualization of kinship and FST from a genotype matrix.

First estimate the kinship matrix `Phi` from the genotypes `X`.
All downstream analysis require `Phi`, none use `X` after this
```R
library(popkin)
Phi <- popkin(X, subpops) # calculate kinship from X and optional subpop labels
```

Plot the kinship matrix, marking the subpopulations.
Note `inbr_diag` replaces the diagonal of `Phi` with inbreeding coefficients
```R
plotPopkin( inbr_diag(Phi), labs=subpops )
```

Extract inbreeding coefficients from `Phi`
```R
inbr <- inbr(Phi)
```

Estimate FST
```R
w <- weightsSubpops(subpops) # weigh individuals so subpopulations are balanced
Fst <- fst(Phi, w) # use kinship matrix and weights to calculate fst
Fst <- fst(inbr, w) # estimate more directly from inbreeding vector (same result)
```

Estimate and visualize the pairwise FST matrix
```R
pwF <- pwfst(Phi) # estimated matrix
legTitle <- expression(paste('Pairwise ', F[ST])) # fancy legend label
plotPopkin(pwF, labs=subpops, legTitle=legTitle) # NOTE no need for inbr_diag() here!
```

Rescale the kinship matrix using different subpopulations (implicitly changes the most recent common ancestor population used as reference)
```R
Phi2 <- rescalePopkin(Phi, subpops2)
```


More details
===

Please see the [popkin vignette](https://github.com/StoreyLab/popkin/raw/master/inst/doc/popkin.pdf) for a description of the key parameters and more detailed examples, including complex plots with multiple kinship matrices and multi-level subpopulation labeling.

Citations
===

Ochoa, Alejandro, and John D. Storey. 2016a. "FST And Kinship for Arbitrary Population Structures I: Generalized Definitions." bioRxiv [doi:10.1101/083915](http://doi.org/10.1101/083915).

Ochoa, Alejandro, and John D. Storey. 2016b. "FST And Kinship for Arbitrary Population Structures II: Method of Moments Estimators." bioRxiv [doi:10.1101/083923](http://doi.org/10.1101/083923).
