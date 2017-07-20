Estimate Kinship and FST under Arbitrary Population Structure with popkin
===

The `popkin` ("population kinship") R package estimates the kinship matrix of individuals and FST from their biallelic genotypes.
Our estimation framework is the first to be practically unbiased under arbitrary population structures.

Installation
===

You can install the package from the GitHub repository using `devtools`:

```R
install.packages("devtools") # if needed
library(devtools)
install_github('StoreyLab/alexviiia/fst/software/popkin')
```

While the repository remains private to Storey Lab members, install [following these instructions](https://github.com/StoreyLab/misc/blob/master/github.md).

Synopsis of commands
===

This is a quick overview for estimating and visualize kinship and FST from a genotype matrix.

We begin assuming the following data are present for `n` individuals and `m` loci:
* The `m`-by-`n` genotype matrix `X`, containing only unphased biallelic variants encoded as 0,1,2 counting a given reference allele per locus.
* The length-`n` vector `subpops` that assigns each individual to a subpopulation.
The second ...

This example assumes input is in BED format and is loaded using BEDMatrix "file" is path to BED file (excluding .bed extension)
```R
library(BEDMatrix)
X <- BEDMatrix(file) # load genotype matrix object
```

Estimate the kinship matrix "Phi" from the genotypes "X"!
All downstream analysis require "Phi", none use "X" after this
```R
Phi <- popkin(X, subpops) # calculate kinship from X and optional subpop labels
```

Plot the kinship matrix, marking the subpopulations
Note inbrDiag replaces the diagonal of Phi with inbreeding coefficients
```R
plotPopkin( inbrDiag(Phi), labels=subpops )
```

Extract inbreeding coefficients from Phi
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
plotPopkin(pwF, labs=subpops, legTitle=legTitle) # NOTE no need for inbrDiag() here!
```

Rescale the kinship matrix using different subpopulations (implicitly changes the MRCA)
```R
Phi2 <- rescalePopkin(Phi, subpops2)
```


More details
===

Please see the [popkin vignette](https://github.com/StoreyLab/alexviiia/blob/master/fst/software/popkin/inst/doc/popkin.pdf) for a description of the key parameters and more detailed examples, including complex plots with multiple kinship matrices and multi-level subpopulation labeling.

Citations
===

Ochoa, Alejandro, and John D. Storey. 2016a. "FST And Kinship for Arbitrary Population Structures I: Generalized Definitions." bioRxiv doi:10.1101/083915. Cold Spring Harbor Labs Journals.

Ochoa, Alejandro, and John D. Storey. 2016b. "FST And Kinship for Arbitrary Population Structures II: Method of Moments Estimators." bioRxiv doi:10.1101/083923. Cold Spring Harbor Labs Journals.
