# popkin <img src="man/figures/logo.png" alt="popkin" align="right" />

The `popkin` ("population kinship") R package estimates the kinship matrix of individuals and FST from their biallelic genotypes.
Our estimation framework is the first to be practically unbiased under arbitrary population structures.

## Installation

The stable version of the package is now on CRAN and can be installed using
```R
install.packages("popkin")
```

The current development version can be installed from the GitHub repository using `devtools`:
```R
install.packages("devtools") # if needed
library(devtools)
install_github('StoreyLab/popkin', build_opts = c())
```

You can see the package vignette, which has more detailed documentation, by typing this into your R session:
```R
vignette('popkin')
```


## Examples

### Input data

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

### `popkin` functions

This is a quick overview of every `popkin` function, covering estimation and visualization of kinship and FST from a genotype matrix.

First estimate the `kinship` matrix from the genotypes `X`.
All downstream analysis require `kinship`, none use `X` after this
```R
library(popkin)
kinship <- popkin(X, subpops) # calculate kinship from X and optional subpop labels
```

Plot the kinship matrix, marking the subpopulations.
Note `inbr_diag` replaces the diagonal of `kinship` with inbreeding coefficients
```R
plot_popkin( inbr_diag(kinship), labs = subpops )
```

Extract inbreeding coefficients from `kinship`
```R
inbreeding <- inbr(kinship)
```

Estimate FST
```R
weights <- weights_subpops(subpops) # weigh individuals so subpopulations are balanced
Fst <- fst(kinship, weights) # use kinship matrix and weights to calculate fst
Fst <- fst(inbreeding, weights) # estimate more directly from inbreeding vector (same result)
```

Estimate and visualize the pairwise FST matrix
```R
pairwise_fst <- pwfst(kinship) # estimated matrix
leg_title <- expression(paste('Pairwise ', F[ST])) # fancy legend label
plot_popkin(pairwise_fst, labs = subpops, leg_title = leg_title) # NOTE no need for inbr_diag() here!
```

Rescale the kinship matrix using different subpopulations (implicitly changes the most recent common ancestor population used as reference)
```R
kinship2 <- rescale_popkin(kinship, subpops2)
```

Estimate the coancestry matrix from a matrix of allele frequencies `P` (useful when `P` comes from an admixture inference model)
```R
coancestry <- popkin_af( P )
```

Please see the `popkin` R vignette for a description of the key parameters and more detailed examples, including complex plots with multiple kinship matrices and multi-level subpopulation labeling.


## Citations

Alejandro Ochoa, John D Storey.  2021.  "Estimating FST and kinship for arbitrary population structures." PLoS Genet 17(1): e1009241. PubMed ID 33465078. [doi:10.1371/journal.pgen.1009241](https://doi.org/10.1371/journal.pgen.1009241). bioRxiv [doi:10.1101/083923](https://doi.org/10.1101/083923) 2016-10-27.

Alejandro Ochoa, John D Storey.  2019.  "New kinship and FST estimates reveal higher levels of differentiation in the global human population." bioRxiv [doi:10.1101/653279](https://doi.org/10.1101/653279).

Alejandro Ochoa, John D Storey.  2016.  "FST And Kinship for Arbitrary Population Structures I: Generalized Definitions." bioRxiv [doi:10.1101/083915](https://doi.org/10.1101/083915).

