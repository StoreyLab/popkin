# 2017-09-11 - popkin 1.0.0.9000

* Public release!

# 2017-11-21 - popkin 1.0.1.9000

* Fix a bug in which genotypes input to popkin via a function (rather than a regular matrix or a BEDMatrix object) caused popkin to die.  Now popkin behaves as expected.  New test unit cases were added to test function inputs (previously this case was untested).

# 2017-11-24 - popkin 1.0.2.9000

* Added option to set colors for the lines that separate subpopulations.

# 2018-01-08 - popkin 1.0.3

* Minor non-code changes for first CRAN submission.

# 2018-01-13 - popkin 1.0.4

* All doc examples are now run (all used to be "dontrun").

* Other minor non-code changes for first CRAN submission.

# 2018-01-26 - popkin 1.0.5

* Updated vignette code to work when suggested package "lfa" is not available (needed for CRAN tests).  This change is not visible in rendered vignette included in package.

# 2018-02-01 - popkin 1.0.5.9000

* plotPopkin now allows NULL elements in input list x, makes empty plots with titles (good for placeholders or other non-existent data)

* Clarified plotPopkin documentation (that marPad is added to xMar values if set)

* README.md now contains instructions for installing from CRAN as well as from GitHub.

# 2018-07-30 - popkin 1.0.6.9000

* Internal function printLabs (used by plotPopkin) is now more flexible in where it places its labels (new args "side1" and "side2")

# 2018-08-08 - popkin 1.0.7.9000

* Added option for continuous colors, off by default.  Default is still to use only the 17 colors given directly by RColorBrewer.

# 2018-09-05 - popkin 1.0.8.9000

* Changed some function parameter defaults from missing to NULL, added more validation tests (affects fst, inbr, plotPopkin).

# 2018-10-19 - popkin 1.1.0.9000

* Added neff function (estimates effective sample size given a kinship matrix and weights; can find optimal weights that are non-negative or sign-unconstrained, yielding maximum neff values)

# 2019-02-13 - popkin 1.1.1.9000

* Now the `popkin` function preserves the individual names if they are present in the input genotype matrix.
These names get copied to the rows and columns of the output kinship matrix.

* Converted the vignette from PDF to HTML

# 2019-02-13 - popkin 1.1.2

* Minor non-code changes for second CRAN submission.

# 2019-04-09 - popkin 1.9.0.9000

* Essentially beta version of 2.0.x series
* Renamed functions to fit tidyverse naming style
  * `inbrDiag` -> `inbr_diag`
* `inbr_diag` now accepts lists of kinship matrices to transform (for easier plotting of multiple matrices)
* Added more input checks to functions, informative error messages
* Added functions: `validate_kinship`
