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

