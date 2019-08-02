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

# 2019-04-10 - popkin 1.2.0.9000

* Renamed functions to fit tidyverse naming style:
  * `inbrDiag` -> `inbr_diag`
  * `neff` -> `n_eff`
  * `plotPopkin` -> `plot_popkin`
  * `rescalePopkin` -> `rescale_popkin`
  * `weightsSubpops` -> `weights_subpops`
  * Several argument names were also updated to be more descriptive (particularly in `plot_popkin`).
  * Functions with old names remain for now as deprecated functions (to be removed in the future); only `plotPopkin` retains the older argument names.
* `inbr_diag` now accepts lists of kinship matrices to transform (for easier plotting of multiple matrices).
* `plot_popkin` now requires its non-NULL inputs to be proper kinship matrices.
  Previously, the code used to somewhat allow for non-square matrices to be visualized, but this case had no guarantees to work.
  The code is cleaner under the assumption of symmetric square matrices.
* Added more input checks to functions, informative error messages.
* Added functions: `validate_kinship`, `mean_kinship`

# 2019-04-24 - popkin 1.2.1.9000

`plot_popkin` bug fixes and enhancements!

* `plot_popkin` now resets graphical parameters when done and after every panel as needed.
  * Fixed a bug where panel margins were not reset per panel.
    In particular, after setting custom margins for one panel, but `NULL` (default) for subsequent panels, the original margins were not reset (instead, the last values were incorrectly propagated).
  * The entire layout (all original `par` values) is now reset after plotting is complete.
  * Updated documentation to reflect new behavior.
* Can now have letters on single-panel plots as long as a single letter is passed to `plot_popkin` option `panel_letters` (default is A-Z, so the default remains to not show letters for a single panel).
* Added `leg_cex` option to `plot_popkin`.

# 2019-05-13 - popkin 1.2.2

* Third CRAN submission.
* Added ORCIDs to authors.
* Added back to `popkin` function the deprecated parameter names `lociOnCols` and `memLim` alongside the new names, to prevent breaking existing code (generate warnings).

# 2019-05-30 - popkin 1.2.3

* `inbr_diag` now handles `NULL` inputs correctly (preserves them as `NULL` without throwing errors).
* `plot_popkin` has a new logical option `null_panel_data`, to change behavior in the presence of `NULL` kinship matrices (whether they must or must not have titles and other parameters).
  * The new default is to not specify any data for `NULL` panels

# 2019-06-05 - popkin 1.2.4

* Non-code changes:
  * Edited .Rbuildignore to stop ignoring README; also removed non-existent files from list
  * Removed unused popkin.Rproj file

# 2019-07-24 - popkin 1.2.4.9000

* Added internal function `solve_m_mem_lim`, which generalizes previous behavior to estimate chunk sizes (in number of loci) given a limited memory and number of individuals for various numbers of matrices (of dimensions (m,n) or (n,n)) and vectors (lengths m or n).
  This function is shared with related projects (such as `popkinsuppl` on GitHub).

# 2019-07-29 - popkin 1.2.5.9000

* Now `solve_m_mem_lim` always returns integer chunk sizes (number of loci).
  Previously the function returned non-integers only if the total matrix size `m` was not provided.

# 2019-08-02 - popkin 1.2.6.9000

* Reorganized internal code, mostly to facilitate use of the internal function `solve_m_mem_lim` in other dependent packages.
  In particular, the internal function `get_mem_lim_m` was removed.
* The `popkin` function accepts the new parameter `mem_factor`.
