# popkin 1.0.0.9000 (2017-09-11)

* Public release!

# popkin 1.0.1.9000 (2017-11-21)

* Fix a bug in which genotypes input to popkin via a function (rather than a regular matrix or a `BEDMatrix` object) caused popkin to die.  Now popkin behaves as expected.  New test unit cases were added to test function inputs (previously this case was untested).

# popkin 1.0.2.9000 (2017-11-24)

* Added option to set colors for the lines that separate subpopulations.

# popkin 1.0.3 (2018-01-08)

* Minor non-code changes for first CRAN submission.

# popkin 1.0.4 (2018-01-13)

* All doc examples are now run (all used to be `dontrun`).

* Other minor non-code changes for first CRAN submission.

# popkin 1.0.5 (2018-01-26)

* Updated vignette code to work when suggested package `lfa` is not available (needed for CRAN tests).  This change is not visible in rendered vignette included in package.

# popkin 1.0.5.9000 (2018-02-01)

* `plotPopkin` now allows NULL elements in input list x, makes empty plots with titles (good for placeholders or other non-existent data)

* Clarified `plotPopkin` documentation (that `marPad` is added to `xMar` values if set)

* `README.md` now contains instructions for installing from CRAN as well as from GitHub.

# popkin 1.0.6.9000 (2018-07-30)

* Internal function `printLabs` (used by `plotPopkin`) is now more flexible in where it places its labels (new args `side1` and `side2`)

# popkin 1.0.7.9000 (2018-08-08)

* Added option for continuous colors, off by default.  Default is still to use only the 17 colors given directly by `RColorBrewer`.

# popkin 1.0.8.9000 (2018-09-05)

* Changed some function parameter defaults from missing to NULL, added more validation tests (affects `fst`, `inbr`, `plotPopkin`).

# popkin 1.1.0.9000 (2018-10-19)

* Added `neff` function (estimates effective sample size given a kinship matrix and weights; can find optimal weights that are non-negative or sign-unconstrained, yielding maximum neff values)

# popkin 1.1.1.9000 (2019-02-13)

* Now the `popkin` function preserves the individual names if they are present in the input genotype matrix.
These names get copied to the rows and columns of the output kinship matrix.

* Converted the vignette from PDF to HTML

# popkin 1.1.2 (2019-02-13)

* Minor non-code changes for second CRAN submission.

# popkin 1.2.0.9000 (2019-04-10)

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

# popkin 1.2.1.9000 (2019-04-24)

`plot_popkin` bug fixes and enhancements!

* `plot_popkin` now resets graphical parameters when done and after every panel as needed.
  * Fixed a bug where panel margins were not reset per panel.
    In particular, after setting custom margins for one panel, but `NULL` (default) for subsequent panels, the original margins were not reset (instead, the last values were incorrectly propagated).
  * The entire layout (all original `par` values) is now reset after plotting is complete.
  * Updated documentation to reflect new behavior.
* Can now have letters on single-panel plots as long as a single letter is passed to `plot_popkin` option `panel_letters` (default is A-Z, so the default remains to not show letters for a single panel).
* Added `leg_cex` option to `plot_popkin`.

# popkin 1.2.2 (2019-05-13)

* Third CRAN submission.
* Added ORCIDs to authors.
* Added back to `popkin` function the deprecated parameter names `lociOnCols` and `memLim` alongside the new names, to prevent breaking existing code (generate warnings).

# popkin 1.2.3 (2019-05-30)

* `inbr_diag` now handles `NULL` inputs correctly (preserves them as `NULL` without throwing errors).
* `plot_popkin` has a new logical option `null_panel_data`, to change behavior in the presence of `NULL` kinship matrices (whether they must or must not have titles and other parameters).
  * The new default is to not specify any data for `NULL` panels

# popkin 1.2.4 (2019-06-05)

* Non-code changes:
  * Edited `.Rbuildignore` to stop ignoring `README`; also removed non-existent files from list
  * Removed unused `popkin.Rproj` file

# popkin 1.2.4.9000 (2019-07-24)

* Added internal function `solve_m_mem_lim`, which generalizes previous behavior to estimate chunk sizes (in number of loci) given a limited memory and number of individuals for various numbers of matrices (of dimensions (m,n) or (n,n)) and vectors (lengths m or n).
  This function is shared with related projects (such as `popkinsuppl` on GitHub).

# popkin 1.2.5.9000 (2019-07-29)

* Now `solve_m_mem_lim` always returns integer chunk sizes (number of loci).
  Previously the function returned non-integers only if the total matrix size `m` was not provided.

# popkin 1.2.6.9000 (2019-08-02)

* Reorganized internal code, mostly to facilitate use of the internal function `solve_m_mem_lim` in other dependent packages.
  In particular, the internal function `get_mem_lim_m` was removed.
* The `popkin` function accepts the new parameter `mem_factor`.

# popkin 1.2.7.9000 (2019-08-08)

`plot_popkin` updates:

* Bug fixes
  * Fixed a minor bug in which labels under `labs_even = TRUE` were not placed correctly.
    The error was most evident for very small samples (i.e. `n = 3` individuals), and was imperceptible otherwise (i.e. `n = 100` or more).
  * Fixed a minor bug in which the diagonal line under `diag_line = TRUE` did not extend fully to extremes.
    This error was again most evident for very small samples, and was imperceptible otherwise.
* New features
  * Added `weights` option, to change width of every individual to highlight individuals with more weight.
  * Added `raster` option, equivalent to `useRaster` option in the `image` function used internally.
    If `weights` are not `NULL`, `raster` is forced to `FALSE` (required for `image` to work in this setting).
    So its only use is to set it when `weights` are null, as needed.
* Other changes
  * Changed default plot range to (0, 1), in terms of boundaries.
    It used to be (1, n) for bin centers, but this setup was too awkward for weighted data.
  * Minor documentation clarifications to some of this function's unchanged arguments.

# popkin 1.2.8.9000 (2019-08-21)

Memory control bugfixes

* Fixed a serious bug that could let memory usage explode when analyzing large datasets.
  * The bug prevented memory usage from being controlled correctly when the genotype matrix did not fit entirely in memory.
  * The bug had no effect if the entire genotype matrix of interest was small enough to fit in memory.
  * The bug was introduced in 2019-08-02, popkin 1.2.6.9000, commit 643a276974171d86ea621df3d25e1937a100d09a
* Tweaked popkin-specific memory formula to better control memory when a `BEDMatrix` object is analyzed
  * This improvement is relative to popkin's formula prior to the aforementioned bug (version 1.2.5.9000 and earlier)
* Internal function `solve_m_mem_lim` now returns memory limit from `get_mem_lim` or user, in addition to the chunk size in both number of loci and in expected memory usage.

Other enhancements

* `n_eff` function now ensures output `n_eff` estimates are in the theoretically valid range of [ 1, 2*n ].
  Numerical issues in small and noisy kinship matrix estimates could lead to out-of-bounds estimates, which are now replaced with their closest boundary values.

# popkin 1.2.9.9000 (2019-09-16)

* In function `plot_popkin`, added option `names_las`
* In internal function `plot_popkin_single`:
  * Changed the default value of `kinship_range` to agree with the default of `plot_popkin` when a single kinship matrix is plotted (as a result, default colors now agree in that case too).
  * Its return value `breaks` is now invisible.
  * No changes through `plot_popkin` are visible, differences are only noticeable calling this internal function `plot_popkin_single` directly.

# popkin 1.2.10.9000 (2019-10-15)

Improvements to function `plot_popkin`:

* Added option `leg_per_panel`, which if true allows each kinship panel to have a different scale (each gets its own legend/color key).
* Extended various `leg_*` options to be able to take on different values per panel.
* Added option `leg_width` to control the width of the legend panels.
  Increased the default width of this legend/color key (from 0.1 to 0.3, as a fraction of the width of the kinship panels), which changes the behavior in the original case when this legend is shared across kinship panels.
  Now the full legend fits in the panel, without needing an outer margin to the right.
* Option `leg_mar` behavior changed.
  Now `leg_mar` can be a scalar, which sets the right margin of the legend panel.
  New default is `leg_mar = 3`, again necessary so the label of the legend fits in the panel.
  Previous behaviors of `leg_mar = NULL` and a full margin specification are retained.

# popkin 1.2.11.9000 (2019-10-17)

More improvements to function `plot_popkin`:

* Added option `oma`, which sets outer margins via `par(oma)` but provides additional useful shortcuts and defaults.
  This changes the default behavior of `plot_popkin` by setting the left outer margin to 1.5 (all other values are zero), whereas before `plot_popkin` did not set any outer margins.
  This new default behavior makes the "Individuals" outer label appear automatically in plots (whereas before, simply calling `plot_popkin` without setting outer margins resulted in this outer-margin y-axis label being hidden from view).
* Extended the behavior of `mar` to accept various shortcuts (scalar values set only bottom and left margin, whereas the second value of a vector of length 2 sets the top margin, which is otherwise zero; in these two cases the right margin is zero).
  Default behavior remains to not change existing margins.
* Vignette was simplified given the new handy defaults and shortcuts.

# popkin 1.3.0 (2019-12-17)

* Fourth CRAN submission.
* Removed deprecated function names: `inbrDiag`, `neff`, `plotPopkin`, `rescalePopkin`, `weightsSubpops`.
  * Removed deprecated parameter names on `popkin` function: `lociOnCols`, `memLim`.
* Updated `class` usage now that matrices return a two-element array in R-devel (required by CRAN).
* Added `calc_leg_width_min` internal function, though it is unfinished and unused.
* Moved logo to `man/figures/`
* Minor Roxygen-related updates.

# popkin 1.3.1.9000 (2020-02-18)

* Added `x_local` parameter to function `fst`, which permits estimation of FST when there is known local inbreeding (estimated from a pedigree or IBD blocks).

# popkin 1.3.2.9000 (2020-06-19)

* `validate_kinship` now tests for symmetry in input kinship matrices too.

# popkin 1.3.3.9000 (2020-07-13)

* Internal (non-exported) function `solve_m_mem_lim` now avoids a rare integer overflow caused when input number of individuals `n` was encoded as an integer and was greater than `sqrt(.Machine$integer.max)`, or 46340.95.

# popkin 1.3.4.9000 (2020-07-23)

* `validate_kinship` now has `sym` option that, if `FALSE`, skips symmetry test (defaults to `TRUE`).
* `plot_popkin` has the same `sym` option passed to `validate_kinship`, but here it defaults to `FALSE` (there is no inherent error caused by plotting non-symmetric matrices).

# popkin 1.3.5.9000 (2020-09-24)

* Function `popkin`
  - Added `want_M` option, which if `TRUE` returns a list containing the `kinship` matrix as well as the pairwise complete count matrix `M`.
  - Added `m_chunk_max` option (default 1000), which sets the maximum number of loci to process at the time.
	The new default behavior reduces memory usage a lot, especially on machines with large memory, without sacrificing speed.
	Original version would use a lot of memory just because it was available, which could be inconvenient when trying to run other processes, and did not result in increased speed, so it was unnecessary at best.
  - Minor documentation improvements.

# popkin 1.3.6.9000 (2020-11-20)

* Added exported functions `popkin_A` (used to be unexported `get_A`) and `popkin_A_min_subpops` (used to be unexported `min_mean_subpops`)
  - These are low-level functions providing intermediate calculations used by the main `popkin` function.
  - Provided for researchers to try to improve the `popkin` method
* Function `validate_kinship` added option `name` (default "kinship") for clear error reports when the matrix being tested is not actually a kinship matrix
  - Internally used with `name = "A"` to validate `A` in `popkin_A_min_subpops`.
* Simplified documentation (most functions) by using markdown roxygen and replacing all LaTeX equations with simpler code equations.

# popkin 1.3.7 (2021-02-09)

* 5th CRAN submission.
* Updated paper citations in `DESCRIPTION`, `README.md` and the vignette, to point to the published method in PLoS Genetics, and also a related preprint of human analysis on bioRxiv.
* Minor spellcheck corrections in documentation.

# popkin 1.3.8 (2021-02-10)

* 6th CRAN submission.
* Removed a warning message triggered on MacOS and other systems (except Linux and Windows) when the `popkin` function is run.
  Free memory is not calculated in these systems and defaults to 1GB, which threw a warning since could cause problems if the actual memory available is less.
  However, since free memory is rarely below 1GB on reasonable systems, throwing this warning had become more problematic than it was useful (it interfered with internal unit testing), so I decided to remove the warning.

# popkin 1.3.8.9000 (2021-02-16)

* Documentation updates:
  - Fixed links to functions, in many cases these were broken because of incompatible mixed Rd and markdown syntax (now markdown is used more fully).


# popkin 1.3.9.9000 (2021-03-16)

- Added function `popkin_af`, which is the analog of `popkin` but for allele frequency matrices instead of genotypes, and as a consequence it estimates coancestry instead of kinship.

# popkin 1.3.10.9000 (2021-05-25)

Overall added tree plotting capabilities and more plotting fine control.

- Added function `plot_phylo` for plotting `phylo` trees.
  This is a wrapper around `ape::plot.phylo` that makes several adjustments so plots agree more with accompanying kinship matrices (package `ape` is now a dependency for this feature).
- Function `plot_popkin` had the following updates:
  - Objects of type `phylo` and `function` are now accepted elements in input list `kinship` (first argument).
    If `phylo`, these trees are plotted via `plot_phylo`.
    If `function`, its code is executed without arguments, which is expected to plot a single panel.
  - Added option `ylab_side` to allow placing labels on x-axis (bottom, but also top, and right side) instead of the default y-axis (left side).
  - Added option `leg_column` for placing legend/color key in any column (default last column, which was the only choice before).
  - Added option `panel_letters_adj` for positioning panel letters more finely, farther into the margin.
    Also, previous hardcoded default of `0` (inside x-axis range) was changed now to `-0.1` (just outside the x-axis range in most cases).
  - Matrix row and column names (when `names = TRUE`) are now always plotted entirely, even if overlapping.
    The old behavior (R's default) plotted names in order and skipped overlapping labels (see `?axis`), which looks prettier but was confusing for this plot as it suggested incorrectly that some individuals or subpopulations were not present.
	The solution is unfortunately a hack, to pass `gap.axis = -1` to `axis` (suggested in `?axis`), which hopefully does not break in the future.
- Function `validate_kinship` now has option `logical = TRUE` to return a logical value instead of throwing errors.

# popkin 1.3.11.9000 (2021-06-04)

- Updated citations in `inst/CITATION` (missed last time I updated them in other locations).

# popkin 1.3.12.9000 (2021-06-09)

- Function `weights_subpops` updates:
  - Now accepts a second optional vector `subsubpops` for calculating weights on two levels.
  - Return value is now a numeric vector (used to be of class `table`).

# popkin 1.3.13 (2021-07-27)

- 7th CRAN submission.
- Removed `LazyData: true` from DESCRIPTION (to avoid a new "NOTE" on CRAN).
- Reformatted this `NEWS.md` slightly to improve its automatic parsing.
- Fixed spelling in documentation.

# popkin 1.3.13.9000 (2021-11-02)

- Added function `avg_kinship_subpops`.
- Function `popkin_A_min_subpops`:
  - Now uses `avg_kinship_subpops` internally to perform the bulk of the calculations
  - When `subpops = NULL`, calculation now returns minimum `A` among off-diagonal elements only (excluding diagonal) rather than the overall minimum of `A`.  There's no difference when `A` is calculated from genotypes (diagonal values are much greater than off-diagonal values), but made the change for consistency when it might differ for arbitrary inputs.
- `README` updated GitHub install instructions for building vignettes.

# popkin 1.3.14.9000 (2021-11-05)

- Function `plot_popkin` fixed a bug when `null_panel_data = TRUE` in which titles that went over panels with `NULL` kinship were incorrectly omitted.

# popkin 1.3.15.9000 (2021-12-02)

- Improved estimation of available memory on linux, fixing a bug where popkin incorrectly believes there is not enough memory.
  - Old: retrieved `MemFree` (from `/proc/meminfo`).  This could underestimate available memory when `Buffers` and `Cached` memory are large (these count as available memory!), and in some cases cause this error:
    ```
    Error in solve_m_mem_lim : 
    The resulting `m_chunk` was negative!  This is because either `mat_n_n` or `vec_n` are non-zero and `n` alone is too large for the available memory (even for `m_chunk == 0`).  The solution is to free more memory (ideal) or to reduce `n` if possible.
    ```
  - New: retrieve `MemAvailable` (still from `/proc/meminfo`), which is ideal but is absent in older linux kernels (<3.14), otherwise fallback into retrieving and returning the sum of `MemFree`, `Buffers`, and `Cached`.  Either way available memory is greater than `MemFree` alone and is also more accurate.
  - Under the hood, cleaned parser considerably and check for several trouble scenarios that were previously taken for granted.

# popkin 1.3.16.9000 (2022-01-15)

- Added function `plot_admix` for making admixture/structure plots with most of the same options as `plot_popkin`!
  - Includes plotting examples in vignette.
  - Cleaned up and extended internal function `print_labels_multi`, `print_labels`.

# popkin 1.3.17 (2022-01-26)

- 8th CRAN submission.
- Fixed some typos.

# popkin 1.3.17.9000 (2022-01-31)

- Function `plot_admix` added options `leg_title_line` and `leg_las`, and changed the default of `leg_mar`, to better accommodate numerous long ancestry labels.

# popkin 1.3.18.9000 (2022-01-31)

- Added functions to complement `plot_admix`:
  - `admix_order_cols`: to automatically order ancestries given ordered individuals.
  - `admix_label_cols`: to automatically assign labels to ancestries given labels to individuals.

# popkin 1.3.19.9000 (2022-04-29)

- Functions `popkin`, `popkin_A`:
  - Added option `mean_of_ratios`, default `FALSE` is original estimator, `TRUE` gives a new estimator that upweighs rare variants, which resembles in this way the standard kinship estimator, and which appears to improve performance in association testing.
  - Fixed minor bug that `M` (one of the return values when `want_M = TRUE`) did not inherit individual names from `X` even though `A` and `kinship` did, and similarly all inherit names when `X` is a function (fixed accidentally when replacing `Rcpp` code with pure R).
  - Internally, for original (`mean_of_ratios = FALSE`) replaced `Rcpp` code with pure R version, which results in large speedups, at a cost of higher memory use (despite my best attempts at improving the original `Rcpp` code, the simpler R code is doing something magically fast I don't understand).
    `Rcpp`, `RcppEigen` dependencies have been dropped as a consequence.
- Tests were separated into more files (contexts), had some cleanups

# popkin 1.3.20.9000 (2022-05-13)

- Internal function `print_labels` fixed bug when `even = TRUE` and the minimum `xb_ind` is not zero, which caused the maximum to be off by `xb_ind`.
  - This bug didn't affect the exported functions that use this function (`plot_popkin` or `plot_admix`) because the minimum `xb_ind` was always zero in those cases.
  - Fixed to handle a new application outside this package.

# popkin 1.3.21.9000 (2022-07-13)

* Function `plot_popkin`
  * Added option `ylab_per_panel` to allow single-panel figures to place y-axis label in inner margin (before that case was forced to use outer margin).
  * Added clarifications to existing options `oma` and `layout_add`, as in some cases you may want to turn off both features to avoid unexpected behaviors (though there are cases where turning off one but not the other also makes sense).

# popkin 1.3.22.9000 (2022-08-04)

- Function `plot_phylo` added option `edge_width`, which defaults to 1.
  - This restores the plotting behavior under the dependency package `ape` version 5.5 and prior, where its function `plot.phylo` (which `popkin::plot_phylo` wraps) had its parameter `edge.width` default to 1.
  - In contrast, starting on `ape` version 5.6 (2021-12-20), `edge.width` defaults to `NULL`, with results in setting it to `par('lwd')`, which had undesirable consequences in my use cases and which is why the old default is overridden in popkin.
  - By extension, `plot_popkin`'s old default edge widths for trees of class `phylo` is also restored.

# popkin 1.3.23 (2023-01-06)

- 9th CRAN submission.
- Added `hgdp_subset` sample data, copied from `lfa`.
  - Removed `lfa` dependency, which was only used for this sample data.
    `lfa` has become unreliable in external testing servers, particularly as it is on Bioconductor and sometimes hard to install on R-devel, so its removal simplifies automatic testing considerably.
- Functions `popkin`, `popkin_A`, and `popkin_af`
  - Removed internal code that obtained available memory on Windows systems, which started to fail recently.  Now only Linux retrieves available memory from system, on all other systems a default of 1GB memory is assumed to be available.
  - Updated documentation.
- Fixed some typos, escaped more function/code names to simplify spellcheck in the future.
- Function `plot_popkin` clarified documentation.

# popkin 1.3.23.9000 (2024-10-02)

- Functions `plot_popkin` and `plot_admix` added option `labs_even_line` to allow location of the end of `labs_even` lines to differ from `labs_line`.

# popkin 1.3.24.9000 (2026-01-02)

- Functions `popkin` and `popkin_af` when `want_M = TRUE` now include `A_min` in their return list.

# popkin 1.3.25.9000 (2026-01-27)

- Function `plot_popkin` added options `xlab`, `xlab_adj`, `xlab_line`, and `xlab_side`.
