## Test environments
* local x86_64-redhat-linux-gnu install, R 3.4.3
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* Possibly mis-spelled words in DESCRIPTION:
    FST (2:29, 8:171)
    biallelic (8:99)

  These two words are spelled correctly.

## Downstream dependencies
Tested 'bnpsd' (only downstream dependency) and found no errors or warnings.

## Comments from last submission

* You are not using Suggests packages conditionally as required in §1.1.3.1 of the manual.
Output from https://cran.r-project.org/web/checks/check_results_popkin.html :
Version: 1.0.4
Check: re-building of vignette outputs
Result: WARN
    Error in re-building vignettes:
     ...
    Quitting from lines 94-98 (popkin.Rmd)
    Error: processing vignette 'popkin.Rmd' failed with diagnostics:
    there is no package called 'lfa'
    Execution halted
Flavors: r-devel-linux-x86_64-fedora-clang, r-devel-linux-x86_64-fedora-gcc

  I have modified the vignette to rebuild without errors when 'lfa' is not available.
