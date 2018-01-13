## Test environments
* local x86_64-redhat-linux-gnu install, R 3.4.3
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* New submission

* Possibly mis-spelled words in DESCRIPTION:
    FST (2:29, 8:171)
    biallelic (8:99)

  These two words are spelled correctly.

## Downstream dependencies
There are currently no downstream dependencies for this package.

## Comments from last submission

* Package has a VignetteBuilder field but no prebuilt vignette index.

  Fixed.

* Is there some reference about the method you can add in the Description field in the form Authors (year) <doi:.....>?

  Added a reference to the Description field of the DESCRIPTION file.

* Please also explain acronyms like FST in your description.

  Done. (FST is not an acronym but is commonly described as "Wright's fixation index").

* Is there any reason why most of your examples are wreapped in \dontrun{}? Why not \donttest{} or completely unwrapped?

  All examples have been modified to be runnable by making use of toy data that was previously absent.

