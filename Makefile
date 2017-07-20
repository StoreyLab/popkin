# build package documentation, run tests, 
# initially copied from http://kbroman.org/pkg_primer/pages/docs.html

# actually runs document and test internally, but doesn't install
# --no-build-vignettes so we use the vignettes that are there already (were manually compressed with vigcomp, which is undone without this option)
all:
	R -e 'devtools::check(build_args="--no-build-vignettes")'

doc:
	R -e 'devtools::document()'

test:
	R -e 'devtools::test()'

vig:
	R -e 'devtools::build_vignettes()'

vigcomp:
	R -e 'tools::compactPDF("inst/doc/", gs_quality = "printer")'

.PHONY: man

man:
	cd ..; if [ -f popkin.pdf ]; then rm popkin.pdf; fi; R CMD Rd2pdf popkin

install:
	R -e 'devtools::install()'
