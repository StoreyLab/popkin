# build package documentation, run tests, 
# initially copied from http://kbroman.org/pkg_primer/pages/docs.html

# actually runs document and test internally, but doesn't install
# --no-build-vignettes so we use the vignettes that are there already (were manually compressed with vigcomp, which is undone without this option)
all:
	R -e 'devtools::check(vignettes=FALSE)'

build:
	R -e 'devtools::build(vignettes=FALSE)'

buildWinR:
	R -e 'devtools::build_win(args="--no-build-vignettes", version="R-release")'

buildWinD:
	R -e 'devtools::build_win(args="--no-build-vignettes")' # R-devel

doc:
	R -e 'devtools::document()'

test:
	R -e 'devtools::test()'

vig:
	R -e 'devtools::build_vignettes()'

vigcomp:
	R -e 'tools::compactPDF("inst/doc/", gs_quality = "printer")'

revdep:
	R -e 'devtools::revdep_check()'

.PHONY: man

.PHONY: build

man:
	cd ..; if [ -f popkin.pdf ]; then rm popkin.pdf; fi; R CMD Rd2pdf popkin

install:
	R -e 'devtools::install()'

### required steps (first time only):
# dnf install aspell aspell-en valgrind
# install.packages(c('rversions', 'hunspell'))
### can't run this way (must be interactive session), but need this trick:
# release:
# 	R -e 'devtools::release(args="--no-build-vignettes")'

### after building, ran this separately:
# R CMD check --use-valgrind popkin_1.0.3.tar.gz 
