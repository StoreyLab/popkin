# build package documentation, run tests, 
# initially copied from http://kbroman.org/pkg_primer/pages/docs.html

# actually runs document and test internally, but doesn't install
all:
	R -e 'devtools::check()'

doc:
	R -e 'devtools::document()'

test:
	R -e 'devtools::test()'

vig:
	R -e 'devtools::build_vignettes()'

.PHONY: man

man:
	cd ..; if [ -f popkin.pdf ]; then rm popkin.pdf; fi; R CMD Rd2pdf popkin

install:
	R -e 'devtools::install()'
