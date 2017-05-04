# build package documentation, run tests, 
# initially copied from http://kbroman.org/pkg_primer/pages/docs.html

# actually runs document and test internally, but doesn't install
all:
	R -e 'devtools::check()'

doc:
	R -e 'devtools::document()'

test:
	R -e 'devtools::test()'

install:
	R -e 'devtools::install()'
