# build package documentation
# copied from http://kbroman.org/pkg_primer/pages/docs.html
all:
	R -e 'devtools::document()'
