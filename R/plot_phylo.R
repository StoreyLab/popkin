#' Plot a `phylo` tree object
#'
#' This is a wrapper around [ape::plot.phylo()] that makes several adjustments so plots agree more with accompanying kinship matrices.
#' In particular, tree is reversed on the y-axis to match matrix orientation, y-axis spacing is more padded for small trees, and an x-axis scale is always added.
#'
#' @param tree A `phylo` object to plot.
#' @param xlab The x-axis label (default "Coancestry").
#' @param xmax X-axis maximum limit.
#' @param leg_n The desired number of ticks in the x-axis (input to [pretty()], see that for more details).
#' @param ... Additional parameters passed to [ape::plot.phylo()].
#' However, these parameters cannot be passed: `x.lim` (controlled via `xmax`), `y.lim` (a better default for small trees is passed and cannot be changed) and `font` (takes the value of `par('font')` instead of `ape`'s default of 3 (italic)).
#'
#' @examples
#' # create a small random tree
#' library(ape)
#' tree <- rtree( 3 )
#'
#' # plot it!
#' plot_phylo( tree )
#'
#' @seealso
#' [plot_popkin()] can create multipanel figures including kinship matrices and trees (calling the present function in the process).
#'
#' @export
plot_phylo <- function(
                       tree,
                       xlab = 'Coancestry',
                       xmax = NULL,
                       leg_n = 5,
                       ...
                       ) {
    # validate inputs
    if ( missing( tree ) )
        stop( '`tree` is required!' )
    # this one catches an earlier bug (hopefully avoided now in general)
    if ( is.na( leg_n ) )
        stop( '`leg_n` cannot be `NA`!' )

    # determine limits to use (shared precalculated range OR per-tree range (get here))
    if ( is.null( xmax ) )
        xmax <- phylo_max_edge( tree )
    # minimum is always zero
    xlim <- c( 0, xmax )
    
    # reverse tree edges
    # (to order like kinship matrices are ordered, descending instead of ascending on the y axis)
    tree <- reverse_edges_phylo( tree )

    # need number of tips to fix ylim
    n_tips <- ape::Ntip( tree )
    
    # plot the tree
    ape::plot.phylo(
             tree,
             # bold labels, use whatever is set globally (for me this is often bold), instead of the default for plot.phylo of 3 (italic)
             font = graphics::par('font'),
             # default limit fills edges, doesn't look good for small trees
             y.lim = c( 0.5, n_tips + 0.5 ),
             x.lim = xlim,
             ...
         )

    ## # default doesn't add axis and labels, add it now
    ## ape::axisPhylo( backward = FALSE )
    # the default `ape::axisPhylo` doesn't work for this application, let's construct our own axis, especially needed for cases where the range is shared
    # decide where labels will go
    lv <- pretty( xlim, n = leg_n )
    # finally ready to add axis!
    graphics::axis( 1, at = as.numeric( lv ), labels = lv )
    
    if ( !is.na( xlab ) )
        graphics::mtext( xlab, side = 1, line = 2 )

}

# reverse edge order
# this is better than ladderize or other options
# (doesn't reverse additive edges, no need for this currently)
reverse_edges_phylo <- function( tree ) {
    # reverse indexes
    indexes <- nrow( tree$edge ) : 1
    # apply reversion to edge matrix
    tree$edge <- tree$edge[ indexes, ]
    # and also to edge values
    tree$edge.length <- tree$edge.length[ indexes ]

    return( tree )
}
