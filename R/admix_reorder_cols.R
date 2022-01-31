#' Reorder admixture matrix columns
#'
#' Returns the order of the columns (ancestries) of an admixture matrix so that they are in their average order of appearance in rows (individuals).
#' More specifically, for each ancestry it calculates its mean row (expected row number weighted by this ancestry's proportion distribution among rows), and returns the order in which these mean row values are increasing.
#' In datasets where the rows/individuals are already ordered in a meaningful way (for example, by distance from the species' geographical origin, and generally grouping the most similar individuals together), this function can lead to a more pleasing automated visualization of the admixture proportions.
#'
#' @param Q The admixture proportions matrix.
#'
#' @return The desired order of the columns (a vector of indexes).
#'
#' @examples
#' # here is a toy admixture proportions matrix with columns in no meaningful order
#' Q <- matrix(
#'     c(
#'         0.1, 0.8, 0.1,
#'         0.1, 0.7, 0.2,
#'         0.0, 0.4, 0.6,
#'         0.0, 0.3, 0.7,
#'         0.9, 0.0, 0.1
#'     ),
#'     nrow = 5,
#'     ncol = 3,
#'     byrow = TRUE
#' )
#' # get nicer order
#' indexes <- admix_order_cols( Q )
#' # apply reordering to columns
#' Q <- Q[ , indexes ]
#' 
#' # notice that now the first columns takes on the highest values initially,
#' # followed by the second column, and lastly the third column.
#'
#' @seealso
#' [admix_label_cols()] to automatically assign labels to ancestries given labels to individuals.
#' 
#' [plot_admix()] for plotting admixture matrices.
#'
#' @export
admix_order_cols <- function( Q ) {
    if ( missing( Q ) )
        stop( 'Admixture matrix `Q` is required!' )
    n <- nrow( Q )
    weights <- ( ( 1 : n ) %*% Q ) / colSums( Q )
    return( order( weights ) )
}
