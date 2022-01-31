#' Label ancestries based on best match to individual labels
#'
#' Returns labels for each ancestry (columns) of an admixture matrix which is the best matching label among the average individual (rows) of each subpopulation.
#' More specifically, each ancestry is associated to the subpopulation label in which its admixture proportion was the highest averaging over all individuals from that subpopulation.
#' If there are two or more ancestries that match to the same label, these are made unique by appending its order of appearance (if the label is "A", then the first column that matches to it is labeled "A1", the next one "A2", etc).
#'
#' @param Q The admixture proportions matrix.
#' @param labs Subpopulation labels for individuals (rows of `Q`).
#'
#' @return The best label assignments for the ancestries (columns of `Q`), made unique by indexes if there are overlaps.
#'
#' @examples
#' # toy admixture matrix with labels for individuals/rows that match well with ancestry/columns
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
#' labs <- c('X', 'X', 'Y', 'Y', 'Z')
#' 
#' # to calculate matches and save as column names, do this:
#' colnames( Q ) <- admix_label_cols( Q, labs )
#'
#' # expected column names: c('Z', 'X', 'Y')
#' 
#' @seealso
#' [admix_order_cols()] to automatically order ancestries given ordered individuals.
#' 
#' [plot_admix()] for plotting admixture matrices.
#'
#' @export
admix_label_cols <- function( Q, labs ) {
    # most tedious part is aggregating rows by labels
    # repeat for each column
    # rownames are unique labs, colnames are empty/missing
    Q_subpops <- apply( Q, 2, admix_mean_subpops_single, labs )
    # now come up with the labelings based on which is max per column
    col_labs <- rownames( Q_subpops )[ apply( Q_subpops, 2, which.max ) ]
    # only potential issue is repeats, let's disambiguate
    # not pretty but will do for now
    col_labs <- make_unique( col_labs )
    return( col_labs )
}

# aggregation procedure for each column separately (internal).
# not sure if this is fastest way to get what we want or not.
admix_mean_subpops_single <- function( q, labs ) {
    # this performs aggregation
    q_subpops <- stats::aggregate( q, by = list( labs = labs ), FUN = mean )
    # flatten to vector but preserve names
    x <- q_subpops$x
    names( x ) <- q_subpops$labs
    return( x )
}

# my version of make.unique that is more visually pleasing in this case (internal)
make_unique <- function( x ) {
    # identify non-unique labels
    x_tab <- table( x )
    x_rep <- names( x_tab[ x_tab > 1 ] )
    # number them without spacers (unlike make.unique), and all get numbered
    for ( x_name in x_rep ) {
        # get number of instances
        n_name <- x_tab[ x_name ]
        # rename as desired!
        x[ x == x_name ] <- paste0( x_name, 1 : n_name )
    }
    return( x )
}
