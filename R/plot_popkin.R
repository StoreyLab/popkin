#' Visualize one or more kinship matrices
#'
#' This function plots one or more kinship matrices and a shared legend for the color key.
#' Many options allow for fine control of individual or subpopulation labeling.
#' This code assumes input matrices are symmetric.
#'
#' \code{plot_popkin} plots the input kinship matrices as-is.
#' For best results, a standard kinship matrix (such as the output of \code{\link{popkin}}) should have its diagonal rescaled to contain inbreeding coefficients (\code{\link{inbr_diag}} does this) before \code{plot_popkin} is used.
#'
#' This function permits the labeling of individuals (from row and column names when \code{names = TRUE}) and of subpopulations (passed through \code{labs}).
#' The difference is that the labels passed through \code{labs} are assumed to be shared by many individuals, and lines (or other optional visual aids) are added to demarcate these subgroups.
#'
#' @param kinship A numeric kinship matrix or a list of matrices.
#' Note \code{kinship} may contain \code{NULL} elements (makes blank plots with titles; good for placeholders or non-existent data)
#' @param titles Titles to add to each matrix panel (default is no titles)
#' @param col Colors for heatmap (default is a red-white-blue palette symmetric about zero constructed using RColorBrewer).
#' @param col_n The number of colors to use in the heatmap (applies if \code{col = NULL}).
#' @param mar Margins for each panel (if a list) or for all panels (if a vector).
#' Margins are in \code{c(bottom,left,top,right)} format that \code{\link[graphics]{par}('mar')} expects.
#' Note the padding \code{mar_pad} below is also added to every margin if set.
#' If \code{NULL}, the original margin values are used without change, and are reset for every panel that has a \code{NULL} value.
#' The original margins are also reset after plotting is complete.
#' @param mar_pad Margin padding added to all panels (\code{mar} above and \code{leg_mar} below).
#' Default 0.2.
#' Must be a scalar or a vector of length 4 to match \code{\link[graphics]{par}('mar')}.
#' @param diag_line If \code{TRUE} adds a line along the diagonal (default no line).
#' May also be a vector of booleans to set per panel (lengths must agree).
#' @param panel_letters Vector of strings for labeling panels (default A-Z).
#' No labels are added when there is only one panel, or if \code{panel_letters = NULL}.
#' @param panel_letters_cex Scaling factor of panel letters (default 1.5).
#'
#' AXIS LABEL OPTIONS
#' 
#' @param ylab The y-axis label (default "Individuals").
#' If \code{length(ylab) == 1}, the label is placed in the outer margin (shared across panels);
#' otherwise \code{length(ylab)} must equal the number of panels and each label is placed in the inner margin of the respective panel.
#' @param ylab_adj The value of "adj" passed to \code{\link[graphics]{mtext}}.
#' If \code{length(ylab) == 1}, only the first value is used, otherwise \code{length(ylab_adj)} must equal the number of panels.
#' @param ylab_line The value of "line" passed to \code{\link[graphics]{mtext}}.
#' If \code{length(ylab) == 1}, only the first value is used, otherwise \code{length(ylab_line)} must equal the number of panels.
#' 
#' LAYOUT OPTIONS
#' 
#' @param layout_add If \code{TRUE} (default) then \code{\link[graphics]{layout}} is called internally with appropriate values for the required number of panels for each matrix, the desired number of rows (see \code{layout_rows} below) plus the color key legend.
#' The original layout is reset when plotting is complete and if \code{layout_add = TRUE}.
#' If a non-standard layout or additional panels (beyond those provided by \code{plot_popkin}) are desired, set to FALSE and call \code{\link[graphics]{layout}} yourself beforehand.
#' @param layout_rows Number of rows in layout, used only if \code{layout_add = TRUE}.
#'
#' LEGEND (COLOR KEY) OPTIONS
#' 
#' @param leg_title The name of the variable that the heatmap colors measure (default "Kinship").
#' @param leg_mar Margin vector (in \code{c(bottom,left,top,right)} format that \code{\link[graphics]{par}('mar')} expects) for the legend panel only.
#' If not provided, the margins used in the last panel are preserved with the exception that the left margin is set to zero (plus the value of \code{mar_pad}, see above).
#' @param leg_n The desired number of ticks in the legend y-axis (input to \code{\link{pretty}}, see that for more details).
#'
#' INDIVIDUAL LABEL OPTIONS
#' 
#' @param names If \code{TRUE}, the column and row names are plotted in the heatmap.
#' @param names_cex Scaling factor for the column and row names.
#' @param names_line Line where column and row names are placed.
#'
#' SUBPOPULATION LABEL OPTIONS
#' 
#' @param labs Subpopulation labels for individuals.
#' Use a matrix of labels to show groupings at more than one level (for a hierarchy or otherwise).
#' If input is a vector or a matrix, the same subpopulation labels are shown for every heatmap panel; the input must be a list of such vectors or matrices if the labels vary per panel.
#' @param labs_cex A vector of label scaling factors for each level of labs, or a list of such vectors if labels vary per panel.
#' @param labs_las A vector of label orientations (in format that \code{\link[graphics]{mtext}} expects) for each level of labs, or a list of such vectors if labels vary per panel.
#' @param labs_line A vector of lines where labels are placed (in format that \code{\link[graphics]{mtext}} expects) for each level of labs, or a list of such vectors if labels vary per panel.
#' @param labs_sep A vector of booleans that specify whether lines separating the subpopulations are drawn for each level of labs, or a list of such vectors if labels vary per panel.
#' @param labs_lwd A vector of line widths for the lines that divide subpopulations (if \code{labs_sep = TRUE}) for each level of labs, or a list of such vectors if labels vary per panel.
#' @param labs_col A vector of colors for the lines that divide subpopulations (if \code{labs_sep = TRUE}) for each level of labs, or a list of such vectors if labels vary per panel.
#' @param labs_ticks A vector of booleans that specify whether ticks separating the subpopulations are drawn for each level of labs, or a list of such vectors if labels vary per panel.
#' @param labs_text A vector of booleans that specify whether the subpopulation labels are shown for each level of labs, or a list of such vectors if labels vary per panel.
#' Useful for including separating lines or ticks without text.
#' @param labs_even A vector of booleans that specify whether the subpopulations labels are drawn with equal spacing for each level of labs, or a list of such vectors if labels vary per panel.
#' When \code{TRUE}, lines mapping the equally-spaced labels to the unequally-spaced subsections of the heatmap are also drawn.
#' 
#' @param ... Additional options passed to \code{\link[graphics]{image}}.
#' These are shared across panels
#'
#' @examples
#' # Construct toy data
#' X <- matrix(c(0,1,2,1,0,1,1,0,2), nrow = 3, byrow = TRUE) # genotype matrix
#' subpops <- c(1,1,2) # subpopulation assignments for individuals
#' 
#' # NOTE: for BED-formatted input, use BEDMatrix!
#' # "file" is path to BED file (excluding .bed extension)
#' ## library(BEDMatrix)
#' ## X <- BEDMatrix(file) # load genotype matrix object
#'
#' # estimate the kinship matrix from the genotypes "X"!
#' kinship <- popkin(X, subpops) # calculate kinship from X and optional subpop labels
#'
#' # simple plot of the kinship matrix, marking the subpopulations only
#' # note inbr_diag replaces the diagonal of kinship with inbreeding coefficients
#' # (see vignette for more elaborate examples)
#' plot_popkin( inbr_diag(kinship), labs = subpops )
#'
#' @export
plot_popkin <- function(
                        kinship,
                        titles = NULL,
                        col = NULL,
                        col_n = 100,
                        mar = NULL,
                        mar_pad = 0.2,
                        diag_line = FALSE,
                        panel_letters = toupper(letters),
                        panel_letters_cex = 1.5,
                        ylab = 'Individuals',
                        ylab_adj = NA,
                        ylab_line = 0,
                        layout_add = TRUE,
                        layout_rows = 1, 
                        leg_title = 'Kinship',
                        leg_mar = NULL,
                        leg_n = 5,
                        names = FALSE,
                        names_cex = 1,
                        names_line = NA,
                        labs = NULL,
                        labs_cex = 1,
                        labs_las = 0,
                        labs_line = 0,
                        labs_sep = TRUE,
                        labs_lwd = 1,
                        labs_col = 'black',
                        labs_ticks = FALSE,
                        labs_text = TRUE,
                        labs_even = FALSE,
                       ...) {
    # wrapper around individual panels and color key
    # does not set PDF output, margins, layout, etc
    # assumes kinship is a list of matrices, if not it is internally turned into one

    # NOTE: we'll let plot_popkin_single validate kinship matrices (no validation here)
    if (!is.list(kinship))
        kinship <- list(kinship) # turn into list so rest works
    
    n <- length(kinship) # number of heatmap panels
    if (!is.null(titles)) {
        if (length(titles) != n)
            stop('titles provided are not the same length as data! Data: ', n, ', titles: ', length(titles))
    } else {
        titles <- rep.int('', n) # make blank titles of the same length as data
    }
    # check label lengths
    if (!is.null(labs) && class(labs) == 'list' && length(labs) != n)
        stop('there are ', n, ' panels but ', length(labs), ' label sets!')
    # expand other things that may vary per panel
    names <- rep_check(names, n)
    names_cex <- rep_check(names_cex, n)
    names_line <- rep_check(names_line, n)
    diag_line <- rep_check(diag_line, n)
    # this is for non-scalars per panel, get turned into lists instead
    mar <- rep_check_list(mar, n)
    labs <- rep_check_list(labs, n)
    labs_cex <- rep_check_list(labs_cex, n)
    labs_las <- rep_check_list(labs_las, n)
    labs_line <- rep_check_list(labs_line, n)
    labs_lwd <- rep_check_list(labs_lwd, n)
    labs_col <- rep_check_list(labs_col, n)
    labs_sep <- rep_check_list(labs_sep, n)
    labs_ticks <- rep_check_list(labs_ticks, n)
    labs_text <- rep_check_list(labs_text, n)
    labs_even <- rep_check_list(labs_even, n)
    # ylab behavior is more dynamic!
    if (length(ylab) > 1) {
        ylab <- rep_check(ylab, n) # just makes sure length is n
        ylab_adj <- rep_check(ylab_adj, n)
        ylab_line <- rep_check(ylab_line, n)
    }

    # figure out which data are non-NULL
    indexes_not_null <- !sapply(kinship, is.null)
    if (!any(indexes_not_null))
        stop('every element of list "kinship" is NULL!')
    
    # code needs two versions of the range
    # - range_real is the real range, used in the end so the color key doesn't show values that weren't actually used
    # - range_sym is a symmetric range, used internally to ensure zero is in the exact middle (set to white in the default)
    # get range and construct symmetric range that helps with plotting nice colors with white at zero
    range_real <- range( unlist(lapply(kinship[indexes_not_null], range, na.rm = TRUE)) )  # range of all non-NULL data plotted
    # these next few lines force symmetry for colors (might look better, not sure)
    max_sym <- max(abs(range_real))
    range_sym <- c(-max_sym, max_sym)

    # save entire original setup, to reset in the end
    # no.readonly is since some parameters cannot be changed (trying to set them results in ugly warnings)
    par_orig <- graphics::par( no.readonly = TRUE )
    # save original margins, which may get reset per panel (separately of final reset)
    mar_orig <- graphics::par('mar')
    
    # figure out layout given a requested number of rows
    if (layout_add)
        plot_popkin_layout(n, layout_rows)

    # breaks of all following plots should match!
    breaks <- NULL
    for (i in 1:n) {
        if (!is.null(mar[[i]])) {
            # change margins if necessary!
            graphics::par(mar = mar[[i]] + mar_pad)
        } else {
            # restore original margins otherwise!
            graphics::par(mar = mar_orig)
        }

        breaks_i <- plot_popkin_single(
            kinship[[i]],
            kinship_range = range_sym,
            col = col,
            col_n = col_n,
            names = names[i],
            names_cex = names_cex[i],
            names_line = names_line[i],
            labs = labs[[i]],
            labs_cex = labs_cex[[i]],
            labs_las = labs_las[[i]],
            labs_line = labs_line[[i]],
            labs_lwd = labs_lwd[[i]],
            labs_sep = labs_sep[[i]],
            labs_ticks = labs_ticks[[i]],
            labs_text = labs_text[[i]],
            labs_col = labs_col[[i]],
            labs_even = labs_even[[i]],
            diag_line = diag_line[i],
            main = titles[i],
            ...
        )
        if (!is.null(breaks_i))
            breaks <- breaks_i # don't overwrite for non-data kinship[[i]] cases
        # add ylab for every panel when there is more than one choice, and provided it was non-NA
        # uses inner rather than outer margin (only choice that makes sense)
        if (length(ylab) > 1 && !is.na(ylab[i])) 
            graphics::mtext(ylab[i], side = 2, adj = ylab_adj[i], line = ylab_line[i])
        # add letters only when there is more than one panel
        # if panel_letters is null then don't add
        if (!is.null(panel_letters) && 1 < n && n <= length(panel_letters)) # for humongous N's we run out of letters, just prevent obvious errors...
            panel_letter(panel_letters[i], cex = panel_letters_cex)
        
    }

    if (!is.null(leg_mar)) {
        graphics::par(mar = leg_mar + mar_pad) # change margins if necessary!
    } else {
        # change the current left margin to the padding value
        marTmp <- graphics::par('mar') # last margins
        marTmp[2] <- mar_pad # replace left margin with zero plus pad
        graphics::par(mar = marTmp) # update margins for legend only!
    }
    heatmap_legend(breaks, kinship_range = range_real, label = leg_title, col = col, col_n = col_n, leg_n = leg_n)

    # add margin only once if there was only one, place in outer margin (only choice that makes sense)
    if (length(ylab) == 1)
        graphics::mtext(ylab, side = 2, adj = ylab_adj, outer = TRUE, line = ylab_line)

    # restore original margins otherwise!
    graphics::par(mar = mar_orig)
    
    # restore original setup when done, but only if we created the default layout
    # otherwise the external layout gets reset, which is bad if we were not done adding panels
    if (layout_add)
        graphics::par( par_orig )
}

# stick deprecated function name here

#' @title Visualize one or more kinship matrices
#' @description Visualize one or more kinship matrices
#' @param ... Params to pass to new plot_popkin
#' @return Nothing
#'
#' @name plotPopkin-deprecated
#' @usage plotPopkin(
#'                       x,
#'                       colCont,
#'                       coln = 100,
#'                       xMar = NULL,
#'                       marPad = 0.2,
#'                       diagLine = FALSE,
#'                       panelLetters = toupper(letters),
#'                       panelLetterCex = 1.5,
#'                       ylabAdj = NA,
#'                       ylabLine = 0,
#'                       addLayout = TRUE,
#'                       nr = 1, 
#'                       legTitle = 'Kinship',
#'                       legMar = NULL,
#'                       nPretty = 5,
#'                       showNames = FALSE,
#'                       namesCex = 1,
#'                       namesLine = NA,
#'                       labsCex = 1,
#'                       labsLas = 0,
#'                       labsLine = 0,
#'                       labsSkipLines = FALSE,
#'                       labsLwd = 1,
#'                       labsCol = 'black',
#'                       labsDoTicks = FALSE,
#'                       labsDoText = TRUE,
#'                       labsEven = FALSE,
#'                       ...
#' )
#' @seealso \code{\link{popkin-deprecated}}
#' @keywords internal
NULL

#' @rdname popkin-deprecated
#' @section \code{plotPopkin}:
#' For \code{plotPopkin}, use \code{\link{plot_popkin}}.
#' Several argument names have also changed!
#'
#' @export
plotPopkin <- function(
                       x,
                       colCont,
                       coln = 100,
                       xMar = NULL,
                       marPad = 0.2,
                       diagLine = FALSE,
                       panelLetters = toupper(letters),
                       panelLetterCex = 1.5,
                       ylabAdj = NA,
                       ylabLine = 0,
                       addLayout = TRUE,
                       nr = 1, 
                       legTitle = 'Kinship',
                       legMar = NULL,
                       nPretty = 5,
                       showNames = FALSE,
                       namesCex = 1,
                       namesLine = NA,
                       labsCex = 1,
                       labsLas = 0,
                       labsLine = 0,
                       labsSkipLines = FALSE,
                       labsLwd = 1,
                       labsCol = 'black',
                       labsDoTicks = FALSE,
                       labsDoText = TRUE,
                       labsEven = FALSE,
                       ...
                       ) {
    # mark as deprecated
    .Deprecated('plot_popkin')
    if (!missing(colCont) && colCont == FALSE)
        warning('`colCont` is deprecated (only `colCont = TRUE` is implemented now)!')
    # return as usual, to not break things just yet
    plot_popkin(
        kinship = x,
        col_n = coln,
        mar = xMar,
        mar_pad = marPad,
        diag_line = diagLine,
        panel_letters = panelLetters,
        panel_letters_cex = panelLetterCex,
        ylab_adj = ylabAdj,
        ylab_line = ylabLine,
        layout_add = addLayout,
        layout_rows = nr,
        leg_title = legTitle,
        leg_mar = legMar,
        leg_n = nPretty,
        names = showNames,
        names_cex = namesCex,
        names_line = namesLine,
        labs_cex = labsCex,
        labs_las = labsLas,
        labs_line = labsLine,
        labs_sep = !labsSkipLines, # gets reversed!
        labs_lwd = labsLwd,
        labs_col = labsCol,
        labs_ticks = labsDoTicks,
        labs_text = labsDoText,
        labs_even = labsEven,
        ...
    )
}


plot_popkin_single <- function (
                                kinship = NULL,
                                kinship_range = range(kinship, na.rm = TRUE),
                                col = NULL,
                                col_n = 100,
                                names = FALSE,
                                names_cex = 1,
                                names_line = NA,
                                xlab = "",
                                ylab = "",
                                labs = NULL,
                                labs_cex = 1,
                                labs_las = 0,
                                labs_line = 0,
                                labs_lwd = 1,
                                labs_sep = TRUE,
                                labs_ticks = FALSE,
                                labs_text = TRUE,
                                labs_col = 'black',
                                labs_even = FALSE,
                                diag_line = FALSE,
                                ...
                                ) {
    # this "raw" version does not plot legend or set margins, best for optimized scenarios...

    # shortcut when there's no data to plot (placeholders)
    # let's still start a new plot with a title, nothing else gets added
    if (is.null(kinship)) {
        graphics::plot.new() # new blank figure
        graphics::title(...) # this works?
        return(NULL) # return null breaks in this case only!
    } else {
        # non-null values must be proper kinship matrices!
        validate_kinship(kinship)
    }
    
    # further data validation
    if (is.null(col))
        col <- plot_popkin_palette(n = col_n) # default coloring
    if (length(col) == 1 && is.character(col)) 
        col <- get(col, mode = "function")
    
    # figure out breaks for colors
    breaks <- seq(kinship_range[1], kinship_range[2], length = length(col) + 1)
    numcols <- length(breaks) - 1
    if (class(col) == "function") 
        col <- col(numcols)
    # here replace potentially forced kinship_range with the true range of the data...
    # but only if data is more extended, don't change anything otherwise!
    rangeRaw <- range(kinship, na.rm = TRUE)
    if ( breaks[1] > rangeRaw[1] )
        breaks[1] <- rangeRaw[1]
    if ( breaks[numcols+1] < rangeRaw[2] )
        breaks[numcols+1] <- rangeRaw[2]
    
    # main plot!
    n <- ncol(kinship)
    # need to reverse the rows of the kinship for plotting!
    # (keep it this way for names below)
    kinship <- kinship[n:1,]
    graphics::image(
                  1:n,
                  1:n,
                  t(kinship),
                  xlim = 0.5 + c(0, n),
                  ylim = 0.5 + c(0, n),
                  axes = FALSE,
                  xlab = xlab,
                  ylab = ylab,
                  col = col,
                  breaks = breaks,
                  useRaster = TRUE,
                  ...
              )
    if (names) {
        # if we want to show labels, use the row/col names to add to picture
        # labels will be perpendicular to axis (las=2)
        graphics::axis(1, 1:n, colnames(kinship), las = 2, cex.axis = names_cex, tick = FALSE, line = names_line)
        graphics::axis(2, 1:n, rownames(kinship), las = 2, cex.axis = names_cex, tick = FALSE, line = names_line)
    }
    if (diag_line)
        graphics::lines(c(1,n),c(n,1)) # diagonal line

    # add subpop labels if present
    if (!is.null(labs))
        print_labels_multi(labs, labs_cex, labs_las, labs_line, labs_lwd, labs_sep, labs_even, labs_ticks, labs_text, labs_col)
    
    breaks # return breaks, since legends need it!
}

print_labels_multi <- function(labs, labs_cex, labs_las, labs_line, labs_lwd, labs_sep, labs_even, labs_ticks, labs_text, labs_col) {
    # normalize so we can loop over cases (assume arbitrary label levels)
    if (!is.matrix(labs))
        labs <- cbind(labs) # a col vector
    
    n <- ncol(labs)
    if (n > 1) {
        # expand scalar options as needed
        labs_cex <- rep_check(labs_cex, n)
        labs_las <- rep_check(labs_las, n)
        labs_line <- rep_check(labs_line, n)
        labs_lwd <- rep_check(labs_lwd, n)
        labs_col <- rep_check(labs_col, n)
        labs_sep <- rep_check(labs_sep, n)
        labs_even <- rep_check(labs_even, n)
        labs_ticks <- rep_check(labs_ticks, n)
        labs_text <- rep_check(labs_text, n)
    }
    # loop through now
    for (i in 1:n) {
        print_labels(
            labs[,i],
            cex = labs_cex[i],
            las = labs_las[i],
            line = labs_line[i],
            lwd = labs_lwd[i],
            sep = labs_sep[i],
            even = labs_even[i],
            ticks = labs_ticks[i],
            text = labs_text[i],
            col = labs_col[i]
        )
    }
}

rep_check <- function(vals, n) {
    # expand scalar options as needed
    if (length(vals) == 1) {
        vals <- rep.int(vals, n)
    } else if (length(vals) != n)
        stop('`', deparse(substitute(vals)), '` does not have the correct length (', length(vals) ,' != ', n, ')')
    vals # return repeated thing or original as needed (guaranteed to be length n)
}

rep_check_list <- function(vals, n) {
    # expand objects as needed
    if (class(vals) != 'list') {
        vals <- rep(list(vals), n) # this is the desired transformation
    } else if (length(vals) != n)
        stop('`', deparse(substitute(vals)), '` does not have the correct length (', length(vals) ,' != ', n, ')')
    vals # return repeated thing or original as needed (guaranteed to be length n)
}

# Plot color key for heatmap
#
# This plots an image with a single column, showing the relationship between colors and values in the heatmap.
# The image fills the entire panel; use \code{\link[graphics]{layout}} to make sure this panel is scaled correctly, usually so it is taller rather than wider.
#
# This function is provided for users that want greater flexibility in creating plot layouts.
# However, \code{\link{plot_popkin}} will be easier to use and should suffice in most cases, please consider using that before calling this function directly.
# 
# Note \code{\link{plot_popkin_single}} construct breaks that are symmetric about zero, which ensures that the middle color (white) corresponds to the zero kinship.
# In contrast, kinship_range need not be symmetric and it is preferably the true range of the data.
#
# @param breaks The vector of \eqn{n+1} values at which colors switch, as returned by \code{\link{plot_popkin_single}}
# @param label The name of the variable that the colors measure (i.e. "Kinship")
# @param col Color vector of length \eqn{n}.  Default colors are a progression from blue to white to red obtained from RColorBrewer.
# @param kinship_range Range of the color key, preferably the range of the data (default is infered from the breaks, but they need not agree)
# @param leg_n The desired number of ticks in the y-axis (input to \code{\link{pretty}}, see that for more details)
#
# @return Nothing
#
# @examples
# # Construct toy data
# X <- matrix(c(0,1,2,1,0,1,1,0,2), nrow=3, byrow=TRUE) # genotype matrix
# subpops <- c(1,1,2) # subpopulation assignments for individuals
# 
# # NOTE: for BED-formatted input, use BEDMatrix!
# # "file" is path to BED file (excluding .bed extension)
# # library(BEDMatrix)
# # X <- BEDMatrix(file) # load genotype matrix object
#
# # estimate the kinship matrix "kinship" from the genotypes "X"!
# kinship <- popkin(X, subpops) # calculate kinship from X and optional subpop labels
# kinship <- inbr_diag(kinship) # transform diagonal permanently
#
# # suppose kinship is the only kinship matrix to plot, then...
#
# # VERSION 1: direct plot through plot_popkin (easiest but restricts layout/margin choices)
# plot_popkin(kinship)
#
# # VERSION 2: detailed construction that calls plot_popkin_single and heatmap_legend directly, setting layout and other features manually
# # equals version 1 above but demonstrates the steps that need to be carried out
# range_real <- range(kinship) # get the real data range first
# # also construct symmetric version of this range, so colors work out better...
# max_sym <- max(abs(range_real))
# range_sym <- c(-max_sym, max_sym)
# # start layout
# layout(1:2, widths=c(1,0.1))
# # plot the heatmap in first panel
# breaks <- plot_popkin_single(kinship, kinship_range=range_sym) # pass symmetric range here
# # add color key to second panel
# heatmap_legend(breaks, kinship_range=range_real) # pass real range here
# # label y-axis on outer margin
# mtext('Individuals', side=2, outer=TRUE)
#
# @export
heatmap_legend <- function(breaks, kinship_range = NULL, label = 'Kinship', col = NULL, col_n = 100, leg_n = 5) {
    # creates a nicer heatmap legend, but it has to be a standalone image (in its own panel, preferably through layout so it's a skinny panel)
    # this function fills panel, so here we don't set margins/etc (it's best left to the end user)

    if (is.null(col))
        col <- plot_popkin_palette(n = col_n) # default coloring
    
    if (!is.null(kinship_range)) {
        # here's a case where we only want to plot things within an altered range than that of breaks/col
        # expected application is kinship_range is real data range, while breaks/col are rigged to be wider (particularly because we forced them to be symmetric, and have a white color at zero)
        # so we'll proceed by assuming kinship_range is contained in breaks, but we might want to remove breaks and colors with them

        # find the bins where our desired range falls
        is <- cut(kinship_range, breaks, labels = FALSE, include.lowest = TRUE)
        # toss breaks that we didn't use
        breaks <- breaks[is[1]:(is[2]+1)] # we need to go one over for top
        col <- col[is[1]:is[2]] # colors don't need one over
    }
    # old processing follows...
    nb <- length(breaks) # length of breaks
    
    # first plot sequence of colors
    # NOTE: breaks[2:length(breaks)] plots the sequence of colors because (from ?image):
    # > intervals are closed on the right and open on the left except for the lowest interval which is closed at both ends.
    # so plotting all top values of breaks works!
    graphics::image(
                  z = matrix(breaks[2:nb], nrow = 1),
                  col = col,
                  breaks = breaks,
                  xaxt = "n",
                  yaxt = "n",
                  useRaster = TRUE
              )

    # now we should add axis numbers/ticks
    # this is a real pain, because 0 and 1 are in the middle of the first and last boxes
    # old code (originally from heatmap.2) assumed 0 and 1 were at the end of the boxes, not their middles!
    # given this "delta"
    delta <- 1/(2*(nb-2))
    # we actually want ends to be at c(-delta, 1+delta)
    lv <- pretty(breaks, n = leg_n)
    xv <- scale_delta(as.numeric(lv), breaks[1], breaks[nb], delta)
    graphics::axis(4, at = xv, labels = lv)
    
    # lastly, add axis label
    graphics::mtext(side = 4, label, line = 2)
}

# spreads data to c(0-delta, 1+delta)
scale_delta <- function(x, x_min = min(x), x_max = max(x), delta) {
    # this is the range we want in the output
    yMin <- -delta
    yMax <- 1 + delta
    # this is the transformed data to return!
    # this holds: (x = x_min -> y = yMin) and (x = x_max -> y = yMax)
    yMin + (yMax - yMin) * (x - x_min) / (x_max - x_min)
}

plot_popkin_layout <- function(n, nr) {
    # figure out layout given a requested number of rows
    # step 1: dimensions
    nr <- as.integer(nr) # in case it's not an integer...
    if (nr < 1) nr <- 1 # just treat as 1 (no complaining)
    if (nr > n) nr <- n # if we asked for more rows than data, set to data (again no complaining)
    nc <- ceiling(n/nr) # this is the correct number of columns (might have blank cells)
    nr <- ceiling(n/nc) # reset backwards in case the nr provided was too large (this will reduce empty space)
    # step 2: fill layout matrix to just work
    layout <- c(1:n, rep.int(0, nr*nc - n)) # first fill in all nr*nc values, including zeroes as needed
    layout <- matrix(layout, nrow = nr, ncol = nc, byrow = TRUE) # turn into matrix
    layout <- cbind(layout, c(n+1, rep.int(0, nr-1))) # add final column for color key and nothing else
    # step 3: set up widths vector too
    widths <- c(rep.int(1, nc), 0.1) # last column for color key is 10% the width of the rest
    # make layout now!
    graphics::layout(layout, widths = widths) # note rows are all equal height
}

plot_popkin_palette <- function(n = 100) {
    #hmcolBR <- rev(RColorBrewer::brewer.pal(11,"RdBu"))
    hmcolR <- RColorBrewer::brewer.pal(9,"Reds")
    hmcolB <- rev(RColorBrewer::brewer.pal(9,"Blues"))
    # create a hybrid...
    hmcolBR <- c(hmcolB[1:8], hmcolR)
    # now interpolate (now always done!)
    # first step creates a new function that performs the mapping
    hmcolBRfun <- grDevices::colorRampPalette(hmcolBR)
    # second step actually creates new palette with n colors
    hmcolBR <- hmcolBRfun(n) # overwrite initial palette
    return(hmcolBR)
}

label_boundaries <- function(labs) {
    # finds change of labels between adjacent positions
    # construct two vectors, one with labels, the other with the boundary points
    n <- length(labs)
    b <- c(1) # first boundary point is always 1
    l <- c(labs[1]) # first label in long list is first label in summary list
    ll <- labs[1] # last label, to detect changes
    # construct the rest of the vectors
    for (i in 2:n) {
        li <- labs[i]
        # keep going unless current label doesn't match last label!
        if (li != ll) {
            l <- c(l, li) # add new label
            b <- c(b, i) # add boundary point
            ll <- li # update last label for next iteration
        }
    }
    # the last boundary is the last point incremented by one (cause that's how we encoded boundaries in the middle)
    b <- c(b, n+1)
    # return as list
    list(labels = l, boundaries = b)
}

line_to_user <- function(line, side) {
    # awesome function from:
    # http://stackoverflow.com/questions/30765866/get-margin-line-locations-in-log-space/30835971#30835971
    lh <- graphics::par('cin')[2] * graphics::par('cex') * graphics::par('lheight')
    x_off <- diff(graphics::grconvertX(c(0, lh), 'inches', 'npc'))
    y_off <- diff(graphics::grconvertY(c(0, lh), 'inches', 'npc'))
    switch(side,
           `1` = graphics::grconvertY(-line * y_off, 'npc', 'user'),
           `2` = graphics::grconvertX(-line * x_off, 'npc', 'user'),
           `3` = graphics::grconvertY(1 + line * y_off, 'npc', 'user'),
           `4` = graphics::grconvertX(1 + line * x_off, 'npc', 'user'),
           stop("Side must be 1, 2, 3, or 4", call. = FALSE))
}

panel_letter <- function(letter, cex = 1.5, line = 0.5, adj = 0) {
    # a wrapper around mtext with useful default values
    graphics::mtext(letter, cex = cex, line = line, adj = adj)
}

print_labels <- function(
                         labs,
                         x = NULL,
                         doMat = TRUE,
                         cex = 1,
                         las = 0,
                         lwd = 1,
                         sep = TRUE,
                         ticks = FALSE,
                         even = FALSE,
                         line = 0,
                         text = TRUE,
                         col = 'black',
                         side1 = 1,
                         side2 = 2
                         ) {
    labsObj <- label_boundaries(labs)
    # extract the data (smaller var names)
    l <- labsObj$labels
    b <- labsObj$boundaries
    n <- length(l) # number of labels (or number of boundaries minus one)
    m <- max(b) # number of individuals plus 1 (usually n+1, but here I messed up notation above, meh)
    # for non-barplots, the sensible default is to use the indexes as coordinates (or do they have to be normalized?)
    if (is.null(x)) {
        x <- 1:m
    } else {
        # a hack necessary for barplots at least
        # in this case x is missing the x[m] element, fill it in!
        x[m] <- x[m-1] + (x[2]-x[1]) # extend by usual gap (assuming it's regular)
    }
    x_min <- min(x)
    x_max <- max(x) # used mostly to reflect coordinates, equals at[n+1]

    # positions of irregular boundaries
    gapX <- (x[2]-x[1])/2 # shared by ticks and lines
    at <- x[b] - gapX # otherwise things are placed in the middle, shift by one half of a pixel

    # add ticks using axis()
    if (ticks) {
        graphics::axis(1, at = at, labels = FALSE, lwd = lwd)
        if (doMat)
            graphics::axis(2, at = x_max-at, labels = FALSE, lwd = lwd)
    }
    
    if (sep) {
        # draw black horizontal lines at every boundary (including ends, looks weird otherwise)
        # assume regular spacing
        graphics::abline(v = at, lwd = lwd, col = col)
        if (doMat)
            graphics::abline(h = x_max-at, lwd = lwd, col = col)
    }

    # label placement, ticky connector line calcs for "even" case
    if (even) {
        # hard case

        # place labels equally spaced on x's range...
        # positions (length of labels)
        gapY <- (x_max - x_min)/n # first divide range into even segments
        y <- x_min + gapY*((1:n) - 0.5) # sequence of middles: length n
        y2 <- x_min + gapY*(0:n) - gapX # sequence of boundaries: length n+1
        # NOTE: y[1] == x_min + gapY/2 # so it starts in middle of first bin
        # NOTE: y[n] == x_max - gapY/2 # so it ends in the middle of last bin
        # NOTE: y2[1] == x_min - gapX and y2[n+1] == x_max - gapX match boundaries with half pixel shift
        
        # connect label boundaries to irregular boundaries in plot
        ysLines <- line_to_user(c(0,line), 1) # shared by every label boundary on x-axis
        xsLines <- line_to_user(c(0,line), 2) # shared by every label boundary on y-axis
        for (i in 1:(n+1)) {
            # middle coordinate is a bit messy: first get ith and ith+1 boundaries, then turn them to coordinates, then average
            # don't forget that end is first element of next group, so we always need to reduce it by one
            xi <- at[i] # coincides with tick positions
            yi <- y2[i] # boundary of words
            graphics::lines(c(xi, yi), ysLines, lwd = lwd, xpd = NA) # draw line for x-axis
            if (doMat)
                graphics::lines(xsLines, x_max - c(xi, yi), lwd = lwd, xpd = NA)
        }
        
    } else {
        # put labels in middle of each range
        y <- (at[1:n] + at[2:(n+1)])/2 # sequence of label positions
    }
    
    if (text) {
        # place labels at "y" (set according to "even = TRUE" or otherwise)
        graphics::mtext(l, side = side1, at = y, cex = cex, las = las, line = line)
        # for matrices, do both ways!
        if (doMat)
            graphics::mtext(l, side = side2, at = x_max-y, cex = cex, las = las, line = line)
    }
}
