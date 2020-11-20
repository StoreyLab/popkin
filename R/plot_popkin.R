#' Visualize one or more kinship matrices
#'
#' This function plots one or more kinship matrices and a shared legend for the color key.
#' Many options allow for fine control of individual or subpopulation labeling.
#' This code assumes input matrices are symmetric.
#'
#' `plot_popkin` plots the input kinship matrices as-is.
#' For best results, a standard kinship matrix (such as the output of `\link{popkin}`) should have its diagonal rescaled to contain inbreeding coefficients (`\link{inbr_diag}` does this) before `plot_popkin` is used.
#'
#' This function permits the labeling of individuals (from row and column names when `names = TRUE`) and of subpopulations (passed through `labs`).
#' The difference is that the labels passed through `labs` are assumed to be shared by many individuals, and lines (or other optional visual aids) are added to demarcate these subgroups.
#'
#' @param kinship A numeric kinship matrix or a list of matrices.
#' Note `kinship` may contain `NULL` elements (makes blank plots with titles; good for placeholders or non-existent data)
#' @param titles Titles to add to each matrix panel (default is no titles)
#' @param col Colors for heatmap (default is a red-white-blue palette symmetric about zero constructed using RColorBrewer).
#' @param col_n The number of colors to use in the heatmap (applies if `col = NULL`).
#' @param mar Margins shared by all panels (if a vector) or for each panel (if a list of such vectors).
#' If the vector has length 1, `mar` corresponds to the shared lower and left margins, while the top and right margins are set to zero.
#' If this length is 2, `mar[1]` is the same as above, while `mar[2]` is the top margin.
#' If this length is 4, then `mar` is a fully-specified margin vector in the standard format `c(bottom, left, top, right)` that `\link[graphics]{par}('mar')` expects.
#' Vectors of invalid lengths produce a warning.
#' Note the padding `mar_pad` below is added to every margin if set.
#' If `NULL`, the original margin values are used without change, and are reset for every panel that has a `NULL` value.
#' The original margins are also reset after plotting is complete.
#' @param mar_pad Margin padding added to all panels (`mar` above and `leg_mar` below).
#' Default 0.2.
#' Must be a scalar or a vector of length 4 to match `\link[graphics]{par}('mar')`.
#' @param oma Outer margin vector.
#' If length 1, the value of `oma` is applied to the left outer margin only (so `ylab` below displays correctly) and zero outer margins elsewhere.
#' If length 4, all outer margins are expected in standard format `\link[graphics]{par}('mar')` expects (see `mar` above).
#' `mar_pad` above is never added to outer margins.
#' If `NULL`, no outer margins are set (previous settings are preserved).
#' Vectors of invalid lengths produce a warning.
#' @param diag_line If `TRUE` adds a line along the diagonal (default no line).
#' May also be a vector of logicals to set per panel (lengths must agree).
#' @param panel_letters Vector of strings for labeling panels (default A-Z).
#' No labels are added if `NULL`, or when there is only one panel except if its set to a single letter in that case (this behavior is useful if goal is to have multiple external panels but popkin only creates one of these panels).
#' @param panel_letters_cex Scaling factor of panel letters (default 1.5).
#' @param null_panel_data
#' If `FALSE` (default), panels with `NULL` kinship matrices must not have titles or other parameters set, and no panel letters are used in these cases.
#' If `TRUE`, panels with `NULL` kinship matrices must have titles and other parameters set.
#' In the latter case, these `NULL` panels also get panel letters.
#' The difference is important when checking that lengths of non-singleton parameters agree.
#' @param weights A vector with weights for every individual, or a list of such vectors if they vary per panel.
#' The width of every individual becomes proportional to their weight.
#' Individuals with zero or negative weights are omitted.
#' @param raster A logical equivalent to `useRaster` option in the `image` function used internally, or a vector of such logicals if the choice varies per panel.
#' If `weights` are non-`NULL` in a given panel, `raster = FALSE` is forced (this is necessary to plot images where columns and rows have variable width).
#' If `weights` are `NULL`, the default is `raster = TRUE`, but in this case the user may override (for example, so panels are visually coherent when some use weights while others do not, as there are small differences in rendering implementation for each value of `raster`).
#' Note that a multipanel figure with a list of `weights` sets `raster = FALSE` to all panels by default, even if the weights were only applied to a subset of panels.
#' @param sym If `FALSE` (default), plots non-symmetric (but square) kinship matrices without issues.  If `TRUE`, stops if any input kinship matrices are not symmetric.
#'
#' AXIS LABEL OPTIONS
#' 
#' @param ylab The y-axis label (default "Individuals").
#' If `length(ylab) == 1`, the label is placed in the outer margin (shared across panels);
#' otherwise `length(ylab)` must equal the number of panels and each label is placed in the inner margin of the respective panel.
#' @param ylab_adj The value of "adj" passed to `\link[graphics]{mtext}`.
#' If `length(ylab) == 1`, only the first value is used, otherwise `length(ylab_adj)` must equal the number of panels.
#' @param ylab_line The value of "line" passed to `\link[graphics]{mtext}`.
#' If `length(ylab) == 1`, only the first value is used, otherwise `length(ylab_line)` must equal the number of panels.
#' 
#' LAYOUT OPTIONS
#' 
#' @param layout_add If `TRUE` (default) then `\link[graphics]{layout}` is called internally with appropriate values for the required number of panels for each matrix, the desired number of rows (see `layout_rows` below) plus the color key legend.
#' The original layout is reset when plotting is complete and if `layout_add = TRUE`.
#' If a non-standard layout or additional panels (beyond those provided by `plot_popkin`) are desired, set to `FALSE` and call `\link[graphics]{layout}` yourself beforehand.
#' @param layout_rows Number of rows in layout, used only if `layout_add = TRUE`.
#'
#' LEGEND (COLOR KEY) OPTIONS
#' 
#' @param leg_per_panel If `TRUE`, every kinship matrix get its own legend/color key (best for matrices with very different scales).
#' If `FALSE` (default), a single legend/color key is shared by all kinship matrix panels.
#' @param leg_title The name of the variable that the heatmap colors measure (default "Kinship"), or a vector of such values if they vary per panel.
#' @param leg_cex Scaling factor for `leg_title` (default 1), or a vector of such values if they vary per panel.
#' @param leg_n The desired number of ticks in the legend y-axis (input to `\link{pretty}`, see that for more details), or a vector of such values if they vary per panel.
#' @param leg_width The width of the legend panel, relative to the width of the kinship panel.
#' This value is passed to `\link[graphics]{layout}` (ignored if `layout_add = FALSE`).
#' @param leg_mar Margin values for the legend panel only, or a list of such values if they vary per panel.
#' A length-4 vector (in `c( bottom, left, top, right )` format that `\link[graphics]{par}('mar')` expects) specifies the full margins, to which `mar_pad` is added.
#' Otherwise, the margins used in the last panel are preserved with the exception that the left margin is set to zero, and if `leg_mar` is length-1, it is used to specify the right margin (plus the value of `mar_pad`, see above).
#'
#' INDIVIDUAL LABEL OPTIONS
#' 
#' @param names If `TRUE`, the column and row names are plotted in the heatmap, or a vector of such values if they vary per panel.
#' @param names_cex Scaling factor for the column and row names, or a vector of such values if they vary per panel.
#' @param names_line Line where column and row names are placed, or a vector of such values if they vary per panel.
#' @param names_las Orientation of labels relative to axis.
#' Default (2) makes labels perpendicular to axis.
#'
#' SUBPOPULATION LABEL OPTIONS
#' 
#' @param labs Subpopulation labels for individuals.
#' Use a matrix of labels to show groupings at more than one level (for a hierarchy or otherwise).
#' If input is a vector or a matrix, the same subpopulation labels are shown for every heatmap panel; the input must be a list of such vectors or matrices if the labels vary per panel.
#' @param labs_cex A vector of label scaling factors for each level of labs, or a list of such vectors if labels vary per panel.
#' @param labs_las A vector of label orientations (in format that `\link[graphics]{mtext}` expects) for each level of labs, or a list of such vectors if labels vary per panel.
#' @param labs_line A vector of lines where labels are placed (in format that `\link[graphics]{mtext}` expects) for each level of labs, or a list of such vectors if labels vary per panel.
#' @param labs_sep A vector of logicals that specify whether lines separating the subpopulations are drawn for each level of labs, or a list of such vectors if labels vary per panel.
#' @param labs_lwd A vector of line widths for the lines that divide subpopulations (if `labs_sep = TRUE`) for each level of labs, or a list of such vectors if labels vary per panel.
#' @param labs_col A vector of colors for the lines that divide subpopulations (if `labs_sep = TRUE`) for each level of labs, or a list of such vectors if labels vary per panel.
#' @param labs_ticks A vector of logicals that specify whether ticks separating the subpopulations are drawn for each level of labs, or a list of such vectors if labels vary per panel.
#' @param labs_text A vector of logicals that specify whether the subpopulation labels are shown for each level of labs, or a list of such vectors if labels vary per panel.
#' Useful for including separating lines or ticks without text.
#' @param labs_even A vector of logicals that specify whether the subpopulations labels are drawn with equal spacing for each level of labs, or a list of such vectors if labels vary per panel.
#' When `TRUE`, lines mapping the equally-spaced labels to the unequally-spaced subsections of the heatmap are also drawn.
#' 
#' @param ... Additional options passed to `\link[graphics]{image}`.
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
                        oma = 1.5,
                        diag_line = FALSE,
                        panel_letters = toupper(letters),
                        panel_letters_cex = 1.5,
                        ylab = 'Individuals',
                        ylab_adj = NA,
                        ylab_line = 0,
                        layout_add = TRUE,
                        layout_rows = 1, 
                        leg_per_panel = FALSE,
                        leg_title = 'Kinship',
                        leg_cex = 1,
                        leg_n = 5,
                        leg_mar = 3,
                        leg_width = 0.3,
                        names = FALSE,
                        names_cex = 1,
                        names_line = NA,
                        names_las = 2,
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
                        null_panel_data = FALSE,
                        weights = NULL,
                        raster = is.null(weights),
                        sym = FALSE, # defaults to plotting non-symmetric matrices too
                        ...
                        ) {
    # wrapper around individual panels and color key
    # does not set PDF output, margins, layout, etc
    # assumes kinship is a list of matrices, if not it is internally turned into one

    # NOTE: we'll let plot_popkin_single validate kinship matrices (no validation here)
    if (!is.list(kinship))
        kinship <- list(kinship) # turn into list so rest works
    
    # number of heatmap panels, including NULL cases
    # this should only be used for layout
    n_all <- length(kinship)
    
    # figure out which data are non-NULL
    indexes_not_null <- !sapply(kinship, is.null)
    if ( !any( indexes_not_null ) )
        stop('Every element of list "kinship" is NULL!')

    # decide which length to use for validations
    if (null_panel_data) {
        # NULL kinship panels require titles and other parameters set
        # (this is original behavior)
        n <- n_all
    } else {
        # NULL kinship panels must NOT have titles and other parameters set
        # most of the parameters have to agree with this value (n rather than n_all)
        n <- sum(indexes_not_null)
    }
    
    # start checking lengths, etc
    if (!is.null(titles)) {
        if (length(titles) != n)
            stop('`titles` provided are not the same length as data! Data: ', n, ', titles: ', length(titles))
    } else {
        titles <- rep.int('', n) # make blank titles of the same length as data
    }
    # check label lengths
    if (!is.null(labs) && is.list(labs) && length(labs) != n)
        stop('there are ', n, ' non-NULL panels but ', length(labs), ' label sets!')
    # expand other things that may vary per panel
    names <- rep_check(names, n)
    names_cex <- rep_check(names_cex, n)
    names_line <- rep_check(names_line, n)
    names_las <- rep_check(names_las, n)
    diag_line <- rep_check(diag_line, n)
    raster <- rep_check(raster, n)
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
    weights <- rep_check_list(weights, n)
    # ylab behavior is more dynamic!
    if (length(ylab) > 1) {
        ylab <- rep_check(ylab, n) # just makes sure length is n
        ylab_adj <- rep_check(ylab_adj, n)
        ylab_line <- rep_check(ylab_line, n)
    }
    # all of these can vary when legend is panel-specific!
    if ( leg_per_panel ) {
        leg_title <- rep_check( leg_title, n )
        leg_cex <- rep_check( leg_cex, n )
        leg_n <- rep_check( leg_n, n )
        leg_mar <- rep_check_list( leg_mar, n )
    }
    
    if ( !leg_per_panel ) {
        # this shared range is for a shared scale
        # if leg_per_panel is TRUE, leave these values NULL so we can catch bugs if they are used improperly 
        
        # code needs two versions of the range
        # - range_real is the real range, used in the end so the color key doesn't show values that weren't actually used
        # - range_sym is a symmetric range, used internally to ensure zero is in the exact middle (set to white in the default)
        # get range and construct symmetric range that helps with plotting nice colors with white at zero
        # range of all non-NULL data plotted
        range_real <- range( unlist( lapply( kinship[ indexes_not_null ], range, na.rm = TRUE ) ) )
        # these next few lines force symmetry for colors (looks better)
        max_sym <- max( abs( range_real ) )
        range_sym <- c( -max_sym, max_sym )
    }

    # save entire original setup, to reset in the end
    # no.readonly is since some parameters cannot be changed (trying to set them results in ugly warnings)
    par_orig <- graphics::par( no.readonly = TRUE )
    # save original margins, which may get reset per panel (separately of final reset)
    mar_orig <- graphics::par('mar')

    # set outer margins as desired
    if ( !is.null( oma ) ) {
        if ( length( oma ) == 1 ) {
            # when a single value, it is the left outer margin (all other values are zero)
            graphics::par( oma = c(0, oma, 0, 0) )
        } else if ( length( oma ) == 4 ) {
            # here the whole vector was specified, set that
            graphics::par( oma = oma )
        } else {
            warning("`oma` has an invalid length (was not 1 or 4): ", length( oma ), "\nOuter margins were unchanged!")
        }
    }
    
    if (layout_add) {
        # ...
        
        # figure out layout given a requested number of rows
        plot_popkin_layout(
            n = n_all,
            nr = layout_rows,
            leg_per_panel = leg_per_panel,
            leg_width = leg_width
        )
    }

    # breaks of all following plots should match!
    breaks <- NULL
    # index of non-null data
    i <- 0 # start here to update in the beginning
    # navigate *_all version (includes NULL cases)
    for (i_all in 1:n_all) {
        if (null_panel_data) {
            # these two are the same in this original version
            i <- i_all
            # panel will be created as it used to, even for NULL cases...
        } else {
            if ( !indexes_not_null[ i_all ] ) {
                # This is a null case
                # Let's just make the blank panel and move on
                # Important: do not add panel letters or titles!
                graphics::plot.new()
                # if every panel has a legend, then we need to skip one more space
                if ( leg_per_panel )
                    graphics::plot.new()
                # force advance to next i_all, not executing any of the rest
                next
            } else {
                # now we can increment this one, for all non-kinship data
                i <- i + 1
            }
        }

        if ( leg_per_panel ) {
            # same issues as before, but applied to the single given matrix
            range_real <- range( kinship[[ i_all ]], na.rm = TRUE )
            max_sym <- max( abs( range_real ) )
            range_sym <- c( -max_sym, max_sym )
        }
        
        mar_i <- mar[[ i ]]
        if ( !is.null( mar_i ) ) {
            # change margins if necessary!
            if ( length( mar_i ) == 1 ) {
                # this value is bottom and left margins (others are zero)
                graphics::par( mar = c(mar_i, mar_i, 0, 0) + mar_pad )
            } else if ( length( mar_i ) == 2 ) {
                # first value is same as above
                # second value is for top margin (to make room for optional titles)
                # right margin remains zero
                graphics::par( mar = c(mar_i[1], mar_i[1], mar_i[2], 0) + mar_pad )
            } else if ( length( mar_i ) == 4 ) {
                # full margin vector has been specified
                graphics::par( mar = mar_i + mar_pad )
            } else {
                warning("`mar` at panel ", i, " has an invalid length (was not 1, 2 or 4): ", length( mar_i ), "\nPrevious margins were reset!")
                # restore original margins otherwise!
                graphics::par( mar = mar_orig )
            }
        } else {
            # restore original margins otherwise!
            graphics::par( mar = mar_orig )
        }

        breaks_i <- plot_popkin_single(
            kinship[[ i_all ]],
            kinship_range = range_sym,
            col = col,
            col_n = col_n,
            names = names[i],
            names_cex = names_cex[i],
            names_line = names_line[i],
            names_las = names_las[i],
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
            weights = weights[[i]],
            raster = raster[i],
            sym = sym,
            ...
        )
        # don't overwrite for non-data kinship[[i]] cases
        if (!is.null(breaks_i))
            breaks <- breaks_i
        
        # add ylab for every panel when there is more than one choice, and provided it was non-NA
        # uses inner rather than outer margin (only choice that makes sense)
        if ( length(ylab) > 1 && !is.na( ylab[i] ) ) 
            graphics::mtext( ylab[i], side = 2, adj = ylab_adj[i], line = ylab_line[i] )
        
        # add letters only when ...
        if (
            # panel_letters is non-NULL
            !is.null(panel_letters) &&
            # and 
            (
                # there is more than one panel
                1 < n ||
                # or if there is one panel but we also set one letter only
                (n == 1 && length(panel_letters) == 1)
            ) &&
            # and if we haven't exceeded the available letters (to prevent stupid errors)
            i <= length(panel_letters)
        )
            panel_letter(panel_letters[i], cex = panel_letters_cex)

        if ( leg_per_panel ) {
            # add panel-specific legend/color key
            heatmap_legend(
                breaks,
                kinship_range = range_real,
                label = leg_title[ i ],
                col = col,
                col_n = col_n,
                leg_n = leg_n[ i ],
                cex = leg_cex[ i ],
                leg_mar = leg_mar[[ i ]],
                mar_pad = mar_pad
            )
        }
    }

    if (! leg_per_panel ) {
        # add shared legend/color key
        heatmap_legend(
            breaks,
            kinship_range = range_real,
            label = leg_title,
            col = col,
            col_n = col_n,
            leg_n = leg_n,
            cex = leg_cex,
            leg_mar = leg_mar,
            mar_pad = mar_pad
        )
    }

    # add margin only once if there was only one, place in outer margin (only choice that makes sense)
    if ( length(ylab) == 1 )
        graphics::mtext( ylab, side = 2, adj = ylab_adj, outer = TRUE, line = ylab_line )

    # restore original setup when done, but only if we created the default layout
    # otherwise the external layout gets reset, which is bad if we were not done adding panels
    if ( layout_add ) {
        graphics::par( par_orig )
    } else {
        # restore original margins only!
        graphics::par( mar = mar_orig )
    }
}

plot_popkin_single <- function (
                                kinship = NULL,
                                kinship_range = NULL,
                                col = NULL,
                                col_n = 100,
                                names = FALSE,
                                names_cex = 1,
                                names_line = NA,
                                names_las = 2,
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
                                weights = NULL,
                                raster = is.null(weights),
                                sym = FALSE, # defaults to plotting non-symmetric matrices too
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
        validate_kinship( kinship, sym = sym )
    }
    
    # further data validation
    if (is.null(col))
        col <- plot_popkin_palette(n = col_n) # default coloring
    if (length(col) == 1 && is.character(col)) 
        col <- get(col, mode = "function")

    # weight validation and processing
    n <- ncol(kinship)
    if (!is.null(weights)) {
        raster <- FALSE # force overwrite, don't let the user change this (it won't work so meh).
        if (n != length(weights))
            stop('length of `weights` (', length(weights), ') and `kinship` dimension (', n, ') must match!')
        # check for zero or negative weights
        indexes_ind_keep <- weights > 0
        if ( !all( indexes_ind_keep ) ) {
            # filter the data!
            # individuals with zero or negative weights will be excluded
            weights <- weights[ indexes_ind_keep ]
            kinship <- kinship[ indexes_ind_keep, indexes_ind_keep ]
            # renormalize weights (not strictly necessary, but nice to have a range of (0, 1)
            weights <- weights / sum( weights )
            # recompute number of individuals!
            n <- ncol(kinship)
        } else {
            # set to null so downstream code knows not to do anything
            indexes_ind_keep <- NULL
        }
        # this sets the boundaries of the data
        # start at zero
        xb <- cumsum( c(0, weights) )
        # y axis has to be reversed in order and in value, this does both
        yb <- 1 - rev(xb)
    } else {
        # NOTE: in this configuration these are bin centers, not boundaries
        ## x <- 1:n
        ## y <- x
        # version for uniform boundaries, c(0, 1) range!
        xb <- ( 0 : n ) / n
        yb <- xb
        # set to null so downstream code knows not to do anything
        indexes_ind_keep <- NULL
    }
    
    # get default range if needed
    if ( is.null( kinship_range ) ) {
        # the version we need is symmetric about zero
        max_sym <- max( abs( kinship ), na.rm = TRUE )
        kinship_range <- c( - max_sym, max_sym )
    }
    # figure out breaks for colors
    breaks <- seq(kinship_range[1], kinship_range[2], length = length(col) + 1)
    numcols <- length(breaks) - 1
    if ( is.function(col) ) 
        col <- col(numcols)
    # here replace potentially forced kinship_range with the true range of the data...
    # but only if data is more extended, don't change anything otherwise!
    rangeRaw <- range(kinship, na.rm = TRUE)
    if ( breaks[1] > rangeRaw[1] )
        breaks[1] <- rangeRaw[1]
    if ( breaks[numcols+1] < rangeRaw[2] )
        breaks[numcols+1] <- rangeRaw[2]
    
    # main plot!
    # need to reverse the rows of the kinship for plotting!
    # (keep it this way for names below)
    kinship <- kinship[ n:1, ]
    graphics::image(
                  xb,
                  yb,
                  t(kinship),
                  xlim = c(0, 1), # 0.5 + c(0, n),
                  ylim = c(0, 1), # 0.5 + c(0, n),
                  axes = FALSE,
                  xlab = xlab,
                  ylab = ylab,
                  col = col,
                  breaks = breaks,
                  useRaster = raster,
                  ...
              )
    if (names) {
        # if we want to show labels, use the row/col names to add to picture
        # compute bin centers from boundaries (whether they were weighted or not)
        xc <- centers_from_boundaries(xb)
        yc <- centers_from_boundaries(yb)
        graphics::axis(1, xc, colnames(kinship), las = names_las, cex.axis = names_cex, tick = FALSE, line = names_line)
        graphics::axis(2, yc, rownames(kinship), las = names_las, cex.axis = names_cex, tick = FALSE, line = names_line)
    }
    if (diag_line)
        # diagonal line, version for c(1, n) range
        # (was wrong as (1,n) were centerpoints, not boundaries)
        ## graphics::lines(c(1,n),c(n,1))
        # diagonal line, version for c(0, 1) range
        graphics::lines( c(0, 1), c(1, 0) )

    # add subpop labels if present
    if (!is.null(labs))
        print_labels_multi(
            labs = labs,
            labs_cex = labs_cex,
            labs_las = labs_las,
            labs_line = labs_line,
            labs_lwd = labs_lwd,
            labs_sep = labs_sep,
            labs_even = labs_even,
            labs_ticks = labs_ticks,
            labs_text = labs_text,
            labs_col = labs_col,
            xb_ind = xb,
            indexes_ind_keep = indexes_ind_keep
        )

    # return breaks, since legends need it!
    # return invisibly so we can run this one standalone and it doesn't print `breaks` awkwardly
    invisible( breaks )
}

centers_from_boundaries <- function(xb) {
    # number of centers equals number of boundaries minus one
    n <- length(xb) - 1
    # indexes facilitate this
    indexes <- 1 : n
    # this is the value we want
    xc <- ( xb[ indexes ] + xb[ indexes + 1 ] ) / 2
    return( xc )
}

print_labels_multi <- function(labs, labs_cex, labs_las, labs_line, labs_lwd, labs_sep, labs_even, labs_ticks, labs_text, labs_col, xb_ind, indexes_ind_keep) {
    # normalize so we can loop over cases (assume arbitrary label levels)
    if (!is.matrix(labs))
        labs <- cbind(labs) # a col vector

    # now filter labs by individuals, as needed
    if ( !is.null( indexes_ind_keep ) )
        labs <- labs[ indexes_ind_keep, , drop = FALSE ]

    # it's good to check numbers of individuals now
    if (nrow(labs) != length(xb_ind) - 1)
        stop('Number of individuals disagrees between `labs` (', nrow(labs), ') and `xb_ind` (', length(xb_ind) - 1, ')!')

    # NOTE: here n is the number of label levels (not number of individuals, as usual)
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
            labs[, i],
            cex = labs_cex[i],
            las = labs_las[i],
            line = labs_line[i],
            lwd = labs_lwd[i],
            sep = labs_sep[i],
            even = labs_even[i],
            ticks = labs_ticks[i],
            text = labs_text[i],
            col = labs_col[i],
            # same for all levels (no [i])
            xb_ind = xb_ind
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
    if ( !is.list(vals) ) {
        vals <- rep(list(vals), n) # this is the desired transformation
    } else if (length(vals) != n)
        stop('`', deparse(substitute(vals)), '` does not have the correct length (', length(vals) ,' != ', n, ')')
    vals # return repeated thing or original as needed (guaranteed to be length n)
}

# Plot color key for heatmap
#
# This plots an image with a single column, showing the relationship between colors and values in the heatmap.
# The image fills the entire panel; use `\link[graphics]{layout}` to make sure this panel is scaled correctly, usually so it is taller rather than wider.
#
# This function is provided for users that want greater flexibility in creating plot layouts.
# However, `\link{plot_popkin}` will be easier to use and should suffice in most cases, please consider using that before calling this function directly.
# 
# Note `\link{plot_popkin_single}` construct breaks that are symmetric about zero, which ensures that the middle color (white) corresponds to the zero kinship.
# In contrast, kinship_range need not be symmetric and it is preferably the true range of the data.
#
# @param breaks The vector of `n+1` values at which colors switch, as returned by `\link{plot_popkin_single}`
# @param label The name of the variable that the colors measure (i.e. "Kinship")
# @param col Color vector of length `n`.  Default colors are a progression from blue to white to red obtained from RColorBrewer.
# @param kinship_range Range of the color key, preferably the range of the data (default is infered from the breaks, but they need not agree)
# @param leg_n The desired number of ticks in the y-axis (input to `\link{pretty}`, see that for more details)
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
# layout(1:2, widths=c(1, 0.3))
# # plot the heatmap in first panel
# breaks <- plot_popkin_single(kinship, kinship_range=range_sym) # pass symmetric range here
# # add color key to second panel
# heatmap_legend(breaks, kinship_range=range_real) # pass real range here
# # label y-axis on outer margin
# mtext('Individuals', side=2, outer=TRUE)
#
# @export
heatmap_legend <- function(
                           breaks,
                           kinship_range = NULL,
                           label = 'Kinship',
                           col = NULL,
                           col_n = 100,
                           leg_n = 5,
                           cex = NA,
                           leg_mar = 3,
                           mar_pad = 0.2
                           ) {
    # creates a nicer heatmap legend, but it has to be a standalone image (in its own panel, preferably through layout so it's a skinny panel)
    # this function fills panel, so here we don't set margins/etc (it's best left to the end user)

    if ( is.null( col ) )
        col <- plot_popkin_palette(n = col_n) # default coloring
    
    if ( !is.null( kinship_range ) ) {
        # here's a case where we only want to plot things within an altered range than that of breaks/col
        # expected application is kinship_range is real data range, while breaks/col are rigged to be wider (particularly because we forced them to be symmetric, and have a white color at zero)
        # so we'll proceed by assuming kinship_range is contained in breaks, but we might want to remove breaks and colors with them

        # find the bins where our desired range falls
        indexes <- cut( kinship_range, breaks, labels = FALSE, include.lowest = TRUE )
        # toss breaks that we didn't use
        breaks <- breaks[ indexes[1] : ( indexes[2] + 1 ) ] # we need to go one over for top
        col <- col[ indexes[1] : indexes[2] ] # colors don't need one over
    }
    # old processing follows...
    nb <- length(breaks) # length of breaks

    # change margins as necessary!
    if ( is.null(leg_mar) || length( leg_mar ) == 1 ) {
        # change the current left margin to the padding value
        # last margins
        mar_tmp <- graphics::par('mar')
        # replace left margin with zero plus pad
        mar_tmp[2] <- mar_pad
        # this value sets the right margin
        if ( length( leg_mar ) == 1 )
            mar_tmp[4] <- leg_mar + mar_pad
        # otherwise the right margin is preserved
        
        # update margins
        graphics::par( mar = mar_tmp )
        
    } else if ( length( leg_mar ) == 4 ) {
        # full specification
        graphics::par( mar = leg_mar + mar_pad )
    } else 
        stop('Invalid length for `leg_mar` (expected NULL, 1, or 4): ', length( leg_mar ) )
    
    # first plot sequence of colors
    # NOTE: breaks[2:length(breaks)] plots the sequence of colors because (from ?image):
    # > intervals are closed on the right and open on the left except for the lowest interval which is closed at both ends.
    # so plotting all top values of breaks works!
    graphics::image(
                  z = matrix( breaks[ 2 : nb ], nrow = 1 ),
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
    delta <- 1 / ( 2 * ( nb - 2 ) )
    # we actually want ends to be at c(-delta, 1+delta)
    lv <- pretty( breaks, n = leg_n )
    xv <- scale_delta( as.numeric(lv), breaks[1], breaks[nb], delta )
    graphics::axis( 4, at = xv, labels = lv )
    
    # lastly, add axis label
    graphics::mtext( side = 4, label, line = 2, cex = cex )
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

calc_leg_width_min <- function
(
    leg_mar = 3,
    mar_pad = 0.2,
    n_col = 1,
    n_leg = 1,
    gap_min_lines = 1
) {
    # to avoid errors and user frustration, let's calculate the minimum leg_width that will result in a nice-looking plot.
    # in particular, we'll assume that the plot area for the legend should be at minimum one "line"
    
    # inner width, equals total width minus outer margins, in inches
    width_inner <- graphics::par('pin')[1]
    # number of lines that each inch has (default 5, is supposed to be read only but meh)
    lines_per_inch <- 1 / graphics::par('csi')
    # only part missing is default margins for legend

    # when margins vary per panel, this won't work... maybe we should disallow this earlier (TODO DOCS)
    if ( is.list( leg_mar ) )
        return( NA )

    # also the default of "preserving the last margin" is crap here
    # in that case the margin could vary per panel too
    if ( is.null(leg_mar) )
        return( NA )
    
    # following code assumes leg_mar is the same across panels
    if ( length( leg_mar ) == 1 ) {
        # value of right margin (left is zero)
        leg_mar_tot <- leg_mar
    } else if ( length( leg_mar ) == 4 ) {
        # full specification
        leg_mar_tot <- leg_mar[2] + leg_mar[4]
    } else 
        stop('Invalid length for `leg_mar` (expected 1 or 4): ', length( leg_mar ) )
    # add padding
    leg_mar_tot <- leg_mar_tot + 2 * mar_pad

    # final calculation
    leg_width_min <- 
        ( leg_mar_tot + gap_min_lines ) * n_col /
        ( width_inner * lines_per_inch - ( leg_mar_tot + gap_min_lines ) * n_leg )

    # this could be negative in an extreme case
    # let's stop with an informative message
    if ( leg_width_min < 0 )
        stop('Figure width is too small to acomodate all panels.  Please increase this width!')
    
    # return that value
    return( leg_width_min )
}

plot_popkin_layout <- function(n, nr = 1, leg_per_panel = FALSE, leg_width = 0.3) {
    # figure out layout given a requested number of rows
    
    # step 1: dimensions
    # in case it's not an integer...
    nr <- as.integer(nr)
    # just treat as 1 (no complaining)
    if (nr < 1)
        nr <- 1
    # if we asked for more rows than data, set to data (again no complaining)
    if (nr > n)
        nr <- n
    # this is the correct number of columns (might have blank cells)
    nc <- ceiling( n / nr )
    # reset backwards in case the nr provided was too large (this will reduce empty space)
    nr <- ceiling( n / nc )
    # number of blank panels
    nb <- nr * nc - n

    # using the given information, get the minimum legend width!
    # get more values from main function
    # TODO!!!  (currently unused)
    leg_width_min <- calc_leg_width_min(
        leg_mar = 3,
        mar_pad = 0.2,
        n_col = nc,
        n_leg = if (leg_per_panel) nc else 1,
        gap_min_lines = 1
    )
    # what to do when the above fails for whatever reason?
    # we only needed with default layout, so is there a scenario where a custom layout allows this all to go away?
    # (i.e. the restrictions on leg_mar)
    
    if ( !leg_per_panel ) {
        # this is the default version
        
        # step 2: fill layout matrix to just work
        # first fill in all nr*nc values, including zeroes as needed
        layout <- c(
            1 : n,
            rep.int( 0, nb )
        )
        # turn into matrix
        layout <- matrix(
            layout,
            nrow = nr,
            ncol = nc,
            byrow = TRUE
        )
        # add final column for color key and nothing else
        layout <- cbind(
            layout,
            c( n + 1 , rep.int( 0, nr - 1 ) )
        )
        
        # step 3: set up widths vector too
        # last column for color key is `leg_width` fraction of the width of the rest
        widths <- c( rep.int(1, nc), leg_width )
    } else {
        # every panel has a legend, so there's 2*n panels actually
        # still need zeroes (twice as many as default too)
        layout <- c(
            1 : ( 2 * n ),
            rep.int( 0, 2 * nb )
        )
        # turn into matrix
        # again note number of columns is double the default
        layout <- matrix(
            layout,
            nrow = nr,
            ncol = 2 * nc,
            byrow = TRUE
        )
        # but we don't add the final column for the color key (all those columns have been added)
        
        # step 3: set up widths vector too
        # every kinship matrix has a full panel followed by leg_width legend/color key panel
        widths <- rep.int( c(1, leg_width), nc)
    }
    
    # make layout now!
    # note rows are all equal height
    graphics::layout( layout, widths = widths )
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
                         x = NULL, # midpoints of individuals (in practice, given for barplots only)
                         xb_ind = NULL, # boundaries between individuals (NEW)
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
    # check only required input
    if (missing(labs))
        stop('`labs` is required!')

    # numbers of individuals, official version is from `labs`
    n_ind <- length(labs)
    
    labsObj <- label_boundaries(labs)
    # extract the data (smaller var names)
    # n_grp labels
    l <- labsObj$labels
    # n_grp + 1 boundary indexes, marking the first index of each group
    indexes_boundaries <- labsObj$boundaries
    # number of groups / labels
    n_grp <- length(l)

    if ( !is.null(x) ) {
        if ( length(x) != n_ind )
            stop('Length of individual centers `x` (', length(x), ') disagrees with `labs` length (', n_ind, ')!')
        # only old barplots trigger this
        # what is passed is centers only, regularly spaced
        xc_ind <- x # copy this
        # get regular spacing interval
        x_delta <- x[ 2L ] - x[ 1L ]
        # infer boundaries...
        # start by copying this again
        xb_ind <- xc_ind
        # then add the middle of the next entry after the current last one
        xb_ind[ n_ind + 1L ] <- xb_ind[ n_ind ] + x_delta
        # lastly, shift everything so that they're boundaries, not centers
        xb_ind <- xb_ind - x_delta / 2
    } else {
        # passing boundary coordinates is most ideal, but if missing we'll have to choose uniform boundaries
        if ( is.null( xb_ind ) ) {
            # (n_ind + 1) boundaries in c(0, 1)
            xb_ind <- ( 0L : n_ind ) / n_ind
        } else {
            # check xb_ind if it was passed
            if ( length(xb_ind) - 1L != n_ind )
                stop('Length of individual boundaries `xb_ind` (', length(xb_ind), ') disagrees with `labs` length (', n_ind, ')!  Expected `n+1` boundaries.')
        }
        # always infer centers for individuals from this
        xc_ind <- centers_from_boundaries(xb_ind)
    }
    
    # group boundary positions
    xb_grp <- xb_ind[ indexes_boundaries ]
    # get extrema, to space things out in `even` case
    xb_min <- xb_grp[ 1L ]
    # ... and this one also used to reflect boundaries
    xb_max <- xb_grp[ n_grp + 1 ]
    if (doMat) {
        # an assumption of our "boundary reflection" is that the first boundary is at zero, let's be sure of that
        if ( xb_min != 0 )
            stop('`xb_min == 0` is required for `doMat` trick to work.  Got xb_min: ', xb_min)
    }
    
    # add ticks to boundaries between groups using axis()
    if (ticks) {
        graphics::axis(1, at = xb_grp, labels = FALSE, lwd = lwd)
        if (doMat)
            graphics::axis(2, at = xb_max - xb_grp, labels = FALSE, lwd = lwd)
    }
    
    if (sep) {
        # draw black horizontal lines at every boundary between groups (including ends, looks weird otherwise)
        # assume regular spacing
        graphics::abline(v = xb_grp, lwd = lwd, col = col)
        if (doMat)
            graphics::abline(h = xb_max - xb_grp, lwd = lwd, col = col)
    }

    # label placement, ticky connector line calcs for "even" case
    if (even) {
        # hard case

        # even sequence of boundaries: length n+1
        xb_grp2 <- xb_min + ( 0 : n_grp ) / n_grp * xb_max
        # even centers calculated from the even boundaries
        # this is where the labels will be placed at the end
        xc_grp <- centers_from_boundaries( xb_grp2 )
            
        # connect label boundaries to irregular boundaries in plot
        ysLines <- line_to_user( c(0, line), 1 ) # shared by every label boundary on x-axis
        xsLines <- line_to_user( c(0, line), 2 ) # shared by every label boundary on y-axis
        for (i in 1 : ( n_grp + 1 ) ) {
            # middle coordinate is a bit messy: first get ith and ith+1 boundaries, then turn them to coordinates, then average
            # don't forget that end is first element of next group, so we always need to reduce it by one
            xi <- xb_grp[i] # coincides with tick positions
            yi <- xb_grp2[i] # boundary of words
            graphics::lines(c(xi, yi), ysLines, lwd = lwd, xpd = NA) # draw line for x-axis
            if (doMat)
                graphics::lines(xsLines, xb_max - c(xi, yi), lwd = lwd, xpd = NA)
        }
        
    } else {
        # put labels in middle of each group
        # calculate group center positions from their boundaries
        xc_grp <- centers_from_boundaries( xb_grp )
    }
    
    if (text) {
        # place labels at "xc_grp" (group centers, or set according to "even = TRUE")
        graphics::mtext(l, side = side1, at = xc_grp, cex = cex, las = las, line = line)
        # for matrices, do both ways!
        if (doMat)
            graphics::mtext(l, side = side2, at = xb_max - xc_grp, cex = cex, las = las, line = line)
    }
}
