#' Make a structure/admixture plot
#'
#' This function facilitates structure plots with options that resemble those of [plot_popkin()] in name and results.
#' The biggest difference is this function plots single panels (technically 2 panels including the legend, unless it is omitted), whereas [plot_popkin()] can plot multiple kinship matrices with a shared legend.
#'
#' @param Q The admixture proportions matrix, with `n` individuals along rows and `K` ancestries along columns.
#' Rows should sum to 1, but this is not enforced.
#' There must be at least 2 ancestries.
#' The ancestry labels used by the legend must be the column names, which are unlabeled if the column names are missing.
#' @param col A vector of at least `K` colors for the ancestries (extra colors are ignored).
#' By default uses the "Paired" palette of `RColorBrewer`, which has at most 12 colors, so please provide colors if `K > 12`.
#' Since the minimum number of colors for "Paired" is 3, when `K = 2` we ask for 3 colors, then remove the middle color internally.
#' @param mar_pad Margin padding used for legend panel only (margins for first/main panel are not altered by this function).
#' @param panel_letters Panel letter to include in first/main panel (default `NA` is no letter).
#' Despite name (matches [plot_popkin()]), must be scalar.
#' @param panel_letters_cex Scaling factor of panel letter (default 1.5).
#' @param panel_letters_adj X-axis adjustment for panel letter (default -0.1).
#' Negative values place the letter into the left margin area.
#' Might need adjustment depending on the size of the left margin.
#'
#' AXIS LABEL OPTIONS
#' 
#' @param axis_lab_cex Scaling factor for x-axis, y-axis, and legend title labels (which can also be set individually, see below).
#' @param xlab X-axis label (default "Individuals").
#' Set to `NA` to omit.
#' @param xlab_line The value of `line` for `xlab` passed to [graphics::mtext()].
#' @param xlab_cex Scaling factor for x-axis label.
#' @param ylab Y-axis label (default "Ancestry").
#' Set to `NA` to omit.
#' @param ylab_line The value of `line` for `ylab` passed to [graphics::mtext()].
#' @param ylab_side The value of `side` for `ylab` passed to [graphics::mtext()] (2 is y-axis, 1 is x-axis, can also place on top (3) or right (4)).
#' @param ylab_cex Scaling factor for y-axis label.
#'
#' LEGEND (COLOR KEY) OPTIONS
#' 
#' @param leg_title The name of the categorical ancestry variable (default "Ancestries").
#' @param leg_title_cex Scaling factor for legend title label.
#' @param leg_title_line The value of `line` for `leg_title` passed to [graphics::mtext()].
#' @param leg_cex Scaling factor for ancestry labels.
#' @param leg_width The width of the legend panel, relative to the width of the main panel.
#' This value is passed to [graphics::layout()] (ignored if `layout_add = FALSE`).
#' @param leg_mar Margin values for the kinship legend panel only.
#' A length-4 vector (in `c( bottom, left, top, right )` format that [graphics::par()] 'mar' expects) specifies the full margins, to which `mar_pad` is added.
#' Otherwise, the margins used in the last panel are preserved with the exception that the left margin is set to `mar_pad`, and if `leg_mar` is length-1 (default), it is added to `mar_pad` to specify the right margin.
#' By default the right margin is large enough to accomodate `leg_title` for the given value of `leg_title_line`.
#' @param leg_las The ancestry label orientations (in format that [graphics::mtext()] expects).
#' @param leg_omit If `TRUE`, no legend (second panel) is produced (default `FALSE` is to include legend).
#' @param layout_add If `TRUE` (default) then [graphics::layout()] is called internally to create two panels: the main panel and the color key legend.
#' The original layout is reset when plotting is complete and if `layout_add = TRUE`.
#' If a non-standard layout or additional panels (beyond those provided by this function) are desired, set to `FALSE` and call [graphics::layout()] yourself beforehand.
#' 
#' INDIVIDUAL LABEL OPTIONS
#' 
#' @param names If `TRUE`, the row (individual) names are plotted in the structure barplot.
#' @param names_cex Scaling factor for the individual names.
#' @param names_line Line where individual names are placed.
#' @param names_las Orientation of labels relative to axis.
#' Default (2) makes labels perpendicular to axis.
#'
#' SUBPOPULATION LABEL OPTIONS
#' 
#' @param labs Subpopulation labels for individuals in the admixture matrix.
#' Use a matrix of labels to show groupings at more than one level (for a hierarchy or otherwise).
#' @param labs_cex A vector of label scaling factors for each level of labs.
#' @param labs_las A vector of label orientations (in format that [graphics::mtext()] expects) for each level of labs.
#' @param labs_line A vector of lines where labels are placed (in format that [graphics::mtext()] expects) for each level of labs.
#' @param labs_sep A vector of logicals that specify whether lines separating the subpopulations are drawn for each level of labs.
#' @param labs_lwd A vector of line widths for the lines that divide subpopulations (if `labs_sep = TRUE`) for each level of labs.
#' @param labs_col A vector of colors for the lines that divide subpopulations (if `labs_sep = TRUE`) for each level of labs.
#' @param labs_ticks A vector of logicals that specify whether ticks separating the subpopulations are drawn for each level of labs.
#' @param labs_text A vector of logicals that specify whether the subpopulation labels are shown for each level of labs.
#' Useful for including separating lines or ticks without text.
#' @param labs_even A vector of logicals that specify whether the subpopulations labels are drawn with equal spacing for each level of labs.
#' When `TRUE`, lines mapping the equally-spaced labels to the unequally-spaced subsections of the heatmap are also drawn.
#' @param ... Additional options passed to [graphics::barplot()].
#'
#' @examples
#' # create random proportions for two ancestries
#' Q <- runif( 10 )
#' Q <- cbind( Q, 1 - Q )
#' # add ancestry names
#' colnames( Q ) <- c('A1', 'A2')
#'
#' # plot this data!
#' plot_admix( Q )
#'
#' # See vignette for more elaborate examples!
#' 
#' @export
plot_admix <- function(
                       Q,
                       col = RColorBrewer::brewer.pal( max( ncol( Q ), 3 ), "Paired" ), # diff default vs plot_popkin
                       mar_pad = 0.2, # used for legend mar hack only
                       panel_letters = NA,
                       panel_letters_cex = 1.5,
                       panel_letters_adj = -0.1,
                       axis_lab_cex = 1, # set several labels to the same value at once
                       xlab = "Individuals", # can be NA to omit
                       xlab_line = 1,
                       xlab_cex = axis_lab_cex,
                       ylab = "Ancestry", # can be NA to omit
                       ylab_line = 2, # diff default vs plot_popkin
                       ylab_side = 2,
                       ylab_cex = axis_lab_cex,
                       leg_title = 'Ancestries', # diff default vs plot_popkin
                       leg_title_cex = axis_lab_cex,
                       leg_title_line = 2,
                       leg_cex = 1, # for labels
                       leg_mar = leg_title_line + 1,
                       leg_width = 0.2, # diff default vs plot_popkin
                       leg_las = 0,
                       leg_omit = FALSE,
                       layout_add = !leg_omit,
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
                       ... # passed to [graphics::barplot()]
                       ) {
    # colors can disagree with number of ancestries because RColorBrewer has limitations, fix or die
    K <- ncol( Q )
    if ( K == 1 )
        stop( '`ncol(Q) > 1` is required!' )
    if ( length( col ) > K ) {
        # expected when K=2, because RColorBrewer has a minimum of 3 colors in most cases
        # force match for simplicity
        if ( K == 2 && missing( col ) ) {
            # for default palette it's nicer to pick colors that are more apart than contiguous colors
            # there's at least one more color so this is guaranteed to work
            col <- col[ c(1, 3) ]
        } else {
            # else pick contiguous colors
            col <- col[ 1:K ]
        }
    } else if ( length( col ) < K )
        stop( 'There were fewer colors (', length( col ), ') than ancestries (', K, ')!' )

    # save entire original setup, to reset in the end
    # no.readonly is since some parameters cannot be changed (trying to set them results in ugly warnings)
    par_orig <- graphics::par( no.readonly = TRUE )
    # save original margins, which may get reset per panel (separately of final reset)
    mar_orig <- graphics::par('mar')

    # figure out layout
    # unlike popkin, here we don't combine several plots
    if ( layout_add )
        plot_popkin_layout(
            n = 1, # n_all,
            nr = 1, # layout_rows,
            leg_per_panel = FALSE, # leg_per_panel,
            leg_width = leg_width,
            leg_column = NA # leg_column
        )

    # main plot of admixture proportions (stacked bars)
    # Q gets its columns reordered so they appear from top to bottom, and transposed for [barplot()]
    xc_ind <- graphics::barplot(
                            t( Q ),
                            col = col,
                            xlab = '',
                            xaxt = 'n',
                            xaxs = 'i',
                            yaxt = 'n', # to be able to put axis on other side
                            border = NA,
                            space = 0,
                            legend.text = FALSE, # never add classic (in-panel) barplot legend
                            ...
                        )
    
    # if we want to show labels, use the row/col names to add to picture
    # gap.axis = -1 # added so names show up overlapping rather than not show up at all (can be misleading/confusing)
    if (names)
        graphics::axis( 1, xc_ind, rownames( Q ), las = names_las, cex.axis = names_cex, tick = FALSE, line = names_line, gap.axis = -1 )
    
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
            # these are the barplot-specific changes
            xc_ind = xc_ind,
            doMat = FALSE
        )
    
    # add y-axis
    graphics::axis( ylab_side )
    if ( !is.na( ylab ) )
        graphics::mtext( ylab, side = ylab_side, line = ylab_line, cex = ylab_cex )
    # add x-axis
    if ( !is.na( xlab ) )
        graphics::mtext( xlab, side = 1, line = xlab_line, cex = xlab_cex )

    # add panel letter (in first panel!) if specified
    if ( !is.na( panel_letters ) )
        panel_letter( panel_letters, cex = panel_letters_cex, adj = panel_letters_adj )

    if ( !leg_omit ) {
        ### COLOR LEG admixture ancestry
        # (adds as new panel)
        
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
        
        x <- 1 : length( col )
        # add color key
        graphics::image( y = x, z = rbind( x ), col = col, xaxt = "n", yaxt = "n" ) # , axes=FALSE
        # add labels if available
        # gap.axis = -1 # added so names show up overlapping rather than not show up at all (can be misleading/confusing)  (NOTE: tested in vignette and appears to have not worked, not sure why, but will leave in just in case)
        if ( !is.null( colnames( Q ) ) )
            graphics::axis( 4, at = x, labels = colnames( Q ), tick = FALSE, cex.axis = leg_cex, las = leg_las, gap.axis = -1 )
        # lastly, add axis label
        graphics::mtext( side = 4, leg_title, line = leg_title_line, cex = leg_title_cex )
    }

    # restore original setup when done, but only if we created the default layout
    # otherwise the external layout gets reset, which is bad if we were not done adding panels
    if ( layout_add ) {
        graphics::par( par_orig )
    } else {
        # restore original margins only!
        graphics::par( mar = mar_orig )
    }
}
