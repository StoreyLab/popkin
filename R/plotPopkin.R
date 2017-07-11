## TODO
## - document pub functions!!!
## - add to package!
## - test in vignette!
## DONE
## - figure out dependencies, what needs to be public or not, etc:
##   - plotPopkin: def pub!
##   - heatLegImage: used internally for maps: MapAlleleFreq2.R, humanOrigins12.mapsCombFig.R
##   - panelLetter: unused!!!
##   - printLabs: used for admixture plots too: tgpPopKin03.admixPlot.R, teraStruc03.figs.R
## DEFERRED
## - add letters to panels automatically, actually use this in my fst paper figures!
## - potentially clean up so rangeS is made internally in plotPopkinSingle (to accept rangeReal instead?)

#' Visualize one or more kinship matrices
#'
#' This function plots one or more kinship matrices and a shared legend for the color key.
#' Many options allow for fine control of individual or subpopulation labeling.
#'
#' \code{plotPopkin} plots the input kinship matrices as-is.
#' For best results, a standard kinship matrix (such as the output of \code{\link{popkin}}) should have its diagonal rescaled to contain inbreeding coefficients (\code{\link{inbrDiag}} does this) before \code{plotPopkin} is used.
#'
#' This function permits the labeling of individuals (from row and column names when \code{showNames=TRUE}) and of subpopulations (passed through \code{labs}).
#' The difference is that the labels passed through \code{labs} are assumed to be shared by many individuals, and lines (or other optional visual aids) are added to demarcate these subgroups.
#'
#' For flexibility, this function will work for non-symmetric and even non-square matrices, even though proper kinship matrices are both.
#' For non-symmetric inputs, differing rownames and colnames will display correctly (if \code{showNames==TRUE}).
#' However, numerous options implicitly assume symmetry.
#' For example, only the y-axis is labeled under the assumption that the x-axis is the same.
#' Also, the same subpopulation labels are reproduced on both axes (for clarity).
#'
#' @param x A numeric kinship matrix or a list of matrices
#' @param titles Titles to add to each matrix panel (default is no titles)
#' @param col Colors for heatmap
#' @param xMar Margins for each panel (if a list) or for all panels (if a vector).  Margins are in \code{c(bottom,left,top,right)} format that \code{\link[graphics]{par}('mar')} expects.  By default the existing margin values are used without change
#' @param diagLine If \code{TRUE} adds a line along the diagonal (default no line).  May also be a vector of booleans to set per panel (lengths must agree)
#' @param marPad Margin padding added to all panels (\code{xMar} above and \code{legMar} below).  Default 0.2.  Must be a scalar or a vector of length 4 to match \code{\link[graphics]{par}('mar')}.
#'
#' AXIS LABEL OPTIONS
#' 
#' @param ylab The y-axis label (default "Individuals").  If \code{length(ylab)==1}, the label is placed in the outer margin (shared across panels); otherwise \code{length(ylab)} must equal the number of panels and each label is placed in the inner margin of the respective panel
#' @param ylabAdj The value of "adj" passed to \code{\link[graphics]{mtext}}.  If \code{length(ylab)==1}, only the first value is used, otherwise \code{length(ylabAdj)} must equal the number of panels
#' @param ylabLine The value of "line" passed to \code{\link[graphics]{mtext}}.  If \code{length(ylab)==1}, only the first value is used, otherwise \code{length(ylabLine)} must equal the number of panels
#' 
#' LAYOUT OPTIONS
#' 
#' @param addLayout If \code{TRUE} (default) then \code{\link[graphics]{layout}} is called internally with appropriate values for the required number of panels for each matrix, the desired number of rows (see \code{nr} below) plus the color key legend.  Set to FALSE and call \code{\link[graphics]{layout}} beforehand if a non-standard layout or additional panels (beyond those provided by \code{plotPopkin}) are desired.
#' @param nr Number of rows in layout, used only if \code{addLayout=TRUE}
#'
#' LEGEND (COLOR KEY) OPTIONS
#' 
#' @param legTitle The name of the variable that the heatmap colors measure (default "Kinship")
#' @param legMar Margin vector (in \code{c(bottom,left,top,right)} format that \code{\link[graphics]{par}('mar')} expects) for the legend panel only.  If not provided, the margins used in the last panel are preserved with the exception that the left margin is set to zero (plus the value of \code{marPad}, see above).
#' @param nPretty The desired number of ticks in the legend y-axis (input to \code{\link{pretty}}, see that for more details)
#'
#' INDIVIDUAL LABEL OPTIONS
#' 
#' @param showNames If \code{TRUE}, the column and row names are plotted in the heatmap
#' @param namesCex Scaling factor for the column and row names
#' @param namesLine Line where column and row names are placed
#'
#' SUBPOPULATION LABEL OPTIONS
#' 
#' @param labs Subpopulation labels for individuals.  Use a matrix of labels to show groupings at more than one level (for a hierarchy or otherwise).  If input is a vector or a matrix, the same subpopulation labels are shown for every heatmap panel; the input must be a list of such vectors or matrices if the labels vary per panel
#' @param labsCex A vector of label scaling factors for each level of labs, or a list of such vectors if labels vary per panel
#' @param labsLas A vector of label orientations (in format that \code{\link[graphics]{mtext}} expects) for each level of labs, or a list of such vectors if labels vary per panel
#' @param labsLine A vector of lines where labels are placed (in format that \code{\link[graphics]{mtext}} expects) for each level of labs, or a list of such vectors if labels vary per panel
#' @param labsSkipLines A vector of booleans that specify whether lines separating the subpopulations are drawn for each level of labs, or a list of such vectors if labels vary per panel
#' @param labsLwd A vector of line widths for the lines that divide subpopulations (if \code{labsSkipLines=FALSE}) for each level of labs, or a list of such vectors if labels vary per panel
#' @param labsDoTicks A vector of booleans that specify whether ticks separating the subpopulations are drawn for each level of labs, or a list of such vectors if labels vary per panel
#' @param labsEven A vector of booleans that specify whether the subpopulations labels are drawn with equal spacing for each level of labs, or a list of such vectors if labels vary per panel.  When \code{TRUE}, lines mapping the equally-spaced labels to the unequally-spaced subsections of the heatmap are also drawn
#' 
#' @param ... Additional options passed to \code{\link[graphics]{image}}.  These are shared across panels
#'
#' @examples
#' \dontrun{
#' ## This example assumes input is in BED format and is loaded using BEDMatrix
#' ## "file" is path to BED file (excluding .bed extension)
#' library(BEDMatrix)
#' X <- BEDMatrix(file) # load genotype matrix object
#' Phi <- popkin(X, subpops) # calculate kinship from genotypes and subpopulation labels
#' ## simple plot without labels or any other non-default options
#' ## (see vignette for more elaborate examples)
#' plotPopkin( inbrDiag(Phi) ) # plot the kinship matrix after transforming the diagonal
#' }
#'
#' @export
plotPopkin <- function(x, titles=NULL, col=NULL, xMar=NULL, marPad=0.2, diagLine=FALSE, ylab='Individuals', ylabAdj=NA, ylabLine=0,
                       addLayout=TRUE, nr=1, 
                       legTitle='Kinship', legMar=NULL, nPretty=5,
                       showNames=FALSE, namesCex=1, namesLine=NA,
                       labs=NULL, labsCex=1, labsLas=0, labsLine=0, labsSkipLines=FALSE, labsLwd=1, labsDoTicks=FALSE, labsEven=FALSE,
                       ...) {
    ## wrapper around individual panels and color key
    ## does not set PDF output, margins, layout, etc
    ## assumes x is a list of matrices, if not it is internally turned into one

    if (class(x) != 'list') {
        if (class(x) == 'matrix') {
            x <- list(x) # turn into list so rest works
        } else {
            stop('Fatal: main input is neither a matrix or a list! Class is ', class(x))
        }
    }
    n <- length(x) # number of heatmap panels
    if (!is.null(titles)) {
        if (length(titles) != n)
            stop('Fatal: titles provided are not the same length as data! Data: ', n, ', titles: ', length(titles))
    } else {
        titles <- rep.int('', n) # make blank titles of the same length as data
    }
    ## check label lengths
    if (!is.null(labs) && class(labs) == 'list' && length(labs) != n)
        stop('Fatal: there are ', n, ' panels but ', length(labs), ' label sets!')
    ## expand other things that may vary per panel
    showNames <- repOrDie(showNames, n)
    namesCex <- repOrDie(namesCex, n)
    namesLine <- repOrDie(namesLine, n)
    diagLine <- repOrDie(diagLine, n)
    ## this is for non-scalars per panel, get turned into lists instead
    xMar <- repOrDieList(xMar, n)
    labs <- repOrDieList(labs, n)
    labsCex <- repOrDieList(labsCex, n)
    labsLas <- repOrDieList(labsLas, n)
    labsLine <- repOrDieList(labsLine, n)
    labsLwd <- repOrDieList(labsLwd, n)
    labsSkipLines <- repOrDieList(labsSkipLines, n)
    labsDoTicks <- repOrDieList(labsDoTicks, n)
    labsEven <- repOrDieList(labsEven, n)
    ## ylab behavior is more dynamic!
    if (length(ylab) > 1) {
        ylab <- repOrDie(ylab, n) # just makes sure length is n
        ylabAdj <- repOrDie(ylabAdj, n)
        ylabLine <- repOrDie(ylabLine, n)
    }
    
    ## code needs two versions of the range
    ## - rangeReal is the real range, used in the end so the color key doesn't show values that weren't actually used
    ## - rangeS is a symmetric range, used internally to ensure zero is in the exact middle (set to white in the default)
    ## get range and construct symmetric range that helps with plotting nice colors with white at zero
    rangeReal <- range( unlist(lapply(x, range, na.rm=TRUE)) )  # range of all data plotted
    ## these next few lines force symmetry for colors (might look better, not sure)
    maxS <- max(abs(rangeReal))
    rangeS <- c(-maxS, maxS)

    ## figure out layout given a requested number of rows
    if (addLayout) plotPopkinLayout(n, nr)
    
    marPre <- graphics::par('mar') # save original margins, in case there are changes
    
    ## breaks of all following plots should match!
    for (i in 1:n) {
        if (!is.null(xMar[[i]])) {
            graphics::par(mar=xMar[[i]]+marPad) # change margins if necessary!
        }
        breaks <- plotPopkinSingle(x[[i]], xRange=rangeS, col=col, showNames=showNames[i], namesCex=namesCex[i], namesLine=namesLine[i], labs=labs[[i]], labsCex=labsCex[[i]], labsLas=labsLas[[i]], labsLine=labsLine[[i]], labsLwd=labsLwd[[i]], labsSkipLines=labsSkipLines[[i]], labsDoTicks=labsDoTicks[[i]], labsEven=labsEven[[i]], diagLine=diagLine[i], main=titles[i], ...)
        ## add ylab for every panel when there is more than one choice, and provided it was non-NA
        ## uses inner rather than outer margin (only choice that makes sense)
        if (length(ylab) > 1 && !is.na(ylab[i])) 
            graphics::mtext(ylab[i], side=2, adj=ylabAdj[i], line=ylabLine[i])
    }

    if (!is.null(legMar)) {
        graphics::par(mar=legMar+marPad) # change margins if necessary!
    } else {
        ## change the current left margin to the padding value
        marTmp <- graphics::par('mar') # last margins
        marTmp[2] <- marPad # replace left margin with zero plus pad
        graphics::par(mar=marTmp) # update margins for legend only!
    }
    heatLegImage(breaks, xRange=rangeReal, varname=legTitle, col=col, nPretty=nPretty)

    ## add margin only once if there was only one, place in outer margin (only choice that makes sense)
    if (length(ylab) == 1) {
        graphics::mtext(ylab, side=2, adj=ylabAdj, outer=TRUE, line=ylabLine)
    }
    
    graphics::par(mar=marPre) # restore margins!
}

plotPopkinSingle <- function (x, xRange=range(x, na.rm=TRUE), col=NULL, showNames=FALSE, namesCex=1, namesLine=NA, xlab = "", ylab = "", labs=NULL, labsCex=1, labsLas=0, labsLine=0, labsLwd=1, labsSkipLines=FALSE, labsDoTicks=FALSE, labsEven=FALSE, diagLine=FALSE, ...) {
    ## this "raw" version does not plot legend or set margins, best for optimized scenarios...
    
    ## data validation
    if (is.null(col)) col <- palPopkin() # default coloring
    if (length(col) == 1 && is.character(col)) 
        col <- get(col, mode = "function")
    if (length(dim(x)) != 2)  #  || !is.numeric(x)
        stop("`x' must be a matrix") # numeric
    nr <- nrow(x)
    nc <- ncol(x)
    
    ## figure out breaks for colors
    breaks <- seq(xRange[1], xRange[2], length = length(col) + 1)
    numcols <- length(breaks) - 1
    if (class(col) == "function") 
        col <- col(numcols)
    ## here replace potentially forced xRange with the true range of the data...
    ## but only if data is more extended, don't change anything otherwise!
    rangeRaw <- range(x, na.rm=TRUE)
    if ( breaks[1] > rangeRaw[1] )
        breaks[1] <- rangeRaw[1]
    if ( breaks[numcols+1] < rangeRaw[2] )
        breaks[numcols+1] <- rangeRaw[2]
    
    ## main plot!
    x <- x[nr:1,]
    graphics::image(1:nc, 1:nr, t(x), xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = xlab, ylab = ylab, col = col, breaks = breaks, useRaster=TRUE, ...)
    if (showNames) {
        ## if we want to show labels, use the row/col names to add to picture
        ## labels will be perpendicular to axis (las=2)
        graphics::axis(1, 1:nc, colnames(x), las=2, cex.axis=namesCex, tick = FALSE, line=namesLine)
        graphics::axis(2, 1:nr, rownames(x), las=2, cex.axis=namesCex, tick = FALSE, line=namesLine)
    }
    if (diagLine) {
        graphics::lines(c(1,nr),c(nc,1)) # diagonal line
    }

    ## add subpo labels if present
    if (!is.null(labs)) {
        ## this function figures out multilevel labeling (loops/repeats as needed)
        printLabsMulti(labs, labsCex, labsLas, labsLine, labsLwd, labsSkipLines, labsEven, labsDoTicks)
    }
    
    breaks # return breaks, since legends need it!
}

printLabsMulti <- function(labs, labsCex, labsLas, labsLine, labsLwd, labsSkipLines, labsEven, labsDoTicks) {
    ## normalize so we can loop over cases (assume arbitrary label levels)
    if (class(labs) != 'matrix') {
        labs <- cbind(labs) # a col vector
    }
    n <- ncol(labs)
    if (n > 1) {
        ## expand scalar options as needed
        labsCex <- repOrDie(labsCex, n)
        labsLas <- repOrDie(labsLas, n)
        labsLine <- repOrDie(labsLine, n)
        labsLwd <- repOrDie(labsLwd, n)
        labsSkipLines <- repOrDie(labsSkipLines, n)
        labsEven <- repOrDie(labsEven, n)
        labsDoTicks <- repOrDie(labsDoTicks, n)
    }
    ## loop through now
    for (i in 1:n) {
        printLabs(labs[,i], cex=labsCex[i], las=labsLas[i], line=labsLine[i], lwd=labsLwd[i], skipLines=labsSkipLines[i], even=labsEven[i], doTicks=labsDoTicks[i])
    }
}

repOrDie <- function(vals, n) {
    ## expand scalar options as needed
    if (length(vals) == 1) {
        vals <- rep.int(vals, n)
    } else if (length(vals) != n)
        stop('Error: ', deparse(substitute(vals)), ' does not have the correct length (', length(vals) ,' != ', n, ')')
    vals # return repeated thing or original as needed (guaranteed to be length n)
}

repOrDieList <- function(vals, n) {
    ## expand objects as needed
    if (class(vals) != 'list') {
        vals <- rep(list(vals), n) # this is the desired transformation
    } else if (length(vals) != n)
        stop('Error: ', deparse(substitute(vals)), ' does not have the correct length (', length(vals) ,' != ', n, ')')
    vals # return repeated thing or original as needed (guaranteed to be length n)
}

## Plot color key for heatmap
##
## This plots an image with a single column, showing the relationship between colors and values in the heatmap.
## The image fills the entire panel; use \code{\link[graphics]{layout}} to make sure this panel is scaled correctly, usually so it is taller rather than wider.
##
## This function is provided for users that want greater flexibility in creating plot layouts.
## However, \code{\link{plotPopkin}} will be easier to use and should suffice in most cases, please consider using that before calling this function directly.
## 
## Note \code{\link{plotPopkinSingle}} construct breaks that are symmetric about zero, which ensures that the middle color (white) corresponds to the zero kinship.
## In contrast, xRange need not be symmetric and it is preferably the true range of the data.
##
## @param breaks The vector of \eqn{n+1} values at which colors switch, as returned by \code{\link{plotPopkinSingle}}
## @param varname The name of the variable that the colors measure (i.e. "Kinship")
## @param col Color vector of length \eqn{n}.  Default colors are a progression from blue to white to red obtained from RColorBrewer.
## @param xRange Range of the color key, preferably the range of the data (default is infered from the breaks, but they need not agree)
## @param nPretty The desired number of ticks in the y-axis (input to \code{\link{pretty}}, see that for more details)
##
## @return Nothing
##
## @examples
## \dontrun{
## ## suppose Phi is the only kinship matrix to plot, then...
##
## ## VERSION 1: direct plot through plotPopkin (easiest but restricts layout/margin choices)
## plotPopkin(Phi)
##
## ## VERSION 2: detailed construction that calls plotPopkinSingle and heatLegImage directly, setting layout and other features manually
## ## equals version 1 above but demonstrates the steps that need to be carried out
## rangeReal <- range(Phi) # get the real data range first
## ## also construct symmetric version of this range, so colors work out better...
## maxS <- max(abs(rangeReal))
## rangeS <- c(-maxS, maxS)
## ## start layout
## layout(1:2, widths=c(1,0.1))
## ## plot the heatmap in first panel
## breaks <- plotPopkinSingle(Phi, xRange=rangeS) # pass symmetric range here
## ## add color key to second panel
## heatLegImage(breaks, xRange=rangeReal) # pass real range here
## ## label y-axis on outer margin
## mtext('Individuals', side=2, outer=TRUE)
## }
##
## @export
heatLegImage <- function(breaks, xRange, varname='Kinship', col=NULL, nPretty=5) {
    ## creates a nicer heatmap legend, but it has to be a standalone image (in its own panel, preferably through layout so it's a skinny panel)
    ## this function fills panel, so here we don't set margins/etc (it's best left to the end user)

    if (is.null(col)) col <- palPopkin() # default coloring
    if (!missing(xRange)) {
        ## here's a case where we only want to plot things within an altered range than that of breaks/col
        ## expected application is xRange is real data range, while breaks/col are rigged to be wider (particularly because we forced them to be symmetric, and have a white color at zero)
        ## so we'll proceed by assuming xRange is contained in breaks, but we might want to remove breaks and colors with them

        ## find the bins where our desired range falls
        is <- cut(xRange, breaks, labels=FALSE, include.lowest=TRUE)
        ## toss breaks that we didn't use
        breaks <- breaks[is[1]:(is[2]+1)] # we need to go one over for top
        col <- col[is[1]:is[2]] # colors don't need one over
    }
    ## old processing follows...
    nb <- length(breaks) # length of breaks
    
    ## first plot sequence of colors
    ## NOTE: breaks[2:length(breaks)] plots the sequence of colors because (from ?image):
    ## > intervals are closed on the right and open on the left except for the lowest interval which is closed at both ends.
    ## so plotting all top values of breaks works!
    graphics::image(z = matrix(breaks[2:nb], nrow = 1), col = col, breaks = breaks, xaxt = "n", yaxt = "n")

    ## now we should add axis numbers/ticks
    ## this is a real pain, because 0 and 1 are in the middle of the first and last boxes
    ## old code (originally from heatmap.2) assumed 0 and 1 were at the end of the boxes, not their middles!
    ## given this "delta"
    delta <- 1/(2*(nb-2))
    ## we actually want ends to be at c(-delta, 1+delta)
    lv <- pretty(breaks, n=nPretty)
    xv <- scaleDelta(as.numeric(lv), breaks[1], breaks[nb], delta)
    graphics::axis(4, at = xv, labels = lv)
    
    ## lastly, add axis label
    graphics::mtext(side = 4, varname, line = 2)
}

## internal use only
## update: spreads data to c(0-delta, 1+delta) instead of c(0,1) (which was incorrect)
scaleDelta <- function(x, xMin = min(x), xMax = max(x), delta) {
    ## this is the range we want in the output
    yMin <- -delta
    yMax <- 1 + delta
    ## this is the transformed data to return!
    ## it's easy to verify that (x=xMin -> y=yMin) and (x=xMax -> y=yMax), as desired!
    yMin + (yMax-yMin) * (x - xMin) / (xMax - xMin)
}

plotPopkinLayout <- function(n, nr) {
    ## figure out layout given a requested number of rows
    ## step 1: dimensions
    nr <- as.integer(nr) # in case it's not an integer...
    if (nr < 1) nr <- 1 # just treat as 1 (no complaining)
    if (nr > n) nr <- n # if we asked for more rows than data, set to data (again no complaining)
    nc <- ceiling(n/nr) # this is the correct number of columns (might have blank cells)
    nr <- ceiling(n/nc) # reset backwards in case the nr provided was too large (this will reduce empty space)
    ## step 2: fill layout matrix to just work
    layoutMat <- c(1:n, rep.int(0, nr*nc-n)) # first fill in all nr*nc values, including zeroes as needed
    layoutMat <- matrix(layoutMat, nrow=nr, ncol=nc, byrow=TRUE) # turn into matrix
    layoutMat <- cbind(layoutMat, c(n+1, rep.int(0, nr-1))) # add final column for color key and nothing else
    ## step 3: set up widths vector too
    layoutWidths <- c(rep.int(1, nc), 0.1) # last column for color key is 10% the width of the rest
    ## make layout now!
    graphics::layout(layoutMat, widths=layoutWidths) # note rows are all equal height
}

## make internal function?
palPopkin <- function() {
    ##hmcolBR <- rev(RColorBrewer::brewer.pal(11,"RdBu"))
    hmcolR <- RColorBrewer::brewer.pal(9,"Reds")
    hmcolB <- rev(RColorBrewer::brewer.pal(9,"Blues"))
    ## create a hybrid...
    hmcolBR <- c(hmcolB[1:8], hmcolR)
}

boundaryLabs <- function(labs) {
    ## kinda retarded, but doesn't assume continuity between equal-label cases, could be broadly useful
    ## construct two vectors, one with labels, the other with the boundary points
    n <- length(labs)
    b <- c(1) # first boundary point is always 1
    l <- c(labs[1]) # first label in long list is first label in summary list
    ll <- labs[1] # last label, to detect changes
    ## construct the rest of the vectors
    for (i in 2:n) {
        li <- labs[i]
        ## keep going unless current label doesn't match last label!
        if (li != ll) {
            l <- c(l, li) # add new label
            b <- c(b, i) # add boundary point
            ll <- li # update last label for next iteration
        }
    }
    ## the last boundary is the last point incremented by one (cause that's how we encoded boundaries in the middle)
    b <- c(b, n+1)
    ## return as list
    list(labels=l, boundaries=b)
}

line2user <- function(line, side) {
    ## awesome function from:
    ## http://stackoverflow.com/questions/30765866/get-margin-line-locations-in-log-space/30835971#30835971
    lh <- graphics::par('cin')[2] * graphics::par('cex') * graphics::par('lheight')
    x_off <- diff(graphics::grconvertX(c(0, lh), 'inches', 'npc'))
    y_off <- diff(graphics::grconvertY(c(0, lh), 'inches', 'npc'))
    switch(side,
           `1` = graphics::grconvertY(-line * y_off, 'npc', 'user'),
           `2` = graphics::grconvertX(-line * x_off, 'npc', 'user'),
           `3` = graphics::grconvertY(1 + line * y_off, 'npc', 'user'),
           `4` = graphics::grconvertX(1 + line * x_off, 'npc', 'user'),
           stop("Side must be 1, 2, 3, or 4", call.=FALSE))
}

## unused so far but useful!
panelLetter <- function(letter) {
    graphics::mtext(letter, line=0.5, adj=0, cex=1.5)
}

printLabs <- function(labs, x, doMat=TRUE, cex=1, las=0, lwd=1, skipLines=FALSE, doTicks=FALSE, even=FALSE, line=0) {
    labsObj <- boundaryLabs(labs)
    ## extract the data (smaller var names)
    l <- labsObj$labels
    b <- labsObj$boundaries
    n <- length(l) # number of labels (or number of boundaries minus one)
    m <- max(b) # number of individuals plus 1 (usually n+1, but here I messed up notation above, meh)
    ## for non-barplots, the sensible default is to use the indexes as coordinates (or do they have to be normalized?)
    if (missing(x)) {
        x <- 1:m
    } else {
        ## a hack necessary for barplots at least
        ## in this case x is missing the x[m] element, fill it in!
        x[m] <- x[m-1] + (x[2]-x[1]) # extend by usual gap (assuming it's regular)
    }
    xMin <- min(x)
    xMax <- max(x) # used mostly to reflect coordinates, equals at[n+1]

    ## positions of irregular boundaries
    gapX <- (x[2]-x[1])/2 # shared by ticks and lines
    at <- x[b] - gapX # otherwise things are placed in the middle, shift by one half of a pixel

    ## add ticks using axis()
    if (doTicks) {
        graphics::axis(1, at=at, labels=FALSE, lwd=lwd)
        if (doMat) graphics::axis(2, at=xMax-at, labels=FALSE, lwd=lwd)
    }
    
    if (!skipLines) {
        ## draw black horizontal lines at every boundary (including ends, looks weird otherwise)
        ## assume regular spacing
        graphics::abline(v=at, lwd=lwd)
        if (doMat) graphics::abline(h=xMax-at, lwd=lwd)
    }

    ## label placement, ticky connector line calcs for "even" case
    if (even) {
        ## hard case

        ## place labels equally spaced on x's range...
        ## positions (length of labels)
        gapY <- (xMax - xMin)/n # first divide range into even segments
        y <- xMin + gapY*((1:n) - 0.5) # sequence of middles: length n
        y2 <- xMin + gapY*(0:n) # sequence of boundaries: length n+1
        ## NOTE: y[1] == xMin + gapY/2 # so it starts in middle of first bin
        ## NOTE: y[n] == xMax - gapY/2 # so it ends in the middle of last bin
        ## NOTE: y2[1] == xMin and y2[n+1] == xMax match boundaries exactly

        ## connect label boundaries to irregular boundaries in plot
        ysLines <- line2user(c(0,line), 1) # shared by every label boundary on x-axis
        xsLines <- line2user(c(0,line), 2) # shared by every label boundary on y-axis
        for (i in 1:(n+1)) {
            ## middle coordinate is a bit messy: first get ith and ith+1 boundaries, then turn them to coordinates, then average
            ## don't forget that end is first element of next group, so we always need to reduce it by one
            xi <- at[i] # coincides with tick positions
            yi <- y2[i] # boundary of words
            graphics::lines(c(xi, yi), ysLines, lwd=lwd, xpd=NA) # draw line for x-axis
            if (doMat) graphics::lines(xsLines, xMax - c(xi, yi), lwd=lwd, xpd=NA)
        }
        
    } else {
        ## put labels in middle of each range
        y <- (at[1:n] + at[2:(n+1)])/2 # sequence of label positions
    }

    ## place labels at "y" (set according to "even=TRUE" or otherwise)
    graphics::mtext(l, side=1, at=y, cex=cex, las=las, line=line)
    ## for matrices, do both ways!
    if (doMat) graphics::mtext(l, side=2, at=xMax-y, cex=cex, las=las, line=line)
}
