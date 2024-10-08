% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_popkin.R
\name{plot_popkin}
\alias{plot_popkin}
\title{Visualize one or more kinship matrices and other related objects}
\usage{
plot_popkin(
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
  panel_letters_adj = -0.1,
  ylab = "Individuals",
  ylab_adj = NA,
  ylab_line = 0,
  ylab_side = 2,
  ylab_per_panel = FALSE,
  layout_add = TRUE,
  layout_rows = 1,
  leg_per_panel = FALSE,
  leg_title = "Kinship",
  leg_cex = 1,
  leg_n = 5,
  leg_mar = 3,
  leg_width = 0.3,
  leg_column = NA,
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
  labs_col = "black",
  labs_ticks = FALSE,
  labs_text = TRUE,
  labs_even = FALSE,
  labs_even_line = labs_line,
  null_panel_data = FALSE,
  weights = NULL,
  raster = is.null(weights),
  sym = FALSE,
  ...
)
}
\arguments{
\item{kinship}{A numeric kinship matrix, a \code{phylo} or \code{function} object, or a list of any such objects (at least one kinship matrix is expected).
This list may contain \code{NULL} elements (makes blank panels with titles; good for placeholders or non-existent data)
\code{phylo} objects are plotted with \code{\link[=plot_phylo]{plot_phylo()}}, which is a wrapper around \code{\link[ape:plot.phylo]{ape::plot.phylo()}} that makes some adjustments so resulting plots are more consistent with accompanying kinship matrices.
\code{function} elements are executed without arguments, and are expected to produce single arbitrary plot panels.}

\item{titles}{Titles to add to each matrix panel (default is no titles).
Applied to kinship and phylo panels only.}

\item{col}{Colors for kinship heatmap (default is a red-white-blue palette symmetric about zero constructed using \code{RColorBrewer}).}

\item{col_n}{The number of colors to use in the heatmap (applies if \code{col = NULL}).}

\item{mar}{Margins shared by all panels (if a vector) or for each panel (if a list of such vectors).
If the vector has length 1, \code{mar} corresponds to the shared lower and left margins, while the top and right margins are set to zero.
If this length is 2, \code{mar[1]} is the same as above, while \code{mar[2]} is the top margin.
If this length is 4, then \code{mar} is a fully-specified margin vector in the standard format \code{c(bottom, left, top, right)} that \code{\link[graphics:par]{graphics::par()}} 'mar' expects.
Vectors of invalid lengths produce a warning.
Note the padding \code{mar_pad} below is added to every margin if set.
If \code{NULL}, the original margin values are used without change, and are reset for every panel that has a \code{NULL} value.
The original margins are also reset after plotting is complete.
Applied to panels of all types (kinship, phylo, and function).}

\item{mar_pad}{Margin padding added to all panels (\code{mar} above and \code{leg_mar} below).
Default 0.2.
Must be a scalar or a vector of length 4 to match \code{\link[graphics:par]{graphics::par()}} 'mar'.
Applied to panels of all types (kinship, phylo, and function).}

\item{oma}{Outer margin vector.
If length 1, the value of \code{oma} is applied to the left outer margin only (so \code{ylab} below displays correctly) and zero outer margins elsewhere.
If length 4, all outer margins are expected in standard format \code{\link[graphics:par]{graphics::par()}} 'mar' expects (see \code{mar} above).
\code{mar_pad} above is never added to outer margins.
If \code{NULL}, no outer margins are set (previous settings are preserved).
Vectors of invalid lengths produce a warning.
Note: if \code{layout_add = FALSE}, this function still (re)sets the outer margins if \code{oma} is not \code{NULL}, which can be convenient if \code{plot_popkin} generates the first few panels, but otherwise a partial multipanel figure will be reset unless \code{oma = NULL} is also set!}

\item{diag_line}{If \code{TRUE} adds a line along the kinship diagonal (default no line).
May also be a vector of logicals to set per panel (lengths must agree).
Has no effect on non-kinship panels.}

\item{panel_letters}{Vector of strings for labeling panels (default A-Z).
No labels are added if \code{NULL}, or when there is only one panel except if its set to a single letter in that case (this behavior is useful if goal is to have multiple external panels but popkin only creates one of these panels).
Applied to panels of all types (kinship, phylo, and function).}

\item{panel_letters_cex}{Scaling factor of panel letters (default 1.5).}

\item{panel_letters_adj}{X-axis adjustment for panel letters (default -0.1).
Negative values place the letter into the left margin area.
Might need adjustment depending on the size of the left margin.}

\item{ylab}{The y-axis label (default "Individuals").
If \code{length(ylab) == 1}, the label is placed in the outer margin (shared across panels);
otherwise \code{length(ylab)} must equal the number of panels and each label is placed in the inner margin of the respective panel.
Applied to panels of all types (kinship, phylo, and function).}

\item{ylab_adj}{The value of \code{adj} for \code{ylab} passed to \code{\link[graphics:mtext]{graphics::mtext()}}.
If \code{length(ylab) == 1}, only the first value is used, otherwise \code{length(ylab_adj)} must equal the number of panels.}

\item{ylab_line}{The value of \code{line} for \code{ylab} passed to \code{\link[graphics:mtext]{graphics::mtext()}}.
If \code{length(ylab) == 1}, only the first value is used, otherwise \code{length(ylab_line)} must equal the number of panels.}

\item{ylab_side}{The value of \code{side} for \code{ylab} passed to \code{\link[graphics:mtext]{graphics::mtext()}} (2 is y-axis, 1 is x-axis, can also place on top (3) or right (4)).
If \code{length(ylab) == 1}, only the first value is used, otherwise \code{length(ylab_side)} must equal the number of panels.}

\item{ylab_per_panel}{Forces y-axis labels to appear for each panel, in the inner margins.
Most useful to cover the case where there is a single panel but no outer margins (\code{oma = NULL}).

LAYOUT OPTIONS}

\item{layout_add}{If \code{TRUE} (default) then \code{\link[graphics:layout]{graphics::layout()}} is called internally with appropriate values for the required number of panels for each matrix, the desired number of rows (see \code{layout_rows} below) plus the color key legend.
The original layout is reset when plotting is complete if \code{layout_add = TRUE}.
If a non-standard layout or additional panels (beyond those provided by \code{plot_popkin}) are desired, set to \code{FALSE} and call \code{\link[graphics:layout]{graphics::layout()}} yourself beforehand; in this case you may want to set \code{oma = NULL} (above) as well!}

\item{layout_rows}{Number of rows in layout, used only if \code{layout_add = TRUE}.

LEGEND (COLOR KEY) OPTIONS}

\item{leg_per_panel}{If \code{TRUE}, every kinship matrix get its own legend/color key (best for matrices with very different scales), and each phylo tree has its own x-axis range.
If \code{FALSE} (default), a single legend/color key is shared by all kinship matrix panels, and also every tree has the same x-axis range (different from the kinship range).}

\item{leg_title}{The name of the variable that the kinship heatmap colors measure (default "Kinship"), or a vector of such values if they vary per panel.}

\item{leg_cex}{Scaling factor for \code{leg_title} (default 1), or a vector of such values if they vary per panel.}

\item{leg_n}{The desired number of ticks in the kinship legend y-axis, and phylo x-axis (input to \code{\link[=pretty]{pretty()}}, see that for more details), or a vector of such values if they vary per panel.}

\item{leg_mar}{Margin values for the kinship legend panel only, or a list of such values if they vary per panel.
A length-4 vector (in \code{c( bottom, left, top, right )} format that \code{\link[graphics:par]{graphics::par()}} 'mar' expects) specifies the full margins, to which \code{mar_pad} is added.
Otherwise, the margins used in the last panel are preserved with the exception that the left margin is set to \code{mar_pad}, and if \code{leg_mar} is length-1, it is added to \code{mar_pad} to specify the right margin.}

\item{leg_width}{The width of the legend panel, relative to the width of a single main panel.
This value is passed to \code{\link[graphics:layout]{graphics::layout()}} (ignored if \code{layout_add = FALSE}).}

\item{leg_column}{The column number in which to place the kinship legend (default \code{NA} is for last column).
Ignored if \code{leg_per_panel = TRUE}.

INDIVIDUAL LABEL OPTIONS}

\item{names}{If \code{TRUE}, the column and row names are plotted in the kinship heatmap, or a vector of such values if they vary per panel.
(\code{names} has no effect on \code{phylo} panels, whose tip labels are always plotted, or other panel types.)}

\item{names_cex}{Scaling factor for the column and row names of a kinship matrix, or the tip labels of a \code{phylo} object, or a vector of such values if they vary per panel.}

\item{names_line}{Line where kinship column and row names are placed, or a vector of such values if they vary per panel.
Has no effect on non-kinship panels.}

\item{names_las}{Orientation of labels relative to axis.
Default (2) makes labels perpendicular to axis.
Has no effect on non-kinship panels.

SUBPOPULATION LABEL OPTIONS}

\item{labs}{Subpopulation labels for individuals in kinship matrices.
Use a matrix of labels (individuals along rows, levels along columns) to show groupings at more than one level (for a hierarchy or otherwise).
If input is a vector or a matrix, the same subpopulation labels are shown for every kinship matrix; the input must be a list of such vectors or matrices if the labels vary per panel.
Has no effect on non-kinship panels.}

\item{labs_cex}{A vector of label scaling factors for each level of labs, or a list of such vectors if labels vary per panel.}

\item{labs_las}{A vector of label orientations (in format that \code{\link[graphics:mtext]{graphics::mtext()}} expects) for each level of labs, or a list of such vectors if labels vary per panel.}

\item{labs_line}{A vector of lines where labels are placed (in format that \code{\link[graphics:mtext]{graphics::mtext()}} expects) for each level of labs, or a list of such vectors if labels vary per panel.}

\item{labs_sep}{A vector of logicals that specify whether lines separating the subpopulations are drawn for each level of labs, or a list of such vectors if labels vary per panel.}

\item{labs_lwd}{A vector of line widths for the lines that divide subpopulations (if \code{labs_sep = TRUE}) for each level of labs, or a list of such vectors if labels vary per panel.}

\item{labs_col}{A vector of colors for the lines that divide subpopulations (if \code{labs_sep = TRUE}) for each level of labs, or a list of such vectors if labels vary per panel.}

\item{labs_ticks}{A vector of logicals that specify whether ticks separating the subpopulations are drawn for each level of labs, or a list of such vectors if labels vary per panel.}

\item{labs_text}{A vector of logicals that specify whether the subpopulation labels are shown for each level of labs, or a list of such vectors if labels vary per panel.
Useful for including separating lines or ticks without text.}

\item{labs_even}{A vector of logicals that specify whether the subpopulations labels are drawn with equal spacing for each level of labs, or a list of such vectors if labels vary per panel.
When \code{TRUE}, lines mapping the equally-spaced labels to the unequally-spaced subsections of the heatmap are also drawn.}

\item{labs_even_line}{A vector of lines where the ends of "even" lines are placed (in format that \code{\link[graphics:mtext]{graphics::mtext()}} expects) for each level of labs, or a list of such vectors if labels vary per panel.
Ignored unless \code{labs_even = TRUE}.}

\item{null_panel_data}{If \code{FALSE} (default), panels with \code{NULL} kinship matrices must not have titles or other parameters set, and no panel letters are used in these cases.
If \code{TRUE}, panels with \code{NULL} kinship matrices must have titles and other parameters set.
In the latter case, these \code{NULL} panels also get panel letters.
The difference is important when checking that lengths of non-singleton parameters agree.}

\item{weights}{A vector with weights for every individual, or a list of such vectors if they vary per panel.
The width of every individual in the kinship matrix becomes proportional to their weight.
Individuals with zero or negative weights are omitted.
Has no effect on non-kinship panels.}

\item{raster}{A logical equivalent to \code{useRaster} option in the \code{image} function used internally, or a vector of such logicals if the choice varies per panel.
If \code{weights} are non-\code{NULL} in a given panel, \code{raster = FALSE} is forced (this is necessary to plot images where columns and rows have variable width).
If \code{weights} are \code{NULL}, the default is \code{raster = TRUE}, but in this case the user may override (for example, so panels are visually coherent when some use weights while others do not, as there are small differences in rendering implementation for each value of \code{raster}).
Note that a multipanel figure with a list of \code{weights} sets \code{raster = FALSE} to all panels by default, even if the weights were only applied to a subset of panels.
Has no effect on non-kinship panels.}

\item{sym}{If \code{FALSE} (default), plots non-symmetric (but square) kinship matrices without issues.  If \code{TRUE}, stops if any input kinship matrices (excluding \code{phylo} or \code{function} objects) are not symmetric.

AXIS LABEL OPTIONS}

\item{...}{Additional options passed to \code{\link[graphics:image]{graphics::image()}}.
These are shared across kinship panels.
Have no effect on non-kinship panels.}
}
\description{
This function plots one or more kinship matrices, trees (class \code{phylo} objects, see \code{ape} package), and arbitrary functions, and a shared legend for the kinship color key.
Many options allow for fine control of individual or subpopulation labeling.
}
\details{
\code{plot_popkin} plots the input kinship matrices as-is.
For best results, a standard kinship matrix (such as the output of \code{\link[=popkin]{popkin()}}) should have its diagonal rescaled to contain inbreeding coefficients using \code{\link[=inbr_diag]{inbr_diag()}} before \code{plot_popkin} is used.

This function permits the labeling of individuals (from row and column names when \code{names = TRUE}) and of subpopulations (passed through \code{labs}).
The difference is that the labels passed through \code{labs} are assumed to be shared by many individuals, and lines (or other optional visual aids) are added to demarcate these subgroups.
}
\examples{
# Construct toy data
X <- matrix(c(0,1,2,1,0,1,1,0,2), nrow = 3, byrow = TRUE) # genotype matrix
subpops <- c(1,1,2) # subpopulation assignments for individuals

# NOTE: for BED-formatted input, use BEDMatrix!
# "file" is path to BED file (excluding .bed extension)
## library(BEDMatrix)
## X <- BEDMatrix(file) # load genotype matrix object

# estimate the kinship matrix from the genotypes "X"!
kinship <- popkin(X, subpops) # calculate kinship from X and optional subpop labels

# simple plot of the kinship matrix, marking the subpopulations only
# note inbr_diag replaces the diagonal of kinship with inbreeding coefficients
# (see vignette for more elaborate examples)
plot_popkin( inbr_diag(kinship), labs = subpops )

}
