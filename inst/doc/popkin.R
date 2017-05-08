## ------------------------------------------------------------------------
library(popkin)
library(lfa) # for hgdp_subset sample data only
X <- hgdp_subset # rename for simplicity
dim(X)

## ------------------------------------------------------------------------
# shorten subpopulation labels
colnames(X)[colnames(X)=='AFRICA'] <- 'AFR'
colnames(X)[colnames(X)=='MIDDLE_EAST'] <- 'MDE'
colnames(X)[colnames(X)=='EUROPE'] <- 'EUR'
colnames(X)[colnames(X)=='CENTRAL_SOUTH_ASIA'] <- 'SAS'
colnames(X)[colnames(X)=='EAST_ASIA'] <- 'EAS'
colnames(X)[colnames(X)=='OCEANIA'] <- 'OCE'
colnames(X)[colnames(X)=='AMERICA'] <- 'AMR'
# order roughly by distance from Africa
popOrder <- c('AFR', 'MDE', 'EUR', 'SAS', 'EAS', 'OCE', 'AMR')
# applies reordering
X <- X[,order(match(colnames(X), popOrder))]
subpops <- colnames(X) # extract subpopulations vector

## ------------------------------------------------------------------------
Phi <- popkin(X, subpops)

## ------------------------------------------------------------------------
# quick and dirty heatmap of the kinship matrix
library(RColorBrewer)
myHeatmap <- function(Phi, subpops=NULL) {
    # place subpop labels in middle of their ranges
    subpopLabs <- NULL
    if (!is.null(subpops)) {
        # mean index per subpop
       	subpopLabs <- aggregate(1:length(subpops), list(subpop=subpops), mean)
    }
    colHM <- brewer.pal(9, 'Reds') # heatmap colors
    par(xaxt='n', yaxt='n') # heatmap doesn't omit axes correctly, force here
    # plot as image, without reordering, dendrograms, etc
    heatmap(Phi, Rowv=NA, Colv=NA, symm=TRUE, col=colHM,
        xlab='individuals', ylab='individuals',
        add.expr={
            # add population labels
	    if (!is.null(subpopLabs)) {
	        mtext(subpopLabs$subpop, 1, at=subpopLabs$x, line=1, cex=0.7)
	        mtext(subpopLabs$subpop, 4, at=subpopLabs$x, line=1, cex=0.7)
	    }
        })
    par(xaxt='s', yaxt='s') # reset to other plots aren't affected
    
}

## ---- fig.width=6, fig.height=6, fig.align='center'----------------------
myHeatmap(Phi, subpops)

## ---- fig.width=6, fig.height=6, fig.align='center'----------------------
myHeatmap(inbrDiag(Phi), subpops)

## ------------------------------------------------------------------------
# get weights
w <- weightsSubpops(subpops)
# compute FST!
# Note: don't use the output to inbrDiag(Phi) or FST will be wrong!
fst(Phi, w)

## ---- fig.width=6, fig.height=6, fig.align='center'----------------------
pwF <- pwfst(Phi) # compute pairwise FST matrix from kinship matrix
myHeatmap(pwF, subpops) # NOTE no need for inbrDiag() here

## ---- fig.width=6, fig.height=3, fig.align='center'----------------------
inbrs <- inbr(Phi) # vector of inbreeding coefficients
par(mar=c(4, 4, 0, 0) + 0.1) # reduce margins
plot(density(inbrs), xlab='inbreeding coefficient', main='') # see their distribution

## ---- fig.width=6, fig.height=6, fig.align='center'----------------------
# filter to only keep individuals within AFR
indexesAfr <- subpops == 'AFR'
PhiAfr <- Phi[indexesAfr,indexesAfr]

# estimate FST before rescaling (this value will be wrong, too high!)
fst(PhiAfr)

# now rescale
# since subpops is missing, minimum Phi value is set to zero
# (no averaging between subpopulations)
PhiAfr <- rescalePopkin(PhiAfr)
# FST is now correct, relative to the MRCA of AFR individuals
fst(PhiAfr)
# kinship matrix visualization
myHeatmap(inbrDiag(PhiAfr))

