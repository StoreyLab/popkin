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

## ---- fig.width=6, fig.height=5, fig.align='center'----------------------
# set outer margin for axis labels (left and right are non-zero)
par(oma=c(0,1.5,0,3))
# set inner margin for subpopulation labels (bottom and left are non-zero), add padding
par(mar=c(1,1,0,0)+0.2)
# now plot!
plotPopkin(Phi, labs=subpops)

## ---- fig.width=6, fig.height=5, fig.align='center'----------------------
par(oma=c(0,1.5,0,3))
par(mar=c(1,1,0,0)+0.2)
plotPopkin(inbrDiag(Phi), labs=subpops)

## ---- fig.width=6, fig.height=5, fig.align='center'----------------------
par(oma=c(0,1.5,0,3))
# increase margins because labels go farther out
par(mar=c(2,2,0,0)+0.2)
plotPopkin(inbrDiag(Phi), labs=subpops, labsEven=TRUE, labsLine=1)

## ------------------------------------------------------------------------
# get weights
w <- weightsSubpops(subpops)
# compute FST!
# Note: don't use the output to inbrDiag(Phi) or FST will be wrong!
fst(Phi, w)

## ---- fig.width=6, fig.height=5, fig.align='center'----------------------
pwF <- pwfst(Phi) # compute pairwise FST matrix from kinship matrix
legTitle <- expression(bold(paste('Pairwise ', F[ST]))) # fancy legend label
par(oma=c(0,1.5,0,3))
par(mar=c(2,2,0,0)+0.2)
# NOTE no need for inbrDiag() here!
plotPopkin(pwF, labs=subpops, labsEven=TRUE, labsLine=1, legTitle=legTitle)

## ---- fig.width=4, fig.height=2, fig.align='center'----------------------
inbrs <- inbr(Phi) # vector of inbreeding coefficients
par(mar=c(4, 4, 0, 0) + 0.1) # reduce margins
plot(density(inbrs), xlab='inbreeding coefficient', main='') # see their distribution

## ---- fig.width=3, fig.height=2, fig.align='center'----------------------
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
par(oma=c(0,1.5,0,3))
# use zero margins because there are no subpopulation labels
par(mar=c(0,0,0,0)+0.2)
plotPopkin(inbrDiag(PhiAfr))

## ---- fig.width=6, fig.height=2.5, fig.align='center'--------------------
par(oma=c(0,1.5,0,3))
par(mar=c(2,2,0,0)+0.2)
# transform both matrices and store in a list
x <- list(inbrDiag(Phi), inbrDiag(PhiAfr))
# dummy labels to have lines in second panel
subpopsAfr <- subpops[indexesAfr]
# pass labels that differ per panel using a list
labs <- list(subpops, subpopsAfr)
# labsEven and labsLine passed as scalars are shared across panels
# labsCex reduces the label size, since there's less space here!
plotPopkin(x, labs=labs, labsEven=TRUE, labsLine=1, labsCex=0.5)

