# to be run when popkin definitely works, to reconstruct test that to compare to when there are changes!

library(popkin)
library(BEDMatrix)

# NOTE: small BED data came from ~/dbs/popgenSim/psd1Ds2/
# cp ~/dbs/popgenSim/psd1Ds2/Xs.{bed,bim,fam} .

# NOTE: this X has missingness! (super important to test for)
X <- BEDMatrix('Xs') # works if script is run from local dir...
X <- t(X[,]) # immediately convert into a regular R int matrix
# data has 10 individuals
stopifnot(ncol(X) == 10) # make sure in case it changes
subpops0 <- 1:ncol(X) # trivial case
subpops <- c(1,1,2,2,3,3,4,4,5,5) # group by pairs, still reasonable!

# generate weights for tests
w0 <- weights_subpops(subpops0)
w <- weights_subpops(subpops)

# high-level estimation (what we expect most users to run)
Phi0 <- popkin(X)
Phi <- popkin(X, subpops)

# make all intermediate data to compare to
obj <- popkin_A(X)
A <- obj$A
M <- obj$M
Amin0 <- popkin_A_min_subpops(A) # minimum value
Amin <- popkin_A_min_subpops(A, subpops) # should also be minimum value
phiMin0 <- popkin_A_min_subpops(Phi0, subpops) # redo to bad kinship estimate

# implied Fst
fst <- fst(Phi)
fstW <- fst(Phi, w)
# inbreeding vectors
inbr <- inbr(Phi)
# pairwise Fst
pwF <- pwfst(Phi)
# and kinship matrices with inbreeding along diagonal
PhiInbr <- inbr_diag(Phi)

# save these values, to use in testing later!
save(X, subpops, subpops0, Phi, Phi0, A, M, Amin, Amin0, phiMin0, w, w0, fst, fstW, inbr, pwF, PhiInbr, file = 'Xs.RData')

# wanted to compare non-NA counts to GCTA:
# ~/bin/gcta_1.93.2beta/gcta64 --bfile Xs --make-grm --out Xs
