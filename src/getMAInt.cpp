// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
Rcpp::List getMAInt (Eigen::Map<Eigen::MatrixXi> & X){
  
  // get dimensions of input matrix
  const int n = X.cols();
  const int m = X.rows();
  int mo = 0; // will count fixed SNPs, to offset M later...
  
  // initialize the desired matrix M with all values m
  Eigen::MatrixXi M = Eigen::MatrixXi::Constant(n, n, m);
  //  std::cout << M << std::endl; // yey looks good!
  
  // initialize per-individual adjustments
  Eigen::VectorXi mj = Eigen::VectorXi::Zero(n);
  Eigen::VectorXi uj = Eigen::VectorXi::Ones(n); // unit vector
  Eigen::VectorXi js = Eigen::VectorXi::Zero(n); // similar, list of NA indexes...
  int nj; // number of NAs, or js slots filled (per SNP)
  int n0; // indicate other cases, to know when a SNP is fixed (no need to count precisely, just 0 vs >0)
  int n1;
  int n2;
  //  std::cout << mj << std::endl; // yey looks good!
  
  // indexes for navigation
  int i; // per SNP
  int j; // per individuals j,k
  int k;
  int b; // indexes of NA lists (j=js[b]; k=js[c];)
  int c;
  for (i=0; i < m; i++) {
    // count things on first pass, turn NAs to 1 but don't update M or mj yet!
    nj=0; // reset at each SNP
    n0=0; // reset at each SNP
    n1=0; // reset at each SNP
    n2=0; // reset at each SNP
    for (j=0; j < n; j++) {
      if ( X(i,j) == NA_INTEGER ) { // this is how we know we have NAs in C++...
	js(nj) = j; // add to list of NA indexes in this SNP
	nj++; // increment in case we get more
	X(i,j) = 1; // set this NA to 1 now! (because we havent't centered yet!)
      } else if (n1 == 0) { // this indicates we have to keep checking
	if (X(i,j)==1) {
	  n1 = 1; // set as heterozygote, definitely polymorphic
	} else if (X(i,j)==0) {
	  n0 = 1; // set this homozygote
	  if (n2 > 0) {
	    n1 = 1; // if both homozygotes have been observed, turn "on" n1 (the overall indicator!)
	  }
	} else { // assumes X(i,j) == 2, no other cases...
	  n2 = 1;
	  if (n0 > 0) {
	    n1 = 1; // ditto both homozygotes were observed, so n1 is "on"
	  }
	}
      }
    }
    if (n1 == 0) {
      // if this is true, SNP was fixed!
      mo++; // increment fixed SNP count
      // turn entire row into 1's (after centering later, these will be zeroes and will be ignored!)
      for (j=0; j < n; j++) {
	X(i,j) = 1;
      }
    } else {
      // continue for polymorphic SNPs
      if (nj > 0) {
	// process NA cases, updating M and mj as needed!
	for (b=0; b < nj; b++) {
	  j = js[b];
	  mj(j)++; // increment number of NAs for this person
	  M(j,j)++; // increment self pair outside loop (or changes are doubled)
	  for (c=0; c<b; c++) { // exclude self case!!! (c!=b)
	    k = js[c];
	    M(j,k)++; // increment these pairs (since we said matrix is symmetric, one way is enough?)
	    M(k,j)++; // have to do both ways it turns out...
	  }
	}
      }
    }
  }

  // now add mj adjustments!  This trick does it!
  M = M.selfadjointView<Eigen::Lower>().rankUpdate(mj, uj, -1);
  M = M.array() - mo; // remove fixed SNPs from overall count!

  X = X.array() - 1;
  
  // second part is fairly trivial now (though it's the slowest, comparatively)
  Eigen::MatrixXi SA = Eigen::MatrixXi(n, n).setZero().selfadjointView<Eigen::Lower>().rankUpdate(X.adjoint());

  return Rcpp::List::create(Rcpp::Named("M") = M, Rcpp::Named("SA") = SA);
}
