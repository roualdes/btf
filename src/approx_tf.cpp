#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

using namespace arma; 
using namespace Rcpp;

List approx_tf(const colvec& y, const colvec& l, const mat& D, const int& n, 
               const int& k, const double& e, const double& tau){

  // author Edward A. Roualdes
  // TODO: use sp_mat
  // initialize some values
  const int J = l.n_rows;       // number of lambdas
  const mat I = eye(n,n);       // make identity
  const int nk = n-k-1;         // dim phi
  const int maxiter = 10000;
  mat bOut(n,J, fill::zeros);   // store output
  colvec iterOut(J);            // store iterations for each lambda
  colvec b = y;                 // initialize beta
  mat M(n,n, fill::zeros);      // initialize linear system
  mat phi(nk,nk, fill::zeros);  // initialize phi
  colvec bOld(n);

  // algorithm
  for (int j=0; j<J; j++) {
    int iter = 0;
    while (iter < maxiter) {
      bOld  = b;                     // store old estimates
      phi.diag() = 1.0/(abs(D*b)+e); // build phi

      // create linear system and solve
      M = I + l(j)*D.t()*phi*D;
      b = solve(M, y);        

      // check convergence
      if (norm(b - bOld, "Inf") < tau) break;
      ++iter;
    }
    // store estimates for each lambda
    bOut.col(j) = b;
    iterOut(j) = iter;          // number iterations at convergence
  }

  return List::create(Named("coefficients") = bOut,
                      Named("lambda") = l,
                      Named("k") = k,
                      Named("iters") = iterOut);
}
