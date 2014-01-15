#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

using namespace arma; 
using namespace Rcpp;

List tf(const colvec& y, const colvec& l, const mat& D, const colvec& P, 
        const int& n, const int& k, const double& e, const double& tau){
  // author Edward A. Roualdes
  // TODO: use sp_mat
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

  // P.print("P:");
  // D.row(1).print("D[1,]:");
    

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
    bOut.col(j) = b;
    iterOut(j) = iter;
  }

  return List::create(Named("coefficients") = bOut,
                      Named("lambda") = l,
                      Named("k") = k,
                      Named("iters") = iterOut);
}
