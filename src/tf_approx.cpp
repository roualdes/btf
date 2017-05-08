#include <RcppEigen.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

template <typename T>
Eigen::SparseMatrix<double> mkDiag(const Eigen::DenseBase<T>& val) {
  int I = val.size();
  Eigen::SparseMatrix<double> W(I,I); W.reserve(I);
  for (int i=0; i<I; i++) {
    W.insert(i,i) = val(i);    // insert along diagonal
  }
  return W;
}

// [[Rcpp::export]]
List tf_approx(const Eigen::Map<Eigen::VectorXd>& y,
               const Eigen::Map<Eigen::VectorXd>& l,
               const Eigen::MappedSparseMatrix<double>& D,
               const int& k, const double& eps,
               const double& tau, const int& max_iter) {
  // author Edward A. Roualdes
  // initialize some values
  const int n = y.size();       // sample size
  const int nk = n-k-1;
  const int J = l.size();        // number of lambdas
  Eigen::VectorXd beta = y;                 // initial values
  Eigen::VectorXd beta_old(n);

  // create output containers
  Eigen::MatrixXd beta_out(n, J);        // beta out matrix
  Eigen::SparseMatrix<double> W(nk, nk);
  W = mkDiag(1.0/((D*beta).cwiseAbs().array() + eps));
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > LLt; // linear system
  LLt.analyzePattern(D.transpose()*W*D); // symbolic decomposition on the sparcity
  LLt.setShift(1.0, 1.0);                // add Identity
  Eigen::VectorXd iter_out(J);              // store number iterations per lambda

  // for each lambda
  for (int j=0; j<J; j++) {
    int iter = 0;               // count iterations until convergence

    // while solution not close enough && under max iterations
    while (iter < max_iter) {
      beta_old = beta;          // store old estimates
      W = mkDiag(1.0/((D*beta).cwiseAbs().array() + eps));

      // create linear system and solve
      LLt.factorize(l(j)*D.transpose()*W*D); // computational decomp on symbolic
      beta = LLt.solve(y);

      // check convergence
      if (beta.isApprox(beta_old, tau)) break;
      ++iter;
    }
    // store estimates for each lambda
    beta_out.col(j) = beta;
    iter_out(j) = iter;
  }

  return List::create(Named("coefficients") = beta_out,
                      Named("iters") = iter_out);
}
