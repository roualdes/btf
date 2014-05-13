#include "individual.cpp"

// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::export]]
Eigen::MatrixXd dexpBTF(const int& iter,
                        const Eigen::Map<Eigen::VectorXd>& y, 
                        const int& k, const Eigen::MappedSparseMatrix<double>& D, 
                        const double& alpha, const double& rho) {
  
  individual *btf;
  btf = new individual(y, k, D, alpha, rho); // initialize btf object
  int P = btf->n + btf->nk + 2; // count number of parameters
  Eigen::MatrixXd history(iter, P); // initialize matrix of posterior draws

  for (int i=0; i<iter; ++i) {
    btf->upBeta();
    btf->upS2();
    btf->upOmega2Lambda2();
    history.row(i) << btf->beta.transpose(), btf->s2, btf->l2, btf->o2.transpose();
  }
  return history;
}
