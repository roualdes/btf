#include "individual.cpp"

// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::export]]
Eigen::MatrixXd dexpBTF(const int& iter,
                        const Eigen::Map<Eigen::VectorXd>& y, 
                        const int& k, const Eigen::MappedSparseMatrix<double>& D, 
                        const double& alpha, const double& rho) {
  
  // initialize btf object
  individual *btf;
  btf = new individual(y, k, D, alpha, rho); 

  // initialize matrix of posterior draws
  int P = btf->n + btf->nk + 3; // count number of parameters
  Eigen::MatrixXd history = Eigen::MatrixXd::Zero(iter, P); 

  // runs sampler
  for (int i=0; i<iter; ++i) {
    btf->upBeta(); btf->upS2();
    btf->upLambda2(); btf->upOmega2(); // btf->upA(); // btf->upR();
    history.row(i) << btf->beta.transpose(), btf->s2, 
      btf->l2, btf->o2.transpose(), btf->alpha;

    for (int j=0; j<P; ++j) {
      if ( std::isnan(history(i,j)) ) {
        Rcpp::Rcout << "whatch out dexp!, nan @ (" <<  i << "," << j << ")" << " with value " << history(i,j) << std::endl;
        return history;
      }
    }
  }
  return history;
}
