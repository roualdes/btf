#include "individual.cpp"

// [[Rcpp::export]]
Eigen::MatrixXd dexpBTF(const int& iter,
                        const Eigen::Map<Eigen::VectorXd>& y, 
                        const Eigen::MappedSparseMatrix<double>& D, 
                        const double& alpha, const double& rho, const bool& debug) {
  
  // initialize btf object
  individual *btf;
  btf = new individual(y, D, alpha, rho); 

  // initialize matrix of posterior draws
  int P = btf->n + btf->nk + 2; // count number of parameters
  Eigen::MatrixXd history = Eigen::MatrixXd::Zero(iter, P); 

  // runs sampler
  for (int i=0; i<iter; ++i) {
    btf->upBeta();
    btf->upS2();
    btf->upLambda2();
    btf->upOmega2();

    history.row(i) << btf->beta.transpose(),
      btf->s2, 
      btf->l2,
      btf->o2.transpose();

    if ( debug ) {
      for (int j=0; j<P; ++j) {
        if ( std::isnan(history(i,j)) ) {
          Rcpp::Rcout << "Warning: watch out dexp!, nan @ ("
                      <<  i << "," << j << ")" << std::endl;
          return history;
        }
      }      
    }
  }
  return history;
}
