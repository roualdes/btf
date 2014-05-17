#include "individual.cpp"

// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::export]]
Eigen::MatrixXd gdPBTF(const int& iter,
                        const Eigen::Map<Eigen::VectorXd>& y, 
                        const int& k, const Eigen::MappedSparseMatrix<double>& D, 
                        const double& alpha, const double& rho) {

  int s = 0;
  int tmp = 0;
  // initialize btf object
  individual *btf; double a;
  if ( alpha < 0.0 ) a = 1.0; else a = alpha;
  btf = new individual(y, k, D, a, rho); 

  // initialize matrix of posterior draws
  int P = btf->n + btf->nk + 3; // count number of parameters
  Eigen::MatrixXd history = Eigen::MatrixXd::Zero(iter, P); 

  // runs sampler
  for (int i=0; i<iter; ++i) {
    btf->upBeta(); btf->upS2();
    btf->upLambda(); btf->upOmega(); 
    if ( alpha < 0.0 ) btf->upA(); // treat alpha as parameter
    // btf->upR();
    history.row(i) << btf->beta.transpose(), btf->s2, 
      btf->l, btf->o2.transpose(), btf->alpha;

    for (int j=0; j<P; ++j) {
      if ( std::isnan(history(i,j)) ) {
        s = 1;
      }
    }
    if (s == 1) {
      Rcpp::Rcout << "stopped at iteration " << i << std::endl;
      tmp = i-1;
      Rcpp::Rcout << history.row(tmp) << std::endl;
      Rcpp::Rcout << history.row(i) << std::endl;
      Rcpp::stop("some nans.");      
    }
  }
  return history;
}
