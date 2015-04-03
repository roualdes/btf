#include "individual.cpp"

// [[Rcpp::export]]
Rcpp::List gdp(const int& iter,
                       const Eigen::Map<Eigen::VectorXd>& y, 
                       const Eigen::MappedSparseMatrix<double>& D, 
                       const double& alpha, const double& rho,
                       const int& m, const bool& debug) {

  // initialize btf object
  individual *btf;
  btf = new individual(y, D, alpha, rho);
  bool broken = false;

  // initialize matrices of posterior draws
  Eigen::MatrixXd beta_draws = Eigen::MatrixXd::Zero(iter, btf->n);
  Eigen::MatrixXd omega_draws = Eigen::MatrixXd::Zero(iter, btf->nk);
  Eigen::VectorXd s2_draws = Eigen::VectorXd::Zero(iter);
  Eigen::VectorXd lambda_draws = Eigen::VectorXd::Zero(iter);
  
  // run sampler
  for (int i=0; i<iter; ++i) {
    if (i % m == 0) {
      btf->upBeta();
      beta_draws.row(i) = btf->beta.transpose();
    }
    btf->upS2();
    s2_draws(i) = btf->s2;
    btf->upLambda();
    lambda_draws(i) = btf->l;
    btf->upOmega(); 
    omega_draws.row(i) = btf->o2.transpose();

    if ( debug ) {
      for (int j=0; j<btf->n; ++j) {
        if (std::isnan(beta_draws(i,j))) {
          Rcpp::Rcout << "Warning: watch out gdp!, nan @ beta("
                      <<  i << "," << j << ")" << std::endl;
          broken = true;
        }
      }
      for (int j=0; j<btf->nk; ++j) {
        if (std::isnan(omega_draws(i,j))) {
          Rcpp::Rcout << "Warning: watch out gdp!, nan @ omega("
                      <<  i << "," << j << ")" << std::endl;
          broken = true;
        }
      }
      if (std::isnan(s2_draws(i)) || std::isnan(lambda_draws(i))) {
          Rcpp::Rcout << "Warning: watch out gdp!, nan @ s2|lambda("
                      <<  i << ")" << std::endl;
          broken = true;
      }
    }
    if (broken) break;
  }
  return Rcpp::List::create(Rcpp::Named("beta") = Rcpp::wrap(beta_draws),
                            Rcpp::Named("s2") = s2_draws,
                            Rcpp::Named("lambda") = lambda_draws,
                            Rcpp::Named("omega") = omega_draws);
}
