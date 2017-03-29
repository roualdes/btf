#include <RcppEigen.h>
#include <vector>
#include <random>

// [[Rcpp::depends(RcppEigen)]]

typedef Eigen::VectorXd Vec;
typedef Eigen::Map<Vec> MVec;
typedef Eigen::MappedSparseMatrix<double> MspMat;
typedef Eigen::SparseMatrix<double> spMat;

class individual {
private:
  /* fields */
  MVec y;
  int max_draws;
  MspMat D;
  spMat sigma;
  spMat I;
  Eigen::SimplicialLLT<spMat > LLt;

  /* random */
  Vec rndNorm(const int& n_) {
    Rcpp::RNGScope scope;
    Rcpp::NumericVector x = Rcpp::rnorm(n_);
    Vec out(Rcpp::as<Vec>(x));
    return out;
  } 
  double rInvGauss(const double& nu_, const double& lambda_) {
    Vec z = rndNorm(1);         // one N(0,1)
    double z2 = z(0)*z(0);
    double nu2 = nu_*nu_;
    double c = 0.5*nu_/lambda_;
    double x = nu_ + c*z2*nu_ - c*std::sqrt(4.0*nu_*lambda_*z2 + nu2*z2*z2);
    Vec u = rndUniform(1);
    double out = (u(0) < nu_/(nu_+x)) ? x : nu2/x;
    return out;
  }
  template <typename T>
  Vec rndInvGauss(const Eigen::DenseBase<T>& nu_, const double& lambda_) {
    Vec out(nk); int draws;
    for (int i=0; i<nk; i++) {
      draws = 0;
      do {
        out(i) = rInvGauss(nu_(i), lambda_);
        ++draws;
      } while ( out(i) < 1e-11 && draws < max_draws );
      if ( draws > 1 ) Rcpp::Rcout << "Note: InvGauss resampled." << std::endl;
    }
    return out;
  }
  // rndMVNorm templated?
  Vec rndMVNorm(const Vec& mu_, const spMat& sqrtCov_, const double& scale) {
    return (scale*(sqrtCov_*rndNorm(sqrtCov_.cols()))+mu_).eval();
  } 
  Vec rndGamma(const int& n_, const double& shape_, const double& scale_) {
    Rcpp::RNGScope scope; Rcpp::NumericVector x;
    int draws = 0; bool toosmall;
    do {                        // cut off tail
      x = Rcpp::rgamma(n_, shape_, scale_);
      toosmall = Rcpp::is_true( Rcpp::any( x < 1e-11) );
      ++draws;
    } while ( toosmall && draws < max_draws );
    if ( draws > 1 ) Rcpp::Rcout << "Note: Gamma resampled." << std::endl;
    Vec out(Rcpp::as<Vec>(x));        // convert to Eigen::VectorXd
    return out;
  }
  Vec rndUniform(const int& n_) {
    Rcpp::RNGScope scope;
    Rcpp::NumericVector x = Rcpp::runif(n_);
    Vec out(Rcpp::as<Vec>(x));
    return out;
  }

  /* initializers */
  void init_l2() {
    Vec L = rndGamma(1, nk+alpha, 2.0/(1.0+rho));
    l2 = L(0);
  }
  void init_l() {
    Vec L = rndGamma(1, nk+alpha, 2.0/(1.0+rho));
    l = L(0);
  }
  void init_o2() {
    Vec O = rndGamma(nk, 1.0, 2.0);
    sigma = D.transpose()*mkDiag(O.cwiseInverse())*D;
    o2 = O;
  }
  void init_s2() {
    Vec S = rndGamma(1, 1.0, 0.5);
    s2 = 1.0/S(0);
  }
  void init_beta() {
    beta = y;
  }

  /* utility */
  // spMat mkDiag(const int& sz) {
  //   spMat W(sz, sz); W.reserve(sz);
  //   for (int i=0; i<sz; i++) {
  //     W.insert(i,i) = 1.0;
  //   }
  //   return W;
  // }
  template <typename T>
  spMat mkDiag(const Eigen::DenseBase<T>& val) {
    int I = val.size();
    spMat W(I,I); W.reserve(I);
    for (int i=0; i<I; i++) {
      W.insert(i,i) = val(i);    // insert along diagonal
    }
    return W;
  }

 public:
 
  /* fields */
  int n, nk;
  Vec beta, o2;
  double l, l2, s2, alpha, rho, Dbl1;
  Vec Db;

  /* constructor */
  individual(const MVec y_, const  MspMat D_, const double alpha_, const double rho_) : y(y_), D(D_), alpha(alpha_), rho(rho_) {

    // general info
    max_draws = 10;
    n = y.size();
    nk = D.rows();

    // initialize
    init_o2();
    init_s2();
    init_beta();
    init_l2();
    init_l(); 

    Db = D*beta;               // initialize D*beta;
    LLt.analyzePattern(sigma); // symbolic decomposition on the sparsity
    LLt.setShift(1.0, 1.0);    // add I to Sigma_f
  }

  /* update parameters */
  void upBeta() {
    LLt.factorize(sigma);
    spMat L(n,n); spMat Ltinv(n,n);
    LLt.matrixL().twistedBy(LLt.permutationPinv()).evalTo(L);
    Ltinv = LLt.solve(L);
    int draws = 0;
    do {
      beta = rndMVNorm(Ltinv*(Ltinv.transpose()*y), Ltinv, std::sqrt(s2));
      Db = D*beta;
      Dbl1 = (Db).lpNorm<1>();
      ++draws;
    } while ( (Db.cwiseAbs().array() < 1e-10).any() && draws < max_draws );
    if (draws > 1) Rcpp::Rcout << "Note: Normal resampled." << std::endl;
  }
  void upS2() {
    double rate = (y-beta).squaredNorm() + beta.transpose()*sigma*beta;
    Vec S = rndGamma(1, n, 2.0/rate);
    s2 = 1.0/S(0);
  }

  // dexp
  void upOmega2() {
    Vec eta = rndInvGauss(Db.cwiseAbs().cwiseInverse()*std::sqrt(l2*s2), l2);
    if ( (o2.array() <= 0.).any() ) {
      Rcpp::Rcout << "Warning: At least one omega <= zero..." << std::endl;
      eta = eta.cwiseAbs();
    }
    o2 = eta.cwiseInverse();    

    // update sigma_f
    sigma = D.transpose()*mkDiag(eta)*D;
    sigma.makeCompressed();
  }
  void upLambda2() {
    double tmp = o2.sum();
    // Rcpp::Rcout << "sum of omegas = " << tmp << std::endl;
    Vec lambda2 = rndGamma(1, nk+alpha, 2.0/(tmp+2*rho));
    l2 = lambda2(0);
  }

  // gdP
  void upOmega() {
    Vec eta = rndInvGauss(Db.cwiseAbs().cwiseInverse()*(l*std::sqrt(s2)), l*l);
    o2 = eta.cwiseInverse();
    if ( (o2.array() < 0.).any() ) {
      Rcpp::Rcout << "Warning: At least one omega <= zero..." << std::endl;
      eta = eta.cwiseAbs();
    }
    
    // update sigma_f
    sigma = D.transpose()*mkDiag(eta)*D;
    sigma.makeCompressed();
  }
  void upLambda() {
    double sig = std::sqrt(s2);
    Vec lambda = rndGamma(1, nk+alpha, sig/((D*beta).lpNorm<1>() + rho*sig));
    l = lambda(0);
  }
};

