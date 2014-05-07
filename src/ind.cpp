#include <RcppEigen.h>

// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vec;
typedef Eigen::MappedSparseMatrix<double> MspMat;
typedef Eigen::SparseMatrix<double> spMat;

class ind {
public:
  // constructor
  ind(const Vec y_, const int k_, const  MspMat D_, const double alpha_, const double rho_, const int cprior_) : y(y_), k(k_), D(D_), cprior(cprior_), alpha(alpha_), rho(rho_) {

    // general info
    n = y.size();
    nk = n-k-1;

    // initialize
    if (cprior == 1) init_l(); else init_l2(); // which conditional prior
    init_o2();
    init_s2();
    init_beta();

    LLt.analyzePattern(sigma); // symbolic decomposition on the sparsity
    LLt.setShift(1.0, 1.0);    // add I to Sigma_f

  }

  // setters
  void set_s2(const double& val) {
    s2 = val;
  }
  void set_l2(const double& val) {
    l2 = val;
  }
  template <typename T>
  void set_l(const Eigen::DenseBase<T>& val) {
    l = val;
  }
  template <typename T>
  void set_o2(const Eigen::DenseBase<T>& val) {
    o2 = val;
  }
  template <typename T>
  void set_beta(const Eigen::DenseBase<T>& val) {
    beta = val;
  }
  template <typename T>
  void set_sigma(const Eigen::DenseBase<T>& val) {
    sigma = D.transpose()*mkDiag(val)*D;
    sigma.makeCompressed();
  }

  // random
  Vec rndNorm(const int& n_) {
    RNGScope scope;
    NumericVector x = rnorm(n_);
    Vec out(as<Vec>(x));
    return out;
  } 
  double rInvGauss(const double& nu_, const double& lambda_) {
    if (std::isnan(nu_) || std::isnan(lambda_)) {
      stop("rInvGauss: Invalid parameters.");
    }
    Vec z = rndNorm(1);          // one N(0,1)
    double z2 = z(0)*z(0);
    double nu2 = nu_*nu_;
    double x = nu_+0.5*z2*nu2/lambda_-(0.5*nu_/lambda_)*std::sqrt(4.0*nu_*lambda_*z2+nu2*z2*z2);
    Vec u = rndUniform(1);
    if ( u(0) < nu_/(nu_+x) ) {
      return x;
    } else {
      return nu2/x;
    }
  }
  template <typename T>
  Vec rndInvGauss(const Eigen::DenseBase<T>& nu_, const double& lambda_) {
    Vec out(nk);
    for (int i=0; i<nk; i++) {
      out(i) = rInvGauss(nu_(i), lambda_);
      if ( std::isnan(out(i)) ) stop("rndInvGauss: Invalid output.");
    }
    return out;
  }
  template <typename T, typename S>
  Vec rndInvGauss(const Eigen::DenseBase<T>& nu_, const Eigen::DenseBase<S>& lambda_) {
    Vec out(nk);
    for (int i=0; i<nk; i++) {
      out(i) = rInvGauss(nu_(i), lambda_(i));
    }
    return out;
  }
  // rndMVNorm templated ? can't figure it out
  Vec rndMVNorm(const Vec& mu_, const spMat& sqrtCov_, const double& scale) {
    return (scale*(sqrtCov_*rndNorm(sqrtCov_.cols()))+mu_).eval();
  } 
  Vec rndGamma(const int& n_, const double& shape_, const double& scale_) {
    if ( std::isnan(shape_) || std::isnan(scale_) ) {
      stop("rndGamma: Invalid parameters.");
    }
    RNGScope scope;
    NumericVector x = rgamma(n_, shape_, scale_);
    Vec out(as<Vec>(x));        // convert to Eigen::VectorXd
    return out;
  }
  Vec rndUniform(const int& n_) {
    RNGScope scope;
    NumericVector x = runif(n_);
    Vec out(as<Vec>(x));
    return out;
  }

  // utility
  template <typename T>
  spMat mkDiag(const Eigen::DenseBase<T>& val) {
    int I = val.size();
    spMat W(I,I); W.reserve(I);
    for (int i=0; i<I; i++) {
      W.insert(i,i) = val(i);    // insert along diagonal
    }
    return W;
  }

  int fact(const int& k_) {
    int out=1;
    for (int i=k_; i>1; i--) out *= i;
    return out;
  }

  // fields
  Vec y, beta, o2, l;
  int n, k, nk;
  spMat D, sigma;
  int cprior;
  double alpha, rho, l2, s2;
  Eigen::SimplicialLLT<spMat > LLt;

private:
  
  // initializers
  void init_l2() {
    Vec L = rndGamma(1, nk+alpha, 2.0/(1.0+rho));
    l2 = L(0);
  }
  void init_l() {
    Vec L = rndGamma(nk, nk+alpha, 1.0/(1.0+rho));
    l = L;
  }
  void init_o2() {
    Vec O = rndGamma(nk, 1.0, 2.0);
    set_sigma(O.cwiseInverse());
    o2 = O;
  }
  void init_s2() {
    Vec S = rndGamma(1, 0.1, 1.0);
    s2 = 1.0/S(0);
  }
  void init_beta() {
    beta = y;
  }
};
