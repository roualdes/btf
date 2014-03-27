#include <RcppEigen.h>
#include <Ziggurat.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]] // this alone doesn't seem to do it

using namespace Rcpp;
typedef Eigen::VectorXd Vec;
typedef Eigen::SparseMatrix<double> spMat;
static Ziggurat::Ziggurat::Ziggurat zigg;

class ind {
public:
  // constructors
  ind(const Vec y_, const double k_, const spMat D_, const double rho_, const double delta_) : y(y_), k(k_), D(D_), rho(rho_), delta(delta_) {

    // general info
    n = y.size();
    nk = n-k-1;

    // initialize
    init_l2();
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
  void set_o2(const Vec& val) {
    o2 = val;
  }
  void set_beta(const Vec& val) {
    beta = val;
  }
  void set_sigma(const Vec& val) {
    spMat W = mkDiag(val);
    sigma = D.transpose()*W*D;
    sigma.makeCompressed();
  }

  // random
  Vec zrnorm(int n) {
    Vec x(n);
    for(int i=0; i<n; i++) {
      x(i) = zigg.norm();
    }
    return x;
  } 
  double rInvGauss(const double& nu_, const double& lambda_) {
    Vec z = zrnorm(1);          // one N(0,1)
    double z2 = z(0)*z(0);
    double nu2 = nu_*nu_;
    double x = nu_+0.5*z2*nu2/lambda_-(0.5*nu_/lambda_)*std::sqrt(4.0*nu_*lambda_*z2+nu2*z2*z2);
    Vec u = Vec::Random(1);
    if ( std::fabs(u(0)) < nu_/(nu_+x) ) {
      return x;
    } else {
      return nu2/x;
    }
  }
  Vec rndInvGauss(const Vec& nu_, const double& lambda_) {
    Vec out(nk);
    for (int i=0; i<nk; i++) {
      out(i) = rInvGauss(nu_(i), lambda_);
    }
    return out;
  }
  // Vec rndMVNorm(const Vec& mu_, const Eigen::MatrixXd& sqrtCov_) {
  //   Vec Y = zrnorm(sqrtCov_.cols());
  //   return mu_+sqrtCov_*Y;
  // }
  Vec rndMVNorm(const Vec& mu_, const spMat& sqrtCov_, const double& scale) {
    return (scale*(sqrtCov_*zrnorm(sqrtCov_.cols()))+mu_).eval();
  }
  Vec rndGamma(const int& n_, const double& shape_, const double& rate_) {
    RNGScope scope;
    NumericVector x = rgamma(n_, shape_, rate_);
    Vec out(as<Vec>(x));        // convert to Eigen::VectorXd
    return out;
  }

  // utility
  spMat mkDiag(const Vec& val) {
    int I = val.size();
    spMat W(I,I); W.reserve(I);
    for (int i=0; i<I; i++) {
      W.insert(i,i) = val(i);    // insert along diagonal
    }
    return W;
  }

  // fields
  Vec y, beta, o2;
  int n, k, nk;
  spMat D, sigma;
  double rho, delta, l2, s2;
  Eigen::SimplicialLLT<spMat > LLt;

private:
  
  // initializers
  void init_l2() {
    Vec L = rndGamma(1, nk+rho, 1.0+delta);
    l2 = L(0);
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
