#include <RcppEigen.h>

typedef Eigen::VectorXd Vec;
typedef Eigen::MappedSparseMatrix<double> MspMat;
typedef Eigen::SparseMatrix<double> spMat;

class individual {
 private:
  /* fields */
  Vec y;
  int k;
  spMat D, sigma;
  double alpha, rho;
  Eigen::SimplicialLLT<spMat > LLt;

  /* random */
  Vec rndNorm(const int& n_) {
    Rcpp::RNGScope scope;
    Rcpp::NumericVector x = Rcpp::rnorm(n_);
    Vec out(Rcpp::as<Vec>(x));
    return out;
  } 
  double rInvGauss(const double& nu_, const double& lambda_) {
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
    Rcpp::RNGScope scope;
    Rcpp::NumericVector x = Rcpp::rgamma(n_, shape_, scale_);
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
    Vec L = rndGamma(nk, nk+alpha, 1.0/(1.0+rho));
    l = L;
  }
  void init_o2() {
    Vec O = rndGamma(nk, 1.0, 2.0);
    sigma = D.transpose()*mkDiag(O.cwiseInverse())*D;
    o2 = O;
  }
  void init_s2() {
    Vec S = rndGamma(1, 0.1, 1.0);
    s2 = 1.0/S(0);
  }
  void init_beta() {
    beta = y;
  }

  /* utility */
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
  Vec beta, o2, l;
  double l2, s2;
  int n, nk;

  /* constructor */
  individual(const Vec y_, const int k_, const  MspMat D_, const double alpha_, const double rho_) : y(y_), k(k_), D(D_), alpha(alpha_), rho(rho_) {

    // general info
    n = y.size();
    nk = n-k-1;

    // initialize
    init_l2();
    init_l(); 
    init_o2();
    init_s2();
    init_beta();

    LLt.analyzePattern(sigma); // symbolic decomposition on the sparsity
    LLt.setShift(1.0, 1.0);    // add I to Sigma_f

  }

  /* update parameters */
  void upBeta() {
    LLt.factorize(sigma);
    spMat L(n,n); spMat Ltinv(n,n);
    LLt.matrixL().twistedBy(LLt.permutationPinv()).evalTo(L);
    Ltinv = LLt.solve(L);
    beta = rndMVNorm(Ltinv*(Ltinv.transpose()*y), Ltinv, std::sqrt(s2));
  }
  void upS2() {
    double rate = (y-beta).squaredNorm() + beta.transpose()*sigma*beta;
    Vec S = rndGamma(1, n, 2.0/rate);
    s2 = 1.0/S(0);
  }
  void upOmega2() {
    Vec eta = rndInvGauss(std::sqrt(l2*s2)/(D*beta).cwiseAbs().array(), l2);
    o2 = eta.cwiseInverse();

    // update sigma_f
    sigma = D.transpose()*mkDiag(eta)*D;
    sigma.makeCompressed();
  }

  void upLambda2() {
    double pnk = std::pow((double)n, 2*k);
    Vec lambda2 = rndGamma(1,nk+alpha, 2.0/(o2.sum()+2*rho/pnk));
    l2 = lambda2(0);
    
  }
  void upOmega() {
    Vec eta = rndInvGauss(l.array()*std::sqrt(s2)/(D*beta).cwiseAbs().array(),l.cwiseProduct(l));
    o2 = eta.cwiseInverse();
    
    // update sigma_f
    sigma = D.transpose()*mkDiag(eta)*D;
    sigma.makeCompressed();    
  }
  void upLambda() {
    Vec lambda(nk); 
    // double pnk = std::pow((double)n, 2*k);
    Vec Db = (D*beta).cwiseAbs();
    double a = 1+alpha; 
    double sig = std::sqrt(s2);
    for (int j=0; j<nk; j++) {
      // maybe change this to one function call
      Vec tmp = rndGamma(1, a, 1.0/((Db(j)/sig + rho)));
      lambda(j) = tmp(0);
    }
    l = lambda;    
  }
};
