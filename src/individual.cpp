#include <RcppEigen.h>
#include <vector>
#include <random>

typedef Eigen::VectorXd Vec;
typedef Eigen::Map<Vec> MVec;
typedef Eigen::MappedSparseMatrix<double> MspMat;
typedef Eigen::SparseMatrix<double> spMat;

class individual {
 private:
  /* fields */
  const int max_draws = 25;
  MVec y;
  int k;
  MspMat D;
  spMat sigma;
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
    double x = nu_+0.5*z2*(nu2/lambda_)-(0.5*(nu_/lambda_))*std::sqrt(4.0*nu_*lambda_*z2+nu2*z2*z2);
    Vec u = rndUniform(1);
    if ( u(0) < nu_/(nu_+x) ) {
      return x;
    } else {
      return nu2/x;
    }
  }
  template <typename T>
  Vec rndInvGauss(const Eigen::DenseBase<T>& nu_, const double& lambda_) {
    Vec out(nk); int draws;
    for (int i=0; i<nk; i++) {
      draws = 0;
      do {
        out(i) = rInvGauss(nu_(i), lambda_);        
        ++draws;
      } while ( (out(i) < 1e-13 || out(i) > 1e13) && draws < max_draws);
      if (draws > 1) Rcpp::Rcout << draws << std::endl;
    }
    return out;
  }
  template <typename T, typename S>
  Vec rndInvGauss(const Eigen::DenseBase<T>& nu_, const Eigen::DenseBase<S>& lambda_) {
    Vec out(nk); int draws;
    for (int i=0; i<nk; i++) {
      draws = 0;
      do {
        out(i) = rInvGauss(nu_(i), lambda_(i));
        ++draws;
      } while ( (out(i) < 1e-13 || out(i) > 1e13) && draws < max_draws);
      if (draws > 1) Rcpp::Rcout << draws << std::endl;
    }
    return out;
  }
  // rndMVNorm templated ? can't figure it out
  Vec rndMVNorm(const Vec& mu_, const spMat& sqrtCov_, const double& scale) {
    return (scale*(sqrtCov_*rndNorm(sqrtCov_.cols()))+mu_).eval();
  } 
  Vec rndGamma(const int& n_, const double& shape_, const double& scale_) {
    Rcpp::RNGScope scope; Rcpp::NumericVector x;
    int draws = 0; bool one; bool two;
    do {
      x = Rcpp::rgamma(n_, shape_, scale_);
      ++draws;
      one = Rcpp::is_true(Rcpp::any(x < 1e-13));
      two = Rcpp::is_true(Rcpp::any(x > 1e13));
    } while ( (one || two) && draws < max_draws);
    Vec out(Rcpp::as<Vec>(x));        // convert to Eigen::VectorXd
    return out;
  }
  template <typename T>
  Vec rndGamma(const double& shape_, const Eigen::DenseBase<T>& scale_) {
    Rcpp::NumericVector x(nk); int draws;
    Rcpp::RNGScope scope; Rcpp::NumericVector tmp;
    for (int i=0; i<nk; ++i) {
      draws = 0;
      do {
        tmp = Rcpp::rgamma(1, shape_, scale_(i));
        ++draws;
      } while ( (tmp(0) < 1e-13 || tmp(0) > 1e13) && draws < max_draws);
      x(i) = tmp(0);
    }
    Vec out(Rcpp::as<Vec>(x)); // convert to Eigen::VectorXd
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
  std::vector<double> pdfA(const Vec& a_) {
    int A = a_.size(); std::vector<double> out;
    double Db = 1.0 + (D*beta).lpNorm<1>()/(std::sqrt(s2)*rho);
    for (int i=0; i<A; ++i) {
      out.push_back(((1-a_(i))/a_(i))*std::pow(Db,-1.0*(nk - 1.0 + 1.0/a_(i))));
    }
    return out;
  }
  std::vector<double> pdfR(const Vec& r_) {
    int R = r_.size(); std::vector<double> out;
    double Db = (D*beta).lpNorm<1>()/std::sqrt(s2);
    for (int i=0; i<R; ++i) {
      out.push_back(r_(i)/(1-r_(i)) * std::pow(1.0+(r_(i)/(1-r_(i)))*Db, -1.0*(nk+alpha)));
    }
    return out;
  }

 public:
 
  /* fields */
  int n, nk;
  Vec beta, o2;
  double l, l2, s2, alpha, rho;

  /* constructor */
  individual(const MVec y_, const int k_, const  MspMat D_, const double alpha_, const double rho_) : y(y_), k(k_), D(D_), alpha(alpha_), rho(rho_) {

    // general info
    n = y.size();
    nk = n-k-1;

    // initialize
    init_o2();
    init_s2();
    init_beta();
    init_l2();
    init_l(); 

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

  // wrong full conditionals for double exponential conditional prior
  void upOmega2() {
    Vec eta = rndInvGauss(std::sqrt(l2*s2)/(D*beta).cwiseAbs().array(), l2);
    o2 = eta.cwiseInverse();

    // update sigma_f
    sigma = D.transpose()*mkDiag(eta)*D;
    sigma.makeCompressed();
  }
  void upLambda2() {
    Vec lambda2 = rndGamma(1, nk+alpha, 2.0/(o2.sum()+2*rho));
    l2 = lambda2(0);
  }

  // full conditionals for generalized double Pareto conditional prior 
  void upOmega() {
    Vec eta = rndInvGauss(l*std::sqrt(s2)/(D*beta).cwiseAbs().array(), l*l);
    o2 = eta.cwiseInverse();
    
    // update sigma_f
    sigma = D.transpose()*mkDiag(eta)*D;
    sigma.makeCompressed();
  }
  void upLambda() {
    double sig = std::sqrt(s2);
    Vec lambda = rndGamma(1, nk+alpha, sig/((D*beta).lpNorm<1>() + rho*sig));
    l = lambda(0);
  }

  // griddy Gibbs sampler
  void upA() {
    Vec uniA = rndUniform(250);
    std::vector<double> A = pdfA(uniA);
    std::random_device rdA; std::mt19937 genA(rdA());
    std::discrete_distribution<> dA(A.begin(), A.end());
    alpha = 1.0/uniA(dA(genA)) - 1.0;
  }
  void upR() {
    Vec uniR = rndUniform(100);
    std::vector<double> R = pdfR(uniR);
    std::random_device rdR; std::mt19937 genR(rdR());
    std::discrete_distribution<> dR(R.begin(), R.end());
    rho = 1.0/uniR(dR(genR)) - 1.0;
  }
};
