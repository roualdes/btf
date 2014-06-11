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
  MVec y;
  Vec Db;
  int k, max_draws;
  MspMat D;
  spMat sigma;
  double c, c2;
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
    // Rcpp::Rcout << "invGauss inputs nu = " << nu_ << " and lambda = " << lambda_ << std::endl;
    for (int i=0; i<nk; i++) {
      draws = 0;
      do {
        out(i) = rInvGauss(nu_(i), lambda_);
        ++draws;
      } while ( (/*out(i) > 1e14 ||*/ out(i) < 1e-11) && draws < max_draws );
      if ( draws > 1 ) Rcpp::Rcout << draws << " for InvGauss." << std::endl;
    }
    return out;
  }
  // rndMVNorm templated ? can't figure it out
  Vec rndMVNorm(const Vec& mu_, const spMat& sqrtCov_, const double& scale) {
    return (scale*(sqrtCov_*rndNorm(sqrtCov_.cols()))+mu_).eval();
  } 
  Vec rndGamma(const int& n_, const double& shape_, const double& scale_) {
    Rcpp::RNGScope scope; Rcpp::NumericVector x;
    // Rcpp::Rcout << "Gamma inputs shape_ = " << shape_ << " and scale = " << scale_ << std::endl;
    int draws = 0; bool toobig; bool toosmall;
    do {                        // cut off tail
      x = Rcpp::rgamma(n_, shape_, scale_);
      toobig = Rcpp::is_true( Rcpp::any( x > 1e11) );
      toosmall = Rcpp::is_true( Rcpp::any( x < 1e-11) );
      if ( toosmall ) Rcpp::Rcout << "too small Gamma draw with shape = " << shape_ << " and scale = " << scale_ << std::endl;
      // if ( toobig ) Rcpp::Rcout << "too big Gamma draw with shape = " << shape_ << " and scale = " << scale_ << std::endl;
      ++draws;
    } while ( (toosmall /*|| toobig*/) && draws < max_draws );
    if ( draws > 1 ) Rcpp::Rcout << draws << " for Gamma." << std::endl;
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
  int fact(const int& k) {
    if (k <= 1) return k;
    return k*fact(k-1);
  }

 public:
 
  /* fields */
  int n, nk;
  Vec beta, o2;
  double l, l2, s2, alpha, rho;

  /* constructor */
  individual(const MVec y_, const int k_, const  MspMat D_, const double alpha_, const double rho_) : y(y_), k(k_), D(D_), alpha(alpha_), rho(rho_) {

    // general info
    max_draws = 10;
    n = y.size();
    nk = n-k-1;
    c = double(fact(k)/std::pow(n, k));
    c2 = double(std::pow(c, 2));

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
      ++draws;
    } while ( (Db.cwiseAbs().array() < 1e-10).any() && draws < max_draws );
    if (draws > 1) Rcpp::Rcout << draws << " for Beta." << std::endl;
  }
  void upS2() {
    double rate = (y-beta).squaredNorm() + beta.transpose()*sigma*beta;
    Vec S = rndGamma(1, n, 2.0/rate);
    s2 = 1.0/S(0);
  }

  // dexp
  void upOmega2() {
    double tmp = (std::sqrt(l2*s2));
    // Rcpp::Rcout << "invGauss scale = " << l2 << " and mean = " << tmp << std::endl;
    Vec eta = rndInvGauss(Db.cwiseAbs().cwiseInverse()*tmp, l2);
    if ( (o2.array() <= 0.).any() ) {
      Rcpp::Rcout << "why is at least one less than omega zero?" << std::endl;
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
    if ( (o2.array() < 0.).any() ) Rcpp::Rcout << "why is at least one less than omega zero?" << std::endl;
    
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
