#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

class individual {              // author Edward A. Roualdes
public:
  // constructor
  individual(colvec y_, double k_, double rho_, double delta_ ) : y(y_), k(k_), rho(rho_), delta(delta_) {

    // general info
    n = y.n_rows;
    nk = n-k-1;
    I = eye(n,n);               // I = speye(n,n); 
    // x = linspace<colvec>(1, 100, n) / n; // TODO as argument of constructor

    // initialize
    init_D(n, nk, k); // init_spD(n, nk, k);    // init_Delta(n, nk, k, x); 
    init_l2(rho_, delta_);
    init_o2(nk);
    init_s2();
    init_beta();
  }

  // setters
  void set_s2(double val) {
    s2 = val;
  }
  void set_l2(double val) {
    l2 = val;
  }
  void set_o2(const colvec& val) {
    o2 = val;
  }
  void set_beta(const colvec& val) {
    beta = val;
  }
  void set_sigma(const colvec& val) {
    mat W = diagmat(val);
    sigma = D.t()*W*D;
  }

  // random
  double rinvGauss(double nu_, double lambda_) {
    double y = randn();         // one N(0,1)
    double y2 = y*y;
    double nu2 = nu_*nu_;
    double x = nu_ + 0.5*y2*nu2/lambda_ -(0.5*nu_/lambda_)*sqrt(4.0*nu_*lambda_*y2+nu2*y2*y2);
    double z = randu();         // one U(0,1)
    if( z < nu_/(nu_+x) ) {
      return x;
    } else {
      return nu2/x;
    } 
  }
  colvec rmvnorm(const colvec& mu_, const mat& sigma_) {
    colvec Y = randn<colvec>(sigma_.n_cols);
    mat b = mu_.t() + Y.t()*chol(sigma_);
    return conv_to<colvec>::from(b);
  }
  colvec rmvnorm2(const colvec& mu_, const mat& retval_) {
    colvec Y = randn<colvec>(retval_.n_cols);
    mat b = mu_.t() + Y.t()*retval_;
    return conv_to<colvec>::from(b);
  }
  NumericVector get_gamma(int n_, double shape_, double rate_) {
    RNGScope scope;
    NumericVector x = rgamma(n_, shape_, rate_);
    return x;
  }

  // fields
  colvec y, o2, beta, dx, x;
  int n, k, nk;
  double rho, delta, l2, s2;
  mat D, sigma, I;
  // sp_mat D, sigma, I;

private:

  // initializers
  void init_l2(double rho_, double delta_) {
    NumericVector L = get_gamma(1, nk+rho_, 1.0+delta_);
    l2 = L[0];
  }
  void init_o2(int nk_) {
    RNGScope scope;
    NumericVector O = rexp(nk_, 0.5);
    colvec O2(O.begin(), O.size(), false);
    set_sigma(1.0/O2);
    o2 = O2;
  }
  void init_s2() {
    NumericVector S = get_gamma(1, 0.1, 1.0); // user specified like lambda?
    s2 = 1/S[0];                              // sigma2 ~ invgamma
  }
  void init_beta() {
    beta = y;                   // reasonable?
  }
  void init_D(int n_, int nk_, int k_) {
    mat d(n_, n_, fill::zeros);
    d.diag() -= 1.0;
    d.diag(1) += 1.0;
    mat d0 = d;
      if (k_ != 0) {
        int i=0;
        while (i<k_) {
          d *= d0;
          i++;
        }
      }
      D = d.rows(0, (nk_-1));
    }
  void init_Delta(int n_, int nk_, int k_, colvec x_) {
    mat d(nk_, n_, fill::zeros);
    double fk = factorial(k_);
    for (int i=0; i<nk_; i++) {
      int k1 = i+k_+1;
      double diff = x_(i) - x_(k1);
      colvec subx = x_.rows(i,k1);
      for (int j=i; j<k1+1; j++) {
        d(i,j) = (diff*fk)/omega(subx, x_(j));
      }
    }
    D = d;
  }
  void init_spD(int n_, int nk_, int k_) {
    // initialize first difference matrix
    int N = n_+(n_-1);          // number of elements to fill d
    colvec v(N, fill::ones);
    v.rows(n_, N-1) *= -1.0;    // values of sparse matrix
    umat loc(2, N);             // locations
    uvec pos(2);                // helps fill loc
    for (int i=0; i<n_; i++) {
      pos.fill(i);     
      loc.col(i) = pos;         // fill diag
      if (i < n_-1) {
        pos(1) += 1.0;
        loc.col(i+n_) = pos;      // fill one above diag        
      }
    }
    sp_mat d(loc, v);
    // compute difference matrix recursively
    sp_mat d0 = d;
    if (k_ != 0) {
      int i=0;
      while (i < k_) {
        d *= d0;
        i++;
      }
    }
    D = d.rows(0, nk_-1);
  }

  // utility
  double factorial(double k_) {
    double out = 1.0;
    if (k_ > 1) {
      out = k_*factorial(k_-1);
    } 
    return out;
  }
  double omega(colvec x_, double xj_) {
    colvec diffs = xj_ - x_;
    double out = 1;
    for (int j=0; j<x_.n_rows; j++) {
      if (diffs[j] != 0.0) {
        out *= diffs[j];
      }
    }
      return(out);
  }
};



// updaters
void upBeta(individual* i) {
  // invert matrix
  mat siginv = inv_sympd(i->I+i->sigma);
  i->set_beta(i->rmvnorm(siginv*i->y, i->s2*siginv));

  // svd
  // doesn't seem to work, but I can't figure out why not
  // mat U;
  // vec s;
  // mat V;
  // mat E = i->I+i->sigma;
  // svd(U,s,V, E);
  // i->set_beta(i->rmvnorm2(U*diagmat(1/s)*V.t()*i->y,
  //                         U*diagmat(i->s2/sqrt(s))*V.t()));

  // sparse version if I can ever use it
  // vec eigval;
  // mat eigvec;
  // sp_mat E = i->I+i->sigma;
  // eigs_sym(eigval, eigvec, E, i->n);
  
  // i->set_beta(i->rmvnorm2((eigvec*diagmat(1/eigval)*eigvec.t())*i->y,
  //                         eigvec*diagmat(i->s2/sqrt(eigval))*eigvec.t()));
}
void upLambda(individual* i) {
  // gamma (rho, delta) prior
  // NumericVector l = i->get_gamma(1, i->nk+i->rho, i->delta+sum(i->o2)/2.0);
  // lambda^{-2} prior: requires n-k > 2
  NumericVector l = i->get_gamma(1, i->nk-1, 2.0/sum(i->o2));
  i->set_l2(l[0]);
}
void upOmega(individual* i) {
  colvec Db = abs(i->D*i->beta);
  double l = i->l2;
  colvec m = sqrt(l*i->s2)/Db;
  colvec eta(i->nk);
  for (int j=0; j<i->nk; j++) {
    eta(j) = i->rinvGauss(m(j), l);
  }
  i->set_sigma(eta);
  i->set_o2(1.0/eta);
}
void upSig(individual* i) {
  colvec res = i->y - i->beta;
  colvec rate = dot(res, res) + i->beta.t()*i->sigma*i->beta;
  NumericVector sig2 = i->get_gamma(1, i->n, 2.0/rate[0]);
  i->set_s2(1/sig2[0]);
}

// to R
RCPP_MODULE(individual) {
  class_<individual>( "individual" )
    .constructor<colvec, double, double, double>()
    .field( "n", &individual::n)
    .field( "y", &individual::y)
    .field( "beta", &individual::beta)
    .field( "k", &individual::k)
    .field( "delta", &individual::delta)
    .field( "rho", &individual::rho)
    .field( "s2", &individual::s2)
    .field( "l2", &individual::l2)
    .field( "o2", &individual::o2)
    .field( "D", &individual::D)
    .field( "sigma", &individual::sigma)
    .field( "nk", &individual::nk)
    .method( "upOmega", &upOmega)
    .method( "upBeta", &upBeta)
    .method( "upSig", &upSig)
    .method( "upLambda", &upLambda)
    .method( "rgam", &individual::get_gamma)
    .method( "rmvnorm", &individual::rmvnorm)
    .method( "rinvGauss", &individual::rinvGauss)
    ;
}
