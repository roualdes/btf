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
    I = speye(n,n);

    // initialize
    init_D(n, k, nk);
    init_l2(rho, delta);
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
  void set_o2(colvec val) {
    o2 = val;
  }
  void set_beta(colvec val) {
    beta = val;
  }
  void set_sigma(colvec eta_) {
    mat W = diagmat(eta_);
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

  NumericVector get_gamma(int n_, double shape_, double rate_) {
    RNGScope scope;
    NumericVector x = rgamma(n_, shape_, rate_);
    return x;
  }

  // fields
  colvec y, o2, beta;
  int n, k, nk;
  double rho, delta, l2, s2;
  mat D, sigma, I;

  // TODO: use sp_mat

private:

  // initializers
  void init_l2(double rho_, double delta_) {
    NumericVector L = get_gamma(1, nk+rho_, 1+delta);
    l2 = L[0];
  }

  void init_o2(int nk_) {
    RNGScope scope;
    NumericVector O = rexp(nk_, l2/2.0);
    colvec O2(O.begin(), O.size(), false);
    set_sigma(1.0/O2);
    o2 = O2;
  }

  void init_s2() {
    NumericVector S = get_gamma(1, 1, 1);
    s2 = 1/S[0];
  }

  void init_beta() {
    beta = y;
  }

  void init_D(int n_, int k_, int nk_) {
    // TODO: use sp_mat
    //   int n1 = n_-1;
    //   int no_vals = n_+n_-1;      // number of values in d^(0)

    //   // create d^(0)
    //   // create values to be inserted into d^(0)
    //   colvec v = ones<colvec>(no_vals); // init values to be inserted into d^(0)
    //   v.rows(0, n1) *= -1.0;

    //   // create locations of values in d^(0)
    //   urowvec l = linspace<urowvec>(0, n1, n_);
    //   umat loc(2, no_vals);
    //   loc.submat( span(0,0), span(0, n1)) = l;
    //   loc.submat( span(1,1), span(0, n1)) = l;
    //   loc.submat( span(0,0), span(n_, no_vals-1)) = l.cols(0, n1-1);
    //   loc.submat( span(1,1), span(n_, no_vals-1)) = l.cols(1, n1);

    //   // fill D^(0)
    //   sp_mat d(loc, v);
    //   sp_mat d0 = d;
    //   if (k_ != 0) {
    //     int i=0;
    //     while (i < k_) {
    //       d *= d0;
    //       i++;
    //     }
    //   }
    //   D = d.rows(0, nk_-1);
    // }
    
    // using mat class until Armadillo has more sp_mat funcitonality
    mat d(n, n, fill::zeros);
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
};


// updaters
void upBeta(individual* i) {
  mat sigma_inv = inv_sympd(i->I + i->sigma);
  i->set_beta(i->rmvnorm(sigma_inv*i->y, i->s2*sigma_inv));
}

void upLambda(individual* i) {
  NumericVector l = i->get_gamma(1, i->nk+i->rho, i->delta+sum(i->o2)/2.0);
  i->set_l2(l[0]);
}

void upOmega(individual* i) {
  colvec Db = abs(i->D*i->beta);
  colvec m = sqrt(i->l2*i->s2)/Db;
  colvec eta(i->nk);
  for (int j=0; j<i->nk; j++) {
    eta(j) = i->rinvGauss(m(j), i->l2);
  }
  i->set_sigma(eta);
  i->set_o2(1.0/eta);
}

void upSig(individual* i) {
  colvec rate = dot(i->y-i->beta, i->y-i->beta) + i->beta.t()*i->sigma*i->beta;
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

