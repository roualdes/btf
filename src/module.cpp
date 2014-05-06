#include <RcppEigen.h>
#include <Ziggurat.h>
#include "ind.cpp"

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins("cpp11")]]

void upParams(ind* i, const double& lambda2) {

  // beta
  i->LLt.factorize(i->sigma);
  spMat L(i->n,i->n); spMat Ltinv(i->n,i->n);
  i->LLt.matrixL().twistedBy(i->LLt.permutationPinv()).evalTo(L);
  Ltinv = i->LLt.solve(L);
  i->set_beta(i->rndMVNorm(Ltinv*(Ltinv.transpose()*i->y), Ltinv, std::sqrt(i->s2)));

  // s2
  double rate = (i->y-i->beta).squaredNorm() + i->beta.transpose()*i->sigma*i->beta;
  Vec S = i->rndGamma(1, i->n, 2.0/rate);
  if ( std::isnan(S(0))) stop("s2 not real.");
  i->set_s2(1.0/S(0));

  Vec eta(i->nk);
  Vec Db = (i->D*i->beta).cwiseAbs();
  double pnk = std::pow((double)i->n, i->k);
  double kf = (double)i->fact(i->k);
  if (i->cprior == 1) {         // generalized double Pareto

    Vec lambda(i->nk); 
    double a = 1+i->alpha; 
    double sig = std::sqrt(i->s2);
    for (int j=0; j<i->nk; j++) {
      // maybe change this to one function call
      Vec tmp = i->rndGamma(1, a, 1.0/((Db(j)/sig + i->rho)));
      lambda(j) = tmp(0);
    }
    i->set_l(lambda);
    eta = i->rndInvGauss(i->l.array()*sig/Db.array(),i->l.cwiseProduct(i->l));

  } else {

    if (lambda2 >= 0) {         // did the user pre-specify lambda?
      i->set_l2(lambda2);
    } else {                    // Laplacian

      Vec lambda2 = i->rndGamma(1,i->nk+i->alpha, pnk*2.0/((i->o2.sum()+2*i->rho)*kf));
      if ( std::isnan(lambda2(0)) ) stop("lambda2 not real.");
      i->set_l2(lambda2(0));

    }
    eta = i->rndInvGauss(std::sqrt(i->l2*i->s2)/Db.array(), i->l2);  
  }
  i->set_sigma(eta);
  i->set_o2(eta.cwiseInverse());
}

RCPP_MODULE(ind) {
  class_<ind>( "ind" )
    .constructor<Vec, double, MspMat, int, double, double>()
    .field( "y", &ind::y)
    .field( "beta", &ind::beta)
    .field( "o", &ind::o2)
    .field( "l", &ind::l)
    .field( "l2", &ind::l2)
    .field( "s2", &ind::s2)
    .field( "D", &ind::D)
    .field( "k", &ind::k)
    .method( "upParams", &upParams)
    ;
}
