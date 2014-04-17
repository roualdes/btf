#include <RcppEigen.h>
#include <Ziggurat.h>
#include "ind.cpp"

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
typedef Eigen::SparseMatrix<double> spMat;

// updater
void upParams(ind* i, const double& lambda2) {

  // beta
  // old: two solves
  // i->LLt.factorize(i->sigma);
  // spMat I(i->n,i->n); I.setIdentity();
  // spMat SigmaInv(i->n, i->n); SigmaInv = i->LLt.solve(I);
  // Eigen::SimplicialLLT<spMat > chol; chol.compute(i->s2*SigmaInv);
  // i->set_beta(i->rndMVNorm(SigmaInv*i->y, chol.matrixL()));

  // new: one solve
  i->LLt.factorize(i->sigma);
  // to get L^{-1}, have to deal with fill-reducing permutation matrices
  // Eigen::MatrixXd Pinv = i->LLt.permutationPinv();
  // Eigen::MatrixXd Ltinv = i->LLt.solve(Pinv*i->LLt.matrixL()*Pinv.transpose());
  // some awkward casting to MatrixXds is going on; I would much prefer spMat

  spMat L(i->n,i->n);
  i->LLt.matrixL().twistedBy(i->LLt.permutationPinv()).evalTo(L);
  spMat Ltinv(i->n,i->n);
  Ltinv = i->LLt.solve(L);
  i->set_beta(i->rndMVNorm(Ltinv*(Ltinv.transpose()*i->y), Ltinv, std::sqrt(i->s2)));

  // lambda^2
  // did the user pre-specify lambda?
  if (lambda2 >= 0) {
    i->set_l2(lambda2);
  } else {
    // [lambda^2] ~ lambda^-2
    Vec l = i->rndGamma(1, i->nk+i->rho-1.0, 2.0/(i->o2.sum()+2*i->delta));
    // uniform prior
    // Vec l = i->rndGamma(1, i->nk+i->rho+1.0, 2.0/(i->o2.sum()+2*i->delta));
    i->set_l2(l(0));    
  }

  // omega^2
  Vec Db = (i->D*i->beta).cwiseAbs();
  Vec nu = std::sqrt(i->l2*i->s2)/Db.array();
  Vec eta = i->rndInvGauss(nu, i->l2);
  i->set_sigma(eta);
  i->set_o2(eta.cwiseInverse());

  // s^2
  double rate = (i->y-i->beta).squaredNorm() + i->beta.transpose()*i->sigma*i->beta;
  Vec s = i->rndGamma(1, i->n, 2.0/rate);
  i->set_s2(1.0/s(0));

}

RCPP_MODULE(ind) {
  class_<ind>( "ind" )
    .constructor<Vec, double, spMat, double, double>()
    .field( "y", &ind::y)
    .field( "beta", &ind::beta)
    .field( "o2", &ind::o2)
    .field( "l2", &ind::l2)
    .field( "s2", &ind::s2)
    .field( "D", &ind::D)
    .field( "k", &ind::k)
    .method( "upParams", &upParams)
    ;
}
