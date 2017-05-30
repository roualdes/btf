#ifndef INDIVIDUAL_H_
#define INDIVIDUAL_H_

#include <RcppEigen.h>
#include <vector>
#include <random>

typedef Eigen::VectorXd Vec;
typedef Eigen::Map<Vec> MVec;
typedef Eigen::MappedSparseMatrix<double> MspMat;
typedef Eigen::SparseMatrix<double> spMat;

class Individual {
 private:
  MVec y;
  int max_draws;
  MspMat D;
  spMat sigma;
  spMat I;
  Eigen::SimplicialLLT<spMat > LLt;

 public:
  int n, nk;
  Vec beta, o2;
  double l, l2, s2, alpha, rho, Dbl1;
  Vec Db;
  void upBeta();
  void upS2();
  void upOmega2();
  void upLambda2();
  void upOmega();
  void upLambda();
  Vec rndNorm(const int& n_);
  double rInvGauss(const double& nu_, const double& lambda_);
  template <typename T>
    Vec rndInvGauss(const Eigen::DenseBase<T>& nu_, const double& lambda_);
  Vec rndMVNorm(const Vec& mu_, const spMat& sqrtCov_, const double& scale);
  Vec rndGamma(const int& n_, const double& shape_, const double& scale_);
  Vec rndUniform(const int& n_);
  void init_l2();
  void init_l();
  void init_o2();
  void init_s2();
  void init_beta();
  template <typename T>
    spMat mkDiag(const Eigen::DenseBase<T>& val);


  Individual(const MVec y_, const  MspMat D_, const double alpha_, const double rho_);
  Individual(const Individual& i);

};

#endif
