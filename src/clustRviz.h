#define ARMA_USE_SUPERLU 0

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#define ARMA_64BIT_WORD

#define CLUSTRVIZ_CHECK_USER_INTERRUPT_RATE 50

// From https://stackoverflow.com/a/26268017
template <typename T_full, typename T_ind>
T_full extract(const T_full& full, const T_ind& ind){
  Eigen::Index num_indices = ind.innerSize();
  T_full target(num_indices);

  for(Eigen::Index i = 0; i < num_indices; i++){
    target[i] = full[ind[i]];
  }

  return target;
}

// Prototypes - arma implementations
arma::vec restride(const arma::vec&, arma::uword k);
double TwoNorm(const arma::colvec&);
arma::vec cv_sparse_solve(const Eigen::SparseMatrix<double>&, const arma::vec&);
arma::colvec DMatOpv2(const arma::colvec&, int, const arma::umat&,
                      const arma::umat&, const arma::umat&);
arma::colvec DtMatOpv2(const arma::colvec&, int, int,
                       const arma::umat&, const arma::umat&, const arma::umat&);
arma::colvec ProxL2(const arma::colvec&, int, const arma::colvec&, const arma::umat&);
arma::colvec ProxL1(const arma::colvec&, int, double, const arma::colvec& weights);

// Prototypes - Eigen implementations
Eigen::VectorXd restride(const Eigen::VectorXd&, arma::uword k);
double TwoNorm(const Eigen::VectorXd&);
Eigen::VectorXd cv_sparse_solve(const Eigen::SparseMatrix<double>&,
                                const Eigen::VectorXd&);
Eigen::VectorXd DMatOpv2(const Eigen::VectorXd&, int, const Eigen::MatrixXi&,
                         const Eigen::MatrixXi&, const Eigen::MatrixXi&);
Eigen::VectorXd DtMatOpv2(const Eigen::VectorXd&, int, int,
                          const Eigen::MatrixXi&, const Eigen::MatrixXi&,
                          const Eigen::MatrixXi&);
Eigen::VectorXd ProxL2(const Eigen::VectorXd&, int,
                       const Eigen::VectorXd&, const Eigen::MatrixXi&);
Eigen::VectorXd ProxL1(const Eigen::VectorXd&, int, double,
                       const Eigen::VectorXd& weights);
