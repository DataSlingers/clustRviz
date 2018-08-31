#include <RcppEigen.h>

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

// Prototypes - Eigen implementations
Eigen::VectorXd restride(const Eigen::VectorXd&, Eigen::Index);
Eigen::VectorXd DMatOpv2(const Eigen::VectorXd&, int, const Eigen::MatrixXi&,
                         const Eigen::MatrixXi&, const Eigen::MatrixXi&);
Eigen::VectorXd DtMatOpv2(const Eigen::VectorXd&, int, int,
                          const Eigen::MatrixXi&, const Eigen::MatrixXi&,
                          const Eigen::MatrixXi&);
Eigen::VectorXd ProxL2(const Eigen::VectorXd&, int,
                       const Eigen::VectorXd&, const Eigen::MatrixXi&);
Eigen::VectorXd ProxL1(const Eigen::VectorXd&, int, double,
                       const Eigen::VectorXd& weights);
