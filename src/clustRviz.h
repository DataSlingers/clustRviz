#include <RcppEigen.h>
#include <vector>
#include <set>
#include "clustRviz_logging.h"

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

// Helper to determine if STL set contains an element
//
// In general, this is not efficient because one wants to do something
// with the element and/or its location, but here we really only need containment
template <typename T>
bool contains(const std::set<T>& container, T element){
  typename std::set<T>::const_iterator it = container.find(element);
  return it != container.end();
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

// Prototypes - utils.cpp
void MatrixProx(const Eigen::MatrixXd&,
                Eigen::MatrixXd&,
                double,
                const Eigen::VectorXd&,
                bool);
