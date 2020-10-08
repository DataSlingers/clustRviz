#ifndef CLUSTRVIZ_BASE_H
#define CLUSTRVIZ_BASE_H 1

#include <RcppEigen.h>
#include <vector>
#include <set>

#define CLUSTRVIZ_STATUS_UPDATE_TIME_SECS 0.1  // Print status to screen every 0.1s
#define CLUSTRVIZ_STATUS_WIDTH_CHECK 20        // Every 20 status updates * 0.1s => every 2s
#define CLUSTRVIZ_DEFAULT_STOP_PRECISION 1e-10 //Stop when cellwise diff between iters < val

// Prototypes - utils.cpp
double soft_thresh(double, double);
std::complex<double> soft_thresh(const std::complex<double>, double);

// Helper to determine if STL set contains an element
//
// In general, this is not efficient because one wants to do something
// with the element and/or its location, but here we really only need containment
template <typename T>
bool contains(const std::set<T>& container, T element){
  typename std::set<T>::const_iterator it = container.find(element);
  return it != container.end();
};

template <typename DataType>
double scaled_squared_norm(const Eigen::MatrixBase<DataType>& X){
  return X.squaredNorm() / X.size();
}

template <typename DataType>
Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic> MatrixRowProx(const Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic>& X,
                                                                      double lambda,
                                                                      const Eigen::VectorXd& weights,
                                                                      bool l1 = true){

  using MatrixXt = Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic>;
  using VectorXt = Eigen::Matrix<DataType, Eigen::Dynamic, 1>;

  Eigen::Index n = X.rows();
  Eigen::Index p = X.cols();

  MatrixXt V(n, p);

  if(l1){
    for(Eigen::Index i = 0; i < n; i++){
      for(Eigen::Index j = 0; j < p; j++){
        V(i, j) = soft_thresh(X(i, j), lambda * weights(i));
      }
    }
  } else {
    for(Eigen::Index i = 0; i < n; i++){
      VectorXt X_i = X.row(i);
      double scale_factor = 1 - lambda * weights(i) / X_i.norm();

      if(scale_factor > 0){
        V.row(i) = X_i * scale_factor;
      } else {
        V.row(i).setZero();
      }
    }
  }

  return V;
}

template <typename DataType>
Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic> MatrixColProx(const Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic>& X,
                                                                      double lambda,
                                                                      const Eigen::VectorXd& weights,
                                                                      bool l1 = true){

  using MatrixXt = Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic>;
  using VectorXt = Eigen::Matrix<DataType, Eigen::Dynamic, 1>;

  Eigen::Index n = X.rows();
  Eigen::Index p = X.cols();

  MatrixXt V(n, p);

  if(l1){
    for(Eigen::Index i = 0; i < n; i++){
      for(Eigen::Index j = 0; j < p; j++){
        V(i, j) = soft_thresh(X(i, j), lambda * weights(j));
      }
    }
  } else {
    for(Eigen::Index j = 0; j < p; j++){
      VectorXt X_j = X.col(j);
      double scale_factor = 1 - lambda * weights(j) / X_j.norm();

      if(scale_factor > 0){
        V.col(j) = X_j * scale_factor;
      } else {
        V.col(j).setZero();
      }
    }
  }

  return V;
}

#endif
