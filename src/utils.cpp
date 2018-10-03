#include "clustRviz.h"

bool is_nan(double x){
  return x != x;
}

int sgn(double x){
  int ret;
  if(x == 0){
    ret = 0;
  } else{
    ret = x/std::abs(x);
  }
  return(ret);
}
double soft_thresh(double x, double lambda){
  if(std::abs(x) < lambda){
    return 0;
  } else if (x > 0){
    return x - lambda;
  } else {
    return x + lambda;
  }
}

// Modify in place version for internal use
void MatrixProx(const Eigen::MatrixXd& X,
                Eigen::MatrixXd& V,
                double lambda,
                const Eigen::VectorXd& weights,
                bool l1 = true){

  Eigen::Index n = X.rows();
  Eigen::Index p = X.cols();

  if(l1){
    for(Eigen::Index i = 0; i < n; i++){
      for(Eigen::Index j = 0; j < p; j++){
        V(i, j) = soft_thresh(X(i, j), lambda * weights(i));
      }
    }
  } else {
    for(Eigen::Index i = 0; i < n; i++){
      Eigen::VectorXd X_i = X.row(i);
      double scale_factor = 1 - lambda * weights(i) / X_i.norm();

      if(scale_factor > 0){
        V.row(i) = X_i * scale_factor;
      } else {
        V.row(i).setZero();
      }
    }
  }
}

// Version for testing C++ code
// [[Rcpp::export]]
Eigen::MatrixXd MatrixProx(const Eigen::MatrixXd& X,
                           double lambda,
                           const Eigen::VectorXd& weights,
                           bool l1 = true){
  Eigen::Index n = X.rows();
  Eigen::Index p = X.cols();

  Eigen::MatrixXd V = Eigen::MatrixXd::Zero(n, p);

  MatrixProx(X, V, lambda, weights, l1);

  return V;
}
