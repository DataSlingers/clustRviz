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

// Apply a row-wise prox operator (with weights) to a matrix
// [[Rcpp::export]]
Eigen::MatrixXd MatrixProx(const Eigen::MatrixXd& X,
                           double lambda,
                           const Eigen::VectorXd& weights,
                           bool l1 = true){
  Eigen::Index n = X.rows();
  Eigen::Index p = X.cols();

  Eigen::MatrixXd V(n, p);

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

  return V;
}

// Some basic cheap checks that a weight
// matrix can lead to a connected graph
//
// Right now, the only check is that every observation
// has a connection (at some weight) to another observation
//
// [[Rcpp::export]]
void check_weight_matrix(const Eigen::MatrixXd& weight_matrix){
  Eigen::Index n = weight_matrix.rows();
  Eigen::Index p = weight_matrix.cols();

  if(n != p){
    ClustRVizLogger::error("Clustering weight matrix is not square.");
  }

  // Start at i = 1, so we don't check first corner of triangular matrix
  // since there are no neighbors considered
  for(int i = 1; i < n; i++){
    bool neighbor_found = false;

    for(int j = 0; j < i; j++){
      if (weight_matrix(i, j) != 0){
        neighbor_found = true;
        break;
      }
    }

    if (!neighbor_found) {
      ClustRVizLogger::error("No neighbor found for observation ") << i + 1 <<
         " -- convex clustering cannot succeed. You may need to rescale your data.";
    }
  }
}

// Tensor projection along the second mode
//
// Given a 3D tensor X in R^{n-by-p-by-q} (observations by features by iterations)
// and a rotation matrix Y in R^{p-by-k} (features by principal components), we
// want to get a projected array in R^{n-by-k-by-q} giving the path of the principal
// components
//
// This is straightforward, but "loopy" so we implement it in Rcpp / RcppEigen for speed
// [[Rcpp::export]]
Rcpp::NumericVector tensor_projection(Rcpp::NumericVector X, const Eigen::MatrixXd& Y){

  // Validate X
  Rcpp::IntegerVector X_dims = X.attr("dim");
  if(X_dims.size() != 3){
    ClustRVizLogger::error("X must be a three rank tensor.");
  }
  int n = X_dims(0);
  int p = X_dims(1);
  int q = X_dims(2);

  // Validate Y
  if(Y.rows() != p){
    ClustRVizLogger::error("The dimensions of X and Y do not match -- ") << p << " != " << Y.rows();
  }

  int k = Y.cols();

  Rcpp::NumericVector result(n * k * q);
  Rcpp::IntegerVector result_dims{n, k, q};
  result.attr("dim") = result_dims;

  for(int i = 0; i < q; i++){
    Eigen::MatrixXd X_slice = Eigen::Map<Eigen::MatrixXd>(&X[n * p * i], n, p);
    Eigen::MatrixXd X_slice_projected = X_slice * Y;
    Eigen::Map<Eigen::MatrixXd>(&result[n * k * i], n, k) = X_slice_projected;
  }

  return result;
}
