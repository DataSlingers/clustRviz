#ifndef CLUSTRVIZ_UTILS_H
#define CLUSTRVIZ_UTILS_H 1
// This header defines template versions of complex utilities
// These are templated on the data type to allow for both real and complex data
// These are not in clustRviz_base.h since some of them depend on ClustRVizLogger...

#include <RcppEigen.h>
#include "clustRviz_logging.h"

// U-smoothing for convex clustering
//
// Given cluster memberships, replace rows of U which belong to the same cluster
// with their mutual mean....
template <typename RcppVector, typename DataType>
RcppVector smooth_u_clustering_impl(RcppVector U_old, Rcpp::List cluster_info_list){

  using MatrixXt = Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic>;
  using VectorXt = Eigen::Matrix<DataType, Eigen::Dynamic, 1>;

  // The first argument is really an array but we pass as a NumericVector
  // The second argument is a list produced by get_cluster_assignments()
  Rcpp::IntegerVector U_dims = U_old.attr("dim");
  if(U_dims.size() != 3){
    ClustRVizLogger::error("U must be a three rank tensor.");
  }
  int N = U_dims(0);
  int P = U_dims(1);
  int Q = U_dims(2);

  // Check length of cluster_info
  if(cluster_info_list.size() != Q){
    ClustRVizLogger::error("Dimensions of U and cluster_info do not match");
  }

  RcppVector U(N * P * Q);
  U.attr("dim") = U_dims;
  Rcpp::rownames(U) = Rcpp::rownames(U_old);
  Rcpp::colnames(U) = Rcpp::colnames(U_old);

  for(int q = 0; q < Q; q++){
    Rcpp::List cluster_info = cluster_info_list[q];
    int n_clusters = Rcpp::as<int>(cluster_info[2]);

    Rcpp::IntegerVector cluster_ids   = cluster_info[0];
    Rcpp::IntegerVector cluster_sizes = cluster_info[1];

    // There's a lot going on on the RHS here, so let's un-pack (inside outwards)
    // First, we get a pointer to the relevant slice of U_old
    //   This works because the RcppVector is a C++ wrapper around a SEXP which
    //   is ultimately just a pointer to the relevant memory
    // We when cast it to an appropriate C++ pointer type
    //   For real data, this is a no-op since both R and Eigen use doubles for real data
    //   For complex data, this matters because we convert from R's homegrown Rcomplex*
    //   to a std::complex* pointer as eigen expects
    // We then use Eigen::Map<Eigen::Matrix<DataType>> to get an Eigen::Matrix<DataType> backed
    //   by R's memory in a read only fashion.
    // The same construct is used below to load the smoothed data into U
    MatrixXt U_old_slice = Eigen::Map<MatrixXt>(reinterpret_cast<DataType*>(&U_old[N * P * q]), N, P);
    MatrixXt U_new(N, P);

    for(int j = 1; j <= n_clusters; j++){ // Cluster IDs are 1-based (per R conventions)
      VectorXt vec(P); vec.setZero();

      // Manually work out new mean
      for(int n = 0; n < N; n++){
        if(cluster_ids[n] == j){
          vec += U_old_slice.row(n);
        }
      }

      vec /= cluster_sizes[j - 1]; // Subtract 1 to adjust to C++ indexing

      // Assign new mean where needed...
      for(int n = 0; n < N; n++){
        if(cluster_ids[n] == j){
          U_new.row(n) = vec;
        }
      }
    }

    Eigen::Map<MatrixXt>(reinterpret_cast<DataType*>(&U[N * P * q]), N, P) = U_new;
  }

  return U;
}

// Tensor projection along the second mode
//
// Given a 3D tensor X in F^{n-by-p-by-q} (observations by features by iterations)
// and a rotation matrix Y in F^{p-by-k} (features by principal components), we
// want to get a projected array in F^{n-by-k-by-q} giving the path of the principal
// components
//
// This is straightforward, but "loopy" so we implement it in Rcpp / RcppEigen for speed
// We use some template magic to support F = R (real) and F = C (complex) data
template <typename RcppVector, typename DataType>
RcppVector tensor_projection_impl(RcppVector X, const Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic>& Y){

  using MatrixXt = Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic>;

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

  RcppVector result(n * k * q);
  Rcpp::IntegerVector result_dims{n, k, q};
  result.attr("dim") = result_dims;

  for(int i = 0; i < q; i++){
    // There's a lot going on on the RHS here, so let's un-pack (inside outwards)
    // First, we get a pointer to the relevant slice of X
    //   This works because the RcppVector is a C++ wrapper around a SEXP which
    //   is ultimately just a pointer to the relevant memory
    // We when cast it to an appropriate C++ pointer type
    //   For real data, this is a no-op since both R and Eigen use doubles for real data
    //   For complex data, this matters because we convert from R's homegrown Rcomplex*
    //   to a std::complex* pointer as eigen expects
    // We then use Eigen::Map<Eigen::Matrix<DataType>> to get an Eigen::Matrix<DataType> backed
    //   by R's memory in a read only fashion.
    // The same construct is used below to load the smoothed data into result
    MatrixXt X_slice = Eigen::Map<MatrixXt>(reinterpret_cast<DataType*>(&X[n * p * i]), n, p);
    MatrixXt X_slice_projected = X_slice * Y;
    Eigen::Map<MatrixXt>(reinterpret_cast<DataType*>(&result[n * k * i]), n, k) = X_slice_projected;
  }

  return result;
}

#endif
