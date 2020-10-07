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

std::complex<double> soft_thresh(const std::complex<double> x, double lambda){
  const double r = std::abs(x);
  if(r > lambda){
    return (r - lambda) * x / r;
  } else {
    return std::complex<double>(0.0, 0.0);
  }
}

// Some basic cheap checks that a weight
// matrix can lead to a connected graph
//
// Right now, the only check is that every observation
// has a connection (at some weight) to another observation
//
// [[Rcpp::export(rng = false)]]
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

// TROUT Alignment
// FIXME - Why doesn't this play nice with overloading? Prototype missing somewhwere?
Eigen::VectorXcd align_phase_v(const Eigen::VectorXcd& u,
                               const Eigen::VectorXcd& x){

  // TODO - Triple check this analytical result!
  std::complex<double> z_hat = u.conjugate().cwiseProduct(x).sum();
  if(std::abs(z_hat) > 0){
    z_hat /= std::abs(z_hat);
  }

  return u * z_hat;
}

// [[Rcpp::export(rng = false)]]
Eigen::MatrixXcd align_phase(const Eigen::MatrixXcd& U,
                             const Eigen::MatrixXcd& X){

  Eigen::Index n = X.rows();
  Eigen::Index p = X.cols();

  Eigen::MatrixXcd V(n, p);

  for(Eigen::Index i = 0; i < n; i++){
    V.row(i) = align_phase_v(U.row(i), X.row(i));
  }

  return V;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix trout_dist_impl(const Eigen::MatrixXcd& X){
  Eigen::Index n = X.rows();
  Rcpp::NumericMatrix distances(n, n);

  for(Eigen::Index i = 0; i < n; i++){
    for(Eigen::Index j = 0; j < n; j++){
      // FIXME - this calculates all distances twice - can simplify...

      // IMPORTANT: We have to transpose the result of align_phase_v to a _row vector_
      // or else it doesn't align with X.row() and silently discards all but the first
      // element
      distances(i, j) = std::sqrt((X.row(i) - align_phase_v(X.row(j), X.row(i)).transpose()).squaredNorm());
    }
  }

  return distances;
}

// TODO - Non-euclidean distances!
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cdist_impl(const Eigen::MatrixXcd& X){
  Eigen::Index n = X.rows();
  Rcpp::NumericMatrix distances(n, n);

  for(Eigen::Index i = 0; i < n; i++){
    for(Eigen::Index j = 0; j < n; j++){
      // FIXME - this calculates all distances twice - can simplify...

      // IMPORTANT: We have to transpose the result of align_phase_v to a _row vector_
      // or else it doesnt align with X.row() and silently discards all but the first
      // element
      distances(i, j) = (X.row(i) - X.row(j)).norm();
    }
  }
  return distances;
}
