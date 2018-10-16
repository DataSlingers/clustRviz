#include "clustRviz.h"

// [[Rcpp::export]]
Rcpp::List CARP_VIZcpp(const Eigen::MatrixXd& X,
                       const Eigen::MatrixXd& D,
                       double epsilon,
                       const Eigen::VectorXd& weights,
                       double rho              = 1,
                       int max_iter            = 10000,
                       int burn_in             = 50,
                       double back             = 0.5,
                       int keep                = 10,
                       int viz_max_inner_iter  = 15,
                       double viz_initial_step = 1.1,
                       double viz_small_step   = 1.01,
                       bool l1                 = false){

  ConvexClustering problem(X, D, weights, rho, l1);
  CARP_VIZ carp_viz(problem,
                    epsilon,
                    max_iter,
                    burn_in,
                    back,
                    keep,
                    viz_max_inner_iter,
                    viz_initial_step,
                    viz_small_step);

  return carp_viz.build_return_object();
}

// [[Rcpp::export]]
Rcpp::List CARPcpp(const Eigen::MatrixXd& X,
                   const Eigen::MatrixXd& D,
                   double epsilon,
                   double t,
                   const Eigen::VectorXd& weights,
                   double rho   = 1,
                   int max_iter = 10000,
                   int burn_in  = 50,
                   int keep     = 10,
                   bool l1      = false){

  ConvexClustering problem(X, D, weights, rho, l1);
  CARP carp(problem, epsilon, t, max_iter, burn_in, keep);

  return carp.build_return_object();
}
