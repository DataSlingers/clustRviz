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
                       bool l1                 = false,
                       bool show_progress      = true){

  ConvexClustering problem(X, D, weights, rho, l1, show_progress);
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
                   bool l1      = false,
                   bool show_progress = true){

  ConvexClustering problem(X, D, weights, rho, l1, show_progress);
  CARP carp(problem, epsilon, t, max_iter, burn_in, keep);

  return carp.build_return_object();
}

// [[Rcpp::export]]
Rcpp::List ConvexClusteringADMMcpp(const Eigen::MatrixXd& X,
                                   const Eigen::MatrixXd& D,
                                   double epsilon,
                                   double t,
                                   const Eigen::VectorXd& weights,
                                   double rho   = 1,
                                   int max_iter = 10000,
                                   bool l1      = false,
                                   bool show_progress = true){

  ConvexClustering problem(X, D, weights, rho, l1, show_progress);
  ConvexClusteringADMM admm(problem, epsilon, t, max_iter);

  return admm.build_return_object();
}

// [[Rcpp::export]]
Rcpp::List CBASS_VIZcpp(const Eigen::MatrixXd& X,
                        const Eigen::MatrixXd& D_row,
                        const Eigen::MatrixXd& D_col,
                        double epsilon,
                        const Eigen::VectorXd& weights_col,
                        const Eigen::VectorXd& weights_row,
                        double rho              = 1,
                        int max_iter            = 10000,
                        int burn_in             = 50,
                        double back             = 0.5,
                        int keep                = 10,
                        int viz_max_inner_iter  = 15,
                        double viz_initial_step = 1.1,
                        double viz_small_step   = 1.01,
                        bool l1                 = false,
                        bool show_progress      = true){

  ConvexBiClustering problem(X, D_row, D_col, weights_row, weights_col, rho, l1, show_progress);
  CBASS_VIZ cbass_viz(problem,
                      epsilon,
                      max_iter,
                      burn_in,
                      back,
                      keep,
                      viz_max_inner_iter,
                      viz_initial_step,
                      viz_small_step);

  return cbass_viz.build_return_object();
}

// [[Rcpp::export]]
Rcpp::List CBASScpp(const Eigen::MatrixXd& X,
                    const Eigen::MatrixXd& D_row,
                    const Eigen::MatrixXd& D_col,
                    double epsilon,
                    double t,
                    const Eigen::VectorXd& weights_row,
                    const Eigen::VectorXd& weights_col,
                    double rho   = 1,
                    int max_iter = 1e4,
                    int burn_in  = 50,
                    int keep     = 10,
                    bool l1      = false,
                    bool show_progress = true){

  ConvexBiClustering problem(X, D_row, D_col, weights_row, weights_col, rho, l1, show_progress);
  CBASS cbass(problem, epsilon, t, max_iter, burn_in, keep);

  return cbass.build_return_object();
}
