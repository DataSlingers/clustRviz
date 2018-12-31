#include "clustRviz.h"

// [[Rcpp::export]]
Rcpp::List CARPcpp(const Eigen::MatrixXd& X,
                   const Eigen::MatrixXd& D,
                   const Eigen::VectorXd& weights,
                   double epsilon,
                   double t,
                   double rho              = 1,
                   int max_iter            = 10000,
                   int burn_in             = 50,
                   double back             = 0.5,
                   int keep                = 10,
                   int viz_max_inner_iter  = 15,
                   double viz_initial_step = 1.1,
                   double viz_small_step   = 1.01,
                   bool l1                 = false,
                   bool show_progress      = true,
                   bool back_track         = false,
                   bool exact              = false){

  ConvexClustering problem(X, D, weights, rho, l1, show_progress);

  if(exact){
    if(back_track){
      ConvexClusteringADMM_VIZ admm_viz(problem,
                                        epsilon,
                                        max_iter,
                                        burn_in,
                                        back,
                                        viz_max_inner_iter,
                                        viz_initial_step,
                                        viz_small_step);

      return admm_viz.build_return_object();
    } else {
      ConvexClusteringADMM admm(problem, epsilon, t, max_iter);
      return admm.build_return_object();
    }
  } else {
    if(back_track){
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

    CARP carp(problem, epsilon, t, max_iter, burn_in, keep);
    return carp.build_return_object();
  }
}

// [[Rcpp::export]]
Rcpp::List CBASScpp(const Eigen::MatrixXd& X,
                    const Eigen::MatrixXd& D_row,
                    const Eigen::MatrixXd& D_col,
                    const Eigen::VectorXd& weights_col,
                    const Eigen::VectorXd& weights_row,
                    double epsilon,
                    double t,
                    double rho              = 1,
                    int max_iter            = 10000,
                    int burn_in             = 50,
                    double back             = 0.5,
                    int keep                = 10,
                    int viz_max_inner_iter  = 15,
                    double viz_initial_step = 1.1,
                    double viz_small_step   = 1.01,
                    bool l1                 = false,
                    bool show_progress      = true,
                    bool back_track         = false,
                    bool exact              = false){

  ConvexBiClustering problem(X, D_row, D_col, weights_row, weights_col, rho, l1, show_progress);

  if(exact){
    if(back_track){
      ConvexBiClusteringADMM_VIZ admm_viz(problem,
                                          epsilon,
                                          max_iter,
                                          burn_in,
                                          back,
                                          viz_max_inner_iter,
                                          viz_initial_step,
                                          viz_small_step);

      return admm_viz.build_return_object();
    } else {
      ConvexBiClusteringADMM admm(problem, epsilon, t, max_iter);
      return admm.build_return_object();
    }
  } else {
    if(back_track){
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

    CBASS cbass(problem, epsilon, t, max_iter, burn_in, keep);
    return cbass.build_return_object();
  }
}
