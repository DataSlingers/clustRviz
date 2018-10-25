#ifndef CLUSTRVIZ_OPTIM_POLICIES_H
#define CLUSTRVIZ_OPTIM_POLICIES_H 1

#include "clustRviz_base.h"
#include "clustRviz_logging.h"

template <class PROBLEM_TYPE>
class ADMMPolicy {
  // This is the full-solution ADMM policy
  //
  // It is pretty straightforwardL just the standard ADMM + a check
  // (based on the norm of the difference in the V and Z variables) for convergence
  //
  // We store the result of each level of the regularization parameter
public:
  ADMMPolicy(PROBLEM_TYPE problem_,
            const double epsilon_,
            const double t_,
            const int max_iter_):

  problem(problem_),
  epsilon(epsilon_),
  t(t_),
  max_iter(max_iter_){};

  void solve(){
    // The PROBLEM_TYPE constructor already stores the gamma = 0 solution,
    // so we start by setting epsilon to gamma and beginning a solve
    problem.gamma = epsilon;

    while( (iter < max_iter) & (!problem.is_complete()) ){
      ClustRVizLogger::info("Starting ADMM with gamma = ") << problem.gamma;

      do {
        problem.save_old_values();
        problem.admm_step();
        iter++;
        problem.tick(iter);

      } while (!problem.admm_converged());

      ClustRVizLogger::info("ADMM converged with gamma = ") << problem.gamma << " after " << iter << " total iterations.";

      problem.store_values();
      problem.gamma *= t;
    }

    if(iter >= max_iter){
      ClustRVizLogger::warning("Clustering ended early -- `max_iter` reached. Treat results with caution.");
    }

    solved = true;
  }

  Rcpp::List build_return_object(){
    if(!solved) solve();

    return problem.build_return_object();
  }
private:
  PROBLEM_TYPE problem;
  const double epsilon;
  const double t;
  const int max_iter = 10000;

  // Algorithm state
  int iter = 0;
  bool solved = false;
};

#endif
