#ifndef CLUSTRVIZ_OPTIM_POLICIES_H
#define CLUSTRVIZ_OPTIM_POLICIES_H 1

#include "clustRviz_base.h"
#include "clustRviz_logging.h"

template <class PROBLEM_TYPE>
class ADMMPolicy {
  // This is the full-solution ADMM policy
  //
  // It is pretty straightforward: just the standard ADMM + a check
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

      int k = 0;

      do {
        problem.save_old_values();
        problem.full_admm_step();
        iter++; k++;

        // problem.tick() will check for interrupts
        problem.tick(iter);

        if(k > 150){
          break; // Avoid infinite loops on a single gamma...
        }

      } while (!problem.admm_converged());

      ClustRVizLogger::info("ADMM converged with gamma = ") << problem.gamma << " after " << iter << " total iterations.";

      problem.store_values();
      problem.reset_aux();
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

template <class PROBLEM_TYPE>
class BackTrackingADMMPolicy {
  // This is the full-solution ADMM policy combined with VIZ-style back-tracking
public:
  BackTrackingADMMPolicy(PROBLEM_TYPE problem_,
                         const double epsilon_,
                         const int max_iter_,
                         const int burn_in_,
                         const double back_,
                         const int viz_max_inner_iter_,
                         const double viz_initial_step_,
                         const double viz_small_step_):

  problem(problem_),
  epsilon(epsilon_),
  max_iter(max_iter_),
  burn_in(burn_in_),
  back(back_),
  viz_max_inner_iter(viz_max_inner_iter_),
  viz_initial_step(viz_initial_step_),
  viz_small_step(viz_small_step_){};

  void solve(){
    // We need to keep an eye on gamma for back-tracking purposes
    double gamma     = epsilon;
    double gamma_old = epsilon;

    // The PROBLEM_TYPE constructor already stores the gamma = 0 solution,
    // so we start by setting gamma to gamma and performing an exact solve.
    problem.gamma = gamma;
    t = viz_initial_step;

    while( (iter < max_iter) & (!problem.is_complete()) ){
      ClustRVizLogger::info("Beginning iteration k = ") << iter + 1;
      ClustRVizLogger::info("gamma = ") << problem.gamma;

      int k = 0;

      // Save fusions to determine if we need to back-track
      problem.save_fusions();

      bool rep_iter = true;
      int try_iter = 0;
      double gamma_upper = gamma;
      double gamma_lower = gamma_old;

      while(rep_iter){

        // Run ADMM till convergence
        // Before running the ADMM, we need to reset the auxiliary variables each time.
        problem.gamma = gamma;
        problem.reset_aux();
        do {
          problem.save_old_values();
          problem.full_admm_step();
          iter++; k++;

          // problem.tick() will check for interrupts
          problem.tick(iter);

          if(k > 150){
            break; // Avoid infinite loops on a single gamma...
          }

        } while (!problem.admm_converged());

        ClustRVizLogger::info("ADMM converged with gamma = ") << problem.gamma << " after " << iter << " total iterations.";

        try_iter++;
        if(try_iter > viz_max_inner_iter){
          break;
        }

        // Now that we've converged, we begin the back-tracking (VIZ) checks
        // This is the core "*-VIZ" logic
        //
        // After running the ADMM, instead of proceeding with a fixed step-size update,
        // we include a back-tracking step
        if( (!problem.is_interesting_iter()) & (try_iter == 1) ){
          // If the sparsity pattern (number of fusions) hasn't changed, we have
          // no need to back-track (we didn't miss any fusions) so we can go immediately
          // to the next iteration.
          rep_iter = false;
          ClustRVizLogger::info("No fusions identified -- continuing to next step.");
        } else if(problem.multiple_fusions()){
          // If we see two (or more) new fusions, we need to back-track and figure
          // out which one occured first
          problem.load_old_fusions();
          if(try_iter == 1){
            gamma = 0.5 * (gamma_lower + gamma_upper);
          } else{
            gamma_upper = gamma;
            gamma = 0.5 * (gamma_lower + gamma_upper);
          }
          ClustRVizLogger::info("Too many fusions -- backtracking.");
        } else if(!problem.is_interesting_iter()){
          // If we don't observe any new fusions, we move our regularization level
          // up to try to find one
          problem.load_old_fusions();
          gamma_lower = gamma;
          gamma = 0.5 * (gamma_lower + gamma_upper);
          ClustRVizLogger::info("Fusion not isolated -- moving forward.");
        } else {
          // If we see exactly one new fusion, we have a good step size and exit
          // the inner back-tracking loop
          rep_iter = false;
          ClustRVizLogger::info("Good iteration - continuing to next step.");
        }
      }

      // We always store the final iterate of a back-tracking search
      problem.store_values();

      gamma_old = gamma; // Save this for future back-tracking iterations

      // If we have gotten to the "lots of fusions" part of the solution space, start
      // taking smaller step sizes.
      if(problem.has_fusions()){
        t = viz_small_step;
      }

      gamma *= t;
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
  const int max_iter = 10000;
  const int burn_in  = 50;
  const double back  = 0.5;
  const int viz_max_inner_iter  = 15;
  const double viz_initial_step = 1.1;
  const double viz_small_step   = 1.01;

  // Algorithm state
  double t;
  int iter = 0;
  bool solved = false;
};

#endif
