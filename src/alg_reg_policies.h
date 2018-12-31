#ifndef CLUSTRVIZ_ALG_REG_POLICIES_H
#define CLUSTRVIZ_ALG_REG_POLICIES_H 1

#include "clustRviz_base.h"
#include "clustRviz_logging.h"

template <class PROBLEM_TYPE>
class AlgorithmicRegularizationFixedStepSizePolicy {
  // This is the fixed step-size version of CARP/CBASS
  //
  // Internally it is pretty straightforward, just the standard ADMM + multiplicative
  // update stuff combined with a bit of book-keeping to store values on "interesting"
  // iterations (i.e., those with a fusion)
public:
  AlgorithmicRegularizationFixedStepSizePolicy(PROBLEM_TYPE problem_,
                                               const double epsilon_,
                                               const double t_,
                                               const int max_iter_,
                                               const int burn_in_,
                                               const int keep_):
  problem(problem_),
  epsilon(epsilon_),
  t(t_),
  max_iter(max_iter_),
  burn_in(burn_in_),
  keep(keep_){};

  void solve(){
    // The PROBLEM_TYPE constructor already stores the gamma = 0 solution,
    // so we start by setting epsilon to gamma and beginning a solve
    problem.gamma = epsilon;

    while( (iter < max_iter) & (!problem.is_complete()) ){
      ClustRVizLogger::info("Beginning iteration k = ") << iter + 1;
      ClustRVizLogger::debug("gamma = ") << problem.gamma;

      problem.save_fusions();
      problem.admm_step();

      // Store interesting iterations, but otherwise ignore the burn-in phase
      if( problem.is_interesting_iter() | ((iter % keep == 0) & (iter > burn_in)) ){
        problem.store_values();
      }

      // The progress bar class also checks for user interrupts on ticks
      iter++;
      problem.tick(iter);

      if (iter >= burn_in) {
        problem.gamma *= t;
      }
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
  const int burn_in  = 50;
  const int keep     = 10;

  // Algorithm state
  int iter = 0;
  bool solved = false;
};

template <class PROBLEM_TYPE>
class AlgorithmicRegularizationBacktrackingPolicy {
  // This is the back-tracking version of CARP/CBASS
  //
  // Internally it is much more complicated than the fixed step-size versions
  // and requires us to repeatedly copy around the optimization variables since
  // the ADMM step only looks at "V" (and we can't set it to look at V_old)
  //
  // See comments below
public:
  AlgorithmicRegularizationBacktrackingPolicy(PROBLEM_TYPE problem_,
                                              const double epsilon_,
                                              const int max_iter_,
                                              const int burn_in_,
                                              const double back_,
                                              const int keep_,
                                              const int viz_max_inner_iter_,
                                              const double viz_initial_step_,
                                              const double viz_small_step_):
  problem(problem_),
  epsilon(epsilon_),
  max_iter(max_iter_),
  burn_in(burn_in_),
  back(back_),
  keep(keep_),
  viz_max_inner_iter(viz_max_inner_iter_),
  viz_initial_step(viz_initial_step_),
  viz_small_step(viz_small_step_){};

  void solve(){
    // We need to keep an eye on gamma for back-tracking purposes
    double gamma     = epsilon;
    double gamma_old = epsilon;

    // The PROBLEM_TYPE constructor already stores the gamma = 0 solution,
    // so we start by setting epsilon to gamma and beginning a solve
    problem.gamma = gamma;
    t = viz_initial_step;

    while( (iter < max_iter) & (!problem.is_complete()) ){
      ClustRVizLogger::info("Beginning iteration k = ") << iter + 1;
      ClustRVizLogger::info("gamma = ") << problem.gamma;

      // Pre-load V_old, Z_old, etc. so we have them for the 'load_old_variables' step below
      problem.save_fusions();
      problem.save_old_values();

      bool rep_iter = true;
      int try_iter = 0;
      double gamma_upper = gamma;
      double gamma_lower = gamma_old;

      while(rep_iter){
        // Re-load old copies to have a true "back-track" instead of a refinement
        problem.load_old_variables();
        problem.gamma = gamma;
        problem.admm_step();

        try_iter++;

        if(try_iter > viz_max_inner_iter){
          break;
        }

        // This is the core "*-VIZ" logic
        //
        // After one iteration of each ADMM step, instead of
        // proceeding with a fixed step-size update, we include a back-tracking step
        // (different conditions described in more detail below)
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

        // The progress bar class also checks for user interrupts on ticks
        // so we do this inside of back-tracking loops so we don't get stuck
        // uninterruptable on very slow calculations
        problem.tick(iter);
      }

      gamma_old = gamma; // Save this for future back-tracking iterations

      // If we have gotten to the "lots of fusions" part of the solution space, start
      // taking smaller step sizes.
      if(problem.has_fusions()){
        t = viz_small_step;
      }

      // If we have seen a fusion or are otherwise interested in keeping this iteration,
      // add values to our storage buffers
      if( problem.is_interesting_iter() | ((iter % keep == 0) & (iter > burn_in)) ){
        problem.store_values();
      }

      iter++;

      if (iter > burn_in) {
        gamma *= t;
      }
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
  const int keep     = 10;
  const int viz_max_inner_iter  = 15;
  const double viz_initial_step = 1.1;
  const double viz_small_step   = 1.01;

  // Algorithm state
  double t;
  int iter = 0;
  bool solved = false;
};

#endif
