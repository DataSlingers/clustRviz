#include "clustRviz.h"

// TODO - Consolidate CARP and CARP-VIZ
// Most of the internal logic is the same (modulo back-tracking vs fixed step size)

// [[Rcpp::export]]
Rcpp::List CARP_VIZ(const Eigen::VectorXd& x,
                    int n,
                    int p,
                    double lambda_init, // TODO: Change to gamma_init
                    const Eigen::VectorXd& weights,
                    const Eigen::VectorXd& uinit, // TODO: Change to u_init
                    const Eigen::VectorXd& vinit, // TODO: Change to v_init
                    const Eigen::SparseMatrix<double>& premat,
                    const Eigen::MatrixXi& IndMat,
                    const Eigen::MatrixXi& EOneIndMat,
                    const Eigen::MatrixXi& ETwoIndMat,
                    double rho      = 1,
                    int max_iter    = 10000,
                    int burn_in     = 50,
                    bool verbose    = false,
                    double back     = 0.5,
                    int ti          = 15,
                    double t_switch = 1.01,
                    int keep        = 10,
                    bool l1         = false){

  // Typically, our weights are "sparse" (i.e., mostly zeros) because we
  // drop small weights to achieve performance.
  Eigen::Index num_edges = EOneIndMat.rows();

  /// Set-up storage for CARP iterates

  // In order to pre-allocate storage arrays, we need to estimate the number of
  // steps with fusions we will encounter. Dendrograms are the common case and
  // we are only clustering observations, so we expect O(n) fusions. It's a bit
  // cheaper to drop observations than to extend the internal buffers of our
  // storage objects, so we use 1.5n for now
  Eigen::Index buffer_size = 1.5 * n;

  // Primal variable (corresponds to u in the notation of Chi & Lange (JCGS, 2015))
  Eigen::VectorXd u_new(n * p);              // Working copies
  Eigen::VectorXd u_old(n * p);
  u_new = uinit;
  Eigen::MatrixXd UPath(n * p, buffer_size); // Storage (for values to return to R)
  UPath.col(0) = uinit;


  // 'Split' variable (corresponds to v in the notation of Chi & Lange (JCGS, 2015))
  Eigen::VectorXd v_new(p * num_edges);              // Working copies
  Eigen::VectorXd v_old(p * num_edges);
  v_new = vinit;
  Eigen::MatrixXd VPath(p * num_edges, buffer_size); // Storage (for values to return to R)
  VPath.col(0) = vinit;

  // (Scaled) dual variable (corresponds to lambda in the notation of Chi and Lange (JCGS, 2015))
  Eigen::VectorXd z_new(p * num_edges); // Working copy
  Eigen::VectorXd z_old(p * num_edges); // No storage needed since these aren't of direct interest
  z_new = v_new;

  // Regularization level
  double gamma     = lambda_init;          // Working copies
  double gamma_old = lambda_init;
  Eigen::VectorXd gamma_path(buffer_size); // Storage (to be returned to R)
  gamma_path(0) = lambda_init;

  // Fusions -- TODO: Confirm the semantics of these objects with JN
  Eigen::VectorXd vZeroIndsnew = Eigen::VectorXd::Zero(num_edges);          // Working copies
  Eigen::VectorXd vZeroIndsold(num_edges);                                  // (we begin with no fusions)
  Eigen::MatrixXd vZeroInds_Path = Eigen::MatrixXd(num_edges, buffer_size); // Storage (to be returned to R)

  /// END Preallocations

  // At each iteration, we need to calculate A^{-1}B_k for some (sparse) A
  // This is a relatively expensive iteration, but the core cost is a sparse LU
  // factorization of A which can be amortized over iterations so we pre-compute it here
  Eigen::SparseLU<Eigen::SparseMatrix<double> > premat_solver;
  premat_solver.compute(premat);

  // Book-keeping variables
  // Number of iterations stored, total iteration count, number of fusions
  Eigen::Index path_iter  = 0;
  Eigen::Index iter       = 0;
  Eigen::Index nzeros_old = 0;
  Eigen::Index nzeros_new = 0;


  // We begin CARP-VIZ by taking relatively large step sizes (t = 1.1)
  // but once we get to an "interesting" part of the path, we switch to
  // smaller step sizes (as determined by t_switch)
  double t = 1.1;

  while( (iter < max_iter) & (nzeros_new < num_edges) ){
    // Begin iteration - move updated values to "_old" values
    //
    // TODO -- Confirm that these use "move semantics" to avoid full copies
    u_old        = u_new;
    v_old        = v_new;
    z_old        = z_new;
    vZeroIndsold = vZeroIndsnew;
    nzeros_old   = nzeros_new;

    bool rep_iter = true;      // Does this iteration need to be repeated?
    Eigen::Index try_iter = 0; // How many times has this iteration already been repeated

    double gamma_upper = gamma;
    double gamma_lower = gamma_old;

    // This is the core CARP-VIZ Logic:
    //
    // We take CARP-type (one iteration of each ADMM step) steps, but instead of
    // proceeding with a fixed step-size update, we include a back-tracking step
    // (described in more detail below)
    //
    //
    while(rep_iter){
      // U-update
      // TODO - Document what is happening here
      Eigen::VectorXd solver_input = DtMatOpv2(rho * v_old - z_old, n, p, IndMat, EOneIndMat, ETwoIndMat);
      solver_input += x;
      solver_input /= rho;
      u_new = premat_solver.solve(solver_input);

      // V-update
      // TODO - Document what is happening here
      Eigen::VectorXd prox_argument = DMatOpv2(u_new, p, IndMat, EOneIndMat, ETwoIndMat) + (1/rho)*z_old;

      if(l1){
        v_new = ProxL1(prox_argument, p, (1/rho) * gamma, weights);
      } else {
        v_new = ProxL2(prox_argument, p, (1/rho) * weights * gamma, IndMat);
      }

      // Z-update
      // TODO - Document what is happening here
      z_new = z_old + rho*(DMatOpv2(u_new, p, IndMat, EOneIndMat, ETwoIndMat) - v_new);

      // TODO- Document this check: what are we trying to do here?
      for(int l = 0; l < num_edges; l++){
        Eigen::VectorXi v_index = IndMat.row(l);
        if(extract(v_new, v_index).sum() == 0){
          vZeroIndsnew(l) = 1;
        }
      }
      nzeros_new = vZeroIndsnew.sum();

      try_iter++; // Increment internal iteration count (used to check for stopping below)

      if( (nzeros_new == nzeros_old) & (try_iter == 1) ){
        // If the sparsity pattern (number of fusions) hasn't changed, we have
        // no need to back-track (we didn't miss any fusions) so we can go immediately
        // to the next iteration.
        rep_iter = false;
      } else if(nzeros_new > nzeros_old + 1){
        // If we see two (or more) new fusions, we need to back-track and figure
        // out which one occured first
        vZeroIndsnew = vZeroIndsold;
        if(try_iter == 1){
          gamma = 0.5 * (gamma_lower + gamma_upper);
        } else{
          gamma_upper = gamma;
          gamma = 0.5 * (gamma_lower + gamma_upper);
        }
      } else if(nzeros_new == nzeros_old){
        // If we don't observe any new fusions, we take another iteration without
        vZeroIndsnew = vZeroIndsold;
        gamma_lower = gamma;
        gamma = 0.5 * (gamma_lower + gamma_upper);
      } else{
        // If we see exactly one new fusion, we have a good step size and exit
        // the inner back-tracking loop
        rep_iter = false;
      }

      // Safety check - only so many iterations of inner loop before we move on
      if(try_iter > ti){
        rep_iter = false;
      }
    }

    // If we have gotten to the "lots of fusions" part of the solution space, start
    // taking smaller step sizes.
    if(nzeros_new > 0){
      t = t_switch;
    }

    // If we have seen a fusion or are otherwise interested in keeping this iteration,
    // add values to our storage buffers
    if( (nzeros_new != nzeros_old) | (iter % keep == 0) ){
      // Before we can store values, we need to make sure we have enough buffer space
      if(path_iter >= buffer_size){
        buffer_size *= 2; // Double our buffer sizes
        UPath.conservativeResize(UPath.rows(), buffer_size);
        VPath.conservativeResize(VPath.rows(), buffer_size);
        gamma_path.conservativeResize(buffer_size);
        vZeroInds_Path.conservativeResize(vZeroInds_Path.rows(), buffer_size);
      }

      // Store values
      UPath.col(path_iter)          = u_new;
      VPath.col(path_iter)          = v_new;
      gamma_path(path_iter)         = gamma;
      vZeroInds_Path.col(path_iter) = vZeroIndsnew;

      path_iter++;

      // Update gamma_old (basis for future back-tracking searches)
      gamma_old = gamma;
    }

    iter++;
    gamma *= t;

    if((iter % CLUSTRVIZ_CHECK_USER_INTERRUPT_RATE) == 0){
      Rcpp::checkUserInterrupt();
    }
  }

  // Now that we are done, we can "drop" unused buffer space before returning to R
  UPath.conservativeResize(UPath.rows(), path_iter - 1);
  VPath.conservativeResize(VPath.rows(), path_iter - 1);
  gamma_path.conservativeResize(path_iter - 1);
  vZeroInds_Path.conservativeResize(vZeroInds_Path.rows(), path_iter - 1);

  // Wrap up our results and pass them to R
  return Rcpp::List::create(Rcpp::Named("u.path")      = UPath,
                            Rcpp::Named("v.path")      = VPath,
                            Rcpp::Named("v.zero.inds") = vZeroInds_Path,
                            Rcpp::Named("lambda.path")  = gamma_path); // TODO - Change lambda -> gamma in R code
}
