#include "clustRviz.h"

// TODO - Consolidate CBASS and CBASS-VIZ
// Most of the internal logic is the same (modulo back-tracking vs fixed step size)

// [[Rcpp::export]]
Rcpp::List CBASS_VIZcpp(const Eigen::VectorXd& x,
                        int n,
                        int p,
                        double lambda_init, // TODO: Change to gamma_init
                        const Eigen::VectorXd& weights_col,
                        const Eigen::VectorXd& weights_row,
                        const Eigen::VectorXd& uinit_row, // TODO: Change to u_init_row
                        const Eigen::VectorXd& uinit_col, // TODO: Change to u_init_col
                        const Eigen::VectorXd& vinit_row, // TODO: Change to v_init_row
                        const Eigen::VectorXd& vinit_col, // TODO: Change to v_init_col
                        const Eigen::SparseMatrix<double>& premat_row,
                        const Eigen::SparseMatrix<double>& premat_col,
                        const Eigen::MatrixXi& IndMat_row,
                        const Eigen::MatrixXi& IndMat_col,
                        const Eigen::MatrixXi& EOneIndMat_row,
                        const Eigen::MatrixXi& EOneIndMat_col,
                        const Eigen::MatrixXi& ETwoIndMat_row,
                        const Eigen::MatrixXi& ETwoIndMat_col,
                        double rho         = 1,
                        int max_iter       = 10000,
                        int burn_in        = 50,
                        bool verbose       = false,
                        int ti             = 15,
                        double t_switch    = 1.01,
                        int keep           = 10,
                        bool l1            = false){

  // Typically, our weights are "sparse" (i.e., mostly zeros) because we
  // drop small weights to achieve performance.
  Eigen::Index num_edges_row = EOneIndMat_row.rows();
  Eigen::Index num_edges_col = EOneIndMat_col.rows();

  /// Set-up storage for CBASS iterates

  // In order to pre-allocate storage arrays, we need to estimate the number of
  // steps with fusions we will encounter. Dendrograms are the common case and
  // we are bi-clustering, so we expect O(n + p) fusions. It's a bit
  // cheaper to drop observations than to extend the internal buffers of our
  // storage objects, so we use 1.5(n + p) for now
  Eigen::Index buffer_size = 1.5 * (n + p);

  // Primal variable
  Eigen::VectorXd u_new(n * p);              // Working copies
  Eigen::VectorXd u_old(n * p);
  u_new = x;
  Eigen::MatrixXd UPath(n * p, buffer_size); // Storage (for values to return to R)
  UPath.col(0) = uinit_col;

  // 'Split' variable for row fusions
  Eigen::VectorXd v_new_row(p * num_edges_row);              // Working copies
  Eigen::VectorXd v_old_row(p * num_edges_row);
  v_new_row = vinit_row;
  Eigen::MatrixXd VPath_row(p * num_edges_row, buffer_size); // Storage (for values to return to R)
  VPath_row.col(0) = vinit_row;

  // 'Split' variable for column fusions
  Eigen::VectorXd v_new_col(p * num_edges_col);              // Working copies
  Eigen::VectorXd v_old_col(p * num_edges_col);
  v_new_col = vinit_col;
  Eigen::MatrixXd VPath_col(p * num_edges_col, buffer_size); // Storage (for values to return to R)
  VPath_col.col(0) = vinit_col;

  // (Scaled) dual variable for row fusions
  Eigen::VectorXd z_new_row(p * num_edges_row);
  Eigen::VectorXd z_old_row(p * num_edges_row);
  z_new_row = v_new_row;

  // (Scaled) dual variable for column fusions
  Eigen::VectorXd z_new_col(p * num_edges_col);
  Eigen::VectorXd z_old_col(p * num_edges_col);
  z_new_col = v_new_col;

  // Regularization level
  double gamma = lambda_init;              // Working copies
  double gamma_old = lambda_init;
  Eigen::VectorXd gamma_path(buffer_size); // Storage (to be returned to R)
  gamma_path(0) = lambda_init;

  // Row Fusions -- TODO: Confirm the semantics of these objects with JN
  Eigen::VectorXd vZeroIndsnew_row = Eigen::VectorXd::Zero(num_edges_row);         // Working copies
  Eigen::VectorXd vZeroIndsold_row(num_edges_row);                                 // (we begin with no fusions)
  Eigen::MatrixXd vZeroIndsPath_row = Eigen::MatrixXd(num_edges_row, buffer_size); // Storage (to be returned to R)
  vZeroIndsPath_row.col(0) = vZeroIndsnew_row;

  // Column Fusions -- TODO: Confirm the semantics of these objects with JN
  Eigen::VectorXd vZeroIndsnew_col = Eigen::VectorXd::Zero(num_edges_col);         // Working copies
  Eigen::VectorXd vZeroIndsold_col(num_edges_col);                                 // (we begin with no fusions)
  Eigen::MatrixXd vZeroIndsPath_col = Eigen::MatrixXd(num_edges_col, buffer_size); // Storage (to be returned to R)
  vZeroIndsPath_col.col(0) = vZeroIndsnew_col;

  // The COBRA algorithm for Convex Bi-Clustering (on which CBASS is based)
  // introduces extra variables Y, P, Q which are necessary to keep the row-
  // and column-fusions appropriately coupled.
  //
  // See Chi, Allen, Baraniuk (2017) for details
  //
  // Since we are working in "vec" form, instead of matrix form, we use the
  // `restride()` function to re-organize a vector from an ordering appropriate
  // for the row-based steps to an ordering for the column-based steps and vice versa
  //
  // Note that we use a different notational convention here (capital P vs lower-case p)
  // to distinguish this "p" from the "p" that gives the problem dimension
  //
  // Y is not carried forward from one iteration to the next, so we declare it in-loop below
  Eigen::VectorXd P_new = Eigen::VectorXd::Zero(n * p);
  Eigen::VectorXd P_old(n * p);
  Eigen::VectorXd Q_new = Eigen::VectorXd::Zero(n * p);
  Eigen::VectorXd Q_old(n * p);

  /// END Preallocations

  // At each iteration, we need to calculate A^{-1}B_k for some (sparse) A [one for rows and one for cols]
  // This is a relatively expensive iteration, but the core cost is a sparse LU
  // factorization of A which can be amortized over iterations so we pre-compute them here
  Eigen::SparseLU<Eigen::SparseMatrix<double> > premat_solver_row;
  premat_solver_row.compute(premat_row);

  Eigen::SparseLU<Eigen::SparseMatrix<double> > premat_solver_col;
  premat_solver_col.compute(premat_col);

  // Book-keeping variables
  // Number of iterations stored, total iteration count, number of column fusions, number of row fusions
  Eigen::Index path_iter  = 1; // path_iter is next column to put something in,
  Eigen::Index iter       = 0; // so we start at 1 since we put data in column 0 above
  Eigen::Index nzeros_old_row = 0;
  Eigen::Index nzeros_new_row = 0;
  Eigen::Index nzeros_old_col = 0;
  Eigen::Index nzeros_new_col = 0;

  // We begin CBASS-VIZ by taking relatively large step sizes (t = 1.1)
  // but once we get to an "interesting" part of the path, we switch to
  // smaller step sizes (as determined by t_switch)
  double t = 1.1;

  while( ((nzeros_new_row < num_edges_row) | (nzeros_new_col < num_edges_col)) & (iter < max_iter)){
    // Begin iteration - move updated values to "_old" values
    //
    // TODO -- Do this as a swap and avoid full copies if possible
    u_old = u_new;                       // Primal
    v_old_row = v_new_row;               // Split
    v_old_col = v_new_col;
    z_old_row = z_new_row;               // Dual
    z_old_col = z_new_col;
    nzeros_old_row = nzeros_new_row;     // Fusion counts
    nzeros_old_col = nzeros_new_col;
    vZeroIndsold_row = vZeroIndsnew_row; // Fusion indices
    vZeroIndsold_col = vZeroIndsnew_col;

    P_old = P_new;
    Q_old = Q_new;
    Eigen::VectorXd P_old_t = restride(P_old, p);
    Eigen::VectorXd u_old_t = restride(u_old, p);

    // VIZ book-keeping
    bool rep_iter = true;
    Eigen::Index try_iter = 0;

    double gamma_upper = gamma;
    double gamma_lower = gamma_old;

    // This is the core CBASS-VIZ Logic:
    //
    // We take CBASS-type (one iteration of each ADMM step) steps, but instead of
    // proceeding with a fixed step-size update, we include a back-tracking step
    // (described in more detail below)
    //
    //
    while(rep_iter){
      /// Row-fusion iterations
      // U-update
      Eigen::VectorXd solver_input_row = DtMatOpv2(rho * v_old_row - z_old_row, p, n, IndMat_row, EOneIndMat_row, ETwoIndMat_row);
      solver_input_row += P_old_t + u_old_t;
      solver_input_row /= rho;
      Eigen::VectorXd Y_t = premat_solver_row.solve(solver_input_row);

      // V-update
      Eigen::VectorXd prox_argument_row = DMatOpv2(Y_t,n, IndMat_row, EOneIndMat_row, ETwoIndMat_row) + (1/rho)*z_old_row;
      if(l1){
        v_new_row = ProxL1(prox_argument_row, n, (1/rho) * gamma, weights_row);
      } else {
        v_new_row = ProxL2(prox_argument_row, n, (1/rho) * weights_row * gamma, IndMat_row);
      }

      // Z-update
      z_new_row = z_old_row + rho*(DMatOpv2(Y_t, n, IndMat_row, EOneIndMat_row, ETwoIndMat_row) - v_new_row);
      /// END Row-fusion iterations

      Eigen::VectorXd Y = restride(Y_t, n);
      P_new = u_old + P_old - Y;

      /// Column-fusion iterations
      // U-update
      Eigen::VectorXd solver_input_col = DtMatOpv2(rho * v_old_col - z_old_col, n, p, IndMat_col, EOneIndMat_col, ETwoIndMat_col);
      solver_input_col += Y + Q_old;
      solver_input_col /= rho;
      u_new = premat_solver_col.solve(solver_input_col);

      // V-update
      Eigen::VectorXd prox_argument_col = DMatOpv2(u_new, p, IndMat_col, EOneIndMat_col, ETwoIndMat_col) + (1/rho)*z_old_col;
      if(l1){
        v_new_col = ProxL1(prox_argument_col, p, (1/rho) * gamma, weights_col);
      } else {
        v_new_col = ProxL2(prox_argument_col, p, (1/rho) * weights_col * gamma, IndMat_col);
      }

      // Z-update
      z_new_col = z_old_col + rho*(DMatOpv2(u_new, p, IndMat_col, EOneIndMat_col, ETwoIndMat_col)-v_new_col);
      /// END Column-fusion iterations

      Q_new = Y + Q_old - u_new;

      // Count number of row fusions
      for(int l = 0; l < num_edges_row; l++){
        Eigen::VectorXi v_index_row = IndMat_row.row(l);
        if(extract(v_new_row, v_index_row).sum() == 0){
          vZeroIndsnew_row(l) = 1;
        }
      }
      nzeros_new_row = vZeroIndsnew_row.sum();

      // Count number of column fusions
      for(int l = 0; l < num_edges_col; l++){
        Eigen::VectorXi v_index_col = IndMat_col.row(l);
        if(extract(v_new_col, v_index_col).sum() == 0){
          vZeroIndsnew_col(l) = 1;
        }
      }
      nzeros_new_col = vZeroIndsnew_col.sum();

      /// END CBASS steps

      try_iter++; // Increment internal iteration count (used to check stopping below)

      if( (nzeros_new_row == nzeros_old_row) & (nzeros_new_col == nzeros_old_col) & (try_iter == 1)){
        // If the sparsity pattern (number of fusions) hasn't changed, we have
        // no need to back-track (we didn't miss any fusions) so we can go immediately
        // to the next iteration.
        rep_iter = false;
      } else if((nzeros_new_row > nzeros_old_row + 1) | (nzeros_new_col > nzeros_old_col + 1) ){
        // If we see two (or more) new fusions, we need to back-track and figure
        // out which one occured first
        // (NB: we don't need a 'cross' ordering of row and global fusions, so we
        //  use a separate check for each, instead of a check on the sum)
        vZeroIndsnew_col = vZeroIndsold_col;
        vZeroIndsnew_row = vZeroIndsold_row;
        if(try_iter == 1){
          gamma = 0.5 * (gamma_lower + gamma_upper);
        } else{
          gamma_upper = gamma;
          gamma = 0.5 * (gamma_lower + gamma_upper);
        }
      } else if( (nzeros_new_row == nzeros_old_row) & (nzeros_new_col == nzeros_old_col) ){
        // If we don't observe any new fusions, we take another iteration without
        vZeroIndsnew_col = vZeroIndsold_col;
        vZeroIndsnew_row = vZeroIndsold_row;
        gamma_lower = gamma;
        gamma = 0.5 * (gamma_lower + gamma_upper);
      } else{
        // If we see exactly one new fusion, we have a good step size and exit
        // the inner back-tracking loop
        rep_iter = false;
      }

      // Safety check - only so many iterations of the inner loop before we move on
      if(try_iter > ti){
        rep_iter = false;
      }
    }

    // If we have gotten to the "lots of fusions" part of the solution space, start
    // taking smaller step sizes.
    if( (nzeros_new_col > 0) | (nzeros_new_row >0) ){
      t = t_switch;
    }

    // If we have seen a fusion or are otherwise interested in keeping this iteration,
    // add values to our storage buffers
    if( (nzeros_new_row != nzeros_old_row) | (nzeros_new_col != nzeros_old_col) | (iter % keep == 0) ){

      // Before we can store values, we need to make sure we have enough buffer space
      if(path_iter >= buffer_size){
        buffer_size *= 2; // Double our buffer sizes
        UPath.conservativeResize(UPath.rows(), buffer_size);
        VPath_row.conservativeResize(VPath_row.rows(), buffer_size);
        VPath_col.conservativeResize(VPath_col.rows(), buffer_size);
        gamma_path.conservativeResize(buffer_size);
        vZeroIndsPath_row.conservativeResize(vZeroIndsPath_row.rows(), buffer_size);
        vZeroIndsPath_col.conservativeResize(vZeroIndsPath_col.rows(), buffer_size);
      }

      // Store values
      UPath.col(path_iter)             = u_new;
      VPath_row.col(path_iter)         = v_new_row;
      VPath_col.col(path_iter)         = v_new_col;
      gamma_path(path_iter)            = gamma;
      vZeroIndsPath_row.col(path_iter) = vZeroIndsnew_row;
      vZeroIndsPath_col.col(path_iter) = vZeroIndsnew_col;

      path_iter++;

      // Update gamma old as the basis for future back-tracking
      gamma_old = gamma;
    }

    iter++;
    if(iter >= burn_in){
      gamma *= t;
    }

    if((iter % CLUSTRVIZ_CHECK_USER_INTERRUPT_RATE) == 0){
      Rcpp::checkUserInterrupt();
    }
  }

  // Now that we are done, we can "drop" unused buffer space before returning to R
  //
  // See explanatory comment in carp_viz.cpp
  UPath.conservativeResize(UPath.rows(), path_iter);
  VPath_row.conservativeResize(VPath_row.rows(), path_iter);
  VPath_col.conservativeResize(VPath_col.rows(), path_iter);
  gamma_path.conservativeResize(path_iter);
  vZeroIndsPath_row.conservativeResize(vZeroIndsPath_row.rows(), path_iter);
  vZeroIndsPath_col.conservativeResize(vZeroIndsPath_col.rows(), path_iter);

  // Wrap up our results and pass them to R
  return Rcpp::List::create(Rcpp::Named("u.path")          = UPath,
                            Rcpp::Named("v.row.path")      = VPath_row,
                            Rcpp::Named("v.col.path")      = VPath_col,
                            Rcpp::Named("v.row.zero.inds") = vZeroIndsPath_row,
                            Rcpp::Named("v.col.zero.inds") = vZeroIndsPath_col,
                            Rcpp::Named("lambda.path")      = gamma_path); // TODO - Change lambda -> gamma in R code
}
