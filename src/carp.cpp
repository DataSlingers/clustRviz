#include "clustRviz.h"

// [[Rcpp::export]]
Rcpp::List CARPcpp(const Eigen::MatrixXd& X,
                   const Eigen::MatrixXd& D,
                   double lambda_init, // TODO: Change to gamma_init
                   double t,
                   const Eigen::VectorXd& weights,
                   double rho   = 1,
                   int max_iter = 10000,
                   int burn_in  = 50,
                   int keep     = 10,
                   bool l1      = false){

  Eigen::Index n = X.rows();
  Eigen::Index p = X.cols();
  // Typically, our weights are "sparse" (i.e., mostly zeros) because we
  // drop small weights to achieve performance.
  Eigen::Index num_edges = D.rows();

  /// Set-up storage for CARP iterates

  // In order to pre-allocate storage arrays, we need to estimate the number of
  // steps with fusions we will encounter. Dendrograms are the common case and
  // we are only clustering observations, so we expect O(n) fusions. It's a bit
  // cheaper to drop observations than to extend the internal buffers of our
  // storage objects, so we use 1.5n for now
  Eigen::Index buffer_size = 1.5 * n;

  // Primal variable
  Eigen::MatrixXd U = X;
  Eigen::MatrixXd UPath(n * p, buffer_size); // Storage (for values to return to R)
  // We cannot directly copy matrices into columns of our return object since
  // a matrix isn't a vector, so we map the same storage into a vector which we
  // can then insert into the storage matrix to be returned to R
  Eigen::Map<Eigen::VectorXd> u_vec = Eigen::Map<Eigen::VectorXd>(U.data(), n * p);

  // 'Split' variable
  Eigen::MatrixXd V = D * U;
  Eigen::MatrixXd VPath(p * num_edges, buffer_size); // Storage (for values to return to R)
  Eigen::Map<Eigen::VectorXd> v_vec = Eigen::Map<Eigen::VectorXd>(V.data(), p * num_edges);

  U.transposeInPlace(); // See note on transpositions below
  UPath.col(0) = u_vec;
  U.transposeInPlace();
  V.transposeInPlace();
  VPath.col(0) = v_vec;
  V.transposeInPlace();

  // (Scaled) dual variable
  Eigen::MatrixXd Z = V;

  // Regularization level
  double gamma = lambda_init;              // Working copy
  Eigen::VectorXd gamma_path(buffer_size); // Storage (to be returned to R)
  gamma_path(0) = lambda_init;

  // Fusions
  Eigen::MatrixXi v_zeros_path(num_edges, buffer_size); // Storage (to be returned to R)
  v_zeros_path.col(0).setZero();

  /// END Preallocations

  // PreCompute chol(I + rho D^TD) for easy inversions in the U update step
  Eigen::MatrixXd IDTD = rho * D.transpose() * D + Eigen::MatrixXd::Identity(n, n);
  Eigen::LLT<Eigen::MatrixXd> u_step_solver; u_step_solver.compute(IDTD);

  // Book-keeping variables
  // Number of iterations stored, total iteration count, number of fusions
  Eigen::Index path_iter  = 1; // path_iter is next column to put something in,
  Eigen::Index iter       = 0; // so we start at 1 since we put data in column 0 above
  Eigen::Index nzeros_old = 0;
  Eigen::Index nzeros_new = 0;

  while( (iter < max_iter) & (nzeros_new < num_edges) ){
    ClustRVizLogger::info("Beginning iteration k = ") << iter + 1;
    ClustRVizLogger::debug("gamma = ") << gamma;

    nzeros_old = nzeros_new;

    // U-update
    U = u_step_solver.solve(X + rho * D.transpose() * (V - Z));
    Eigen::MatrixXd DU = D * U;
    ClustRVizLogger::debug("U = ") << U;

    // V-update
    Eigen::MatrixXd DUZ = DU + Z;

    if(l1){
      for(Eigen::Index i = 0; i < num_edges; i++){
        for(Eigen::Index j = 0; j < p; j++){
          V(i, j) = soft_thresh(DUZ(i, j), gamma / rho * weights(i));
        }
      }
    } else {
      for(Eigen::Index i = 0; i < num_edges; i++){
        Eigen::VectorXd DUZ_i = DUZ.row(i);
        double scale_factor = 1 - gamma * weights(i) / (rho * DUZ_i.norm());

        if(scale_factor > 0){
          V.row(i) = DUZ_i * scale_factor;
        } else {
          V.row(i).setZero();
        }
      }
    }

    ClustRVizLogger::debug("V = ") << V;

    // Z-update
    Z += DU - V;

    ClustRVizLogger::debug("Z = ") << Z;

    // Identify cluster fusions (rows of V which have gone to zero)
    Eigen::VectorXd v_norms = V.rowwise().squaredNorm();
    Eigen::ArrayXi  v_zeros(num_edges);

    for(Eigen::Index i = 0; i < num_edges; i++){
      v_zeros(i) = v_norms(i) == 0;
    }
    nzeros_new = v_zeros.sum();

    ClustRVizLogger::debug("Number of fusions identified ") << nzeros_new;

    // If we have seen a fusion or are otherwise interested in keeping this iteration,
    // add values to our storage buffers
    if( (nzeros_new != nzeros_old) | (iter % keep == 0) ) {
      // Before we can store values, we need to make sure we have enough buffer space
      if(path_iter >= buffer_size){
        ClustRVizLogger::info("Resizing storage from ") << buffer_size << " to " << 2 * buffer_size << " iterations.";
        buffer_size *= 2; // Double our buffer sizes
        UPath.conservativeResize(UPath.rows(), buffer_size);
        VPath.conservativeResize(VPath.rows(), buffer_size);
        gamma_path.conservativeResize(buffer_size);
        v_zeros_path.conservativeResize(v_zeros_path.rows(), buffer_size);
      }

      // Store values

      // FIXME -- The post-processing code assumes output in the form of
      //          Chi and Lange (JCGS, 2015) which more or less is equivalent
      //          to vec(U^T) and vec(V^T) instead of what we're using internally
      //          here.
      //
      //          It should be re-written, but until it is, we can achieve the
      //          same result by transposing and un-transposing the data internally
      //          before copying it to our storage buffers
      //
      //          Obviously, this burns some cycles so it would be better to avoid this
      U.transposeInPlace(); V.transposeInPlace();

      UPath.col(path_iter)          = u_vec;
      VPath.col(path_iter)          = v_vec;
      gamma_path(path_iter)         = gamma;
      v_zeros_path.col(path_iter)   = v_zeros;

      U.transposeInPlace(); V.transposeInPlace();

      path_iter++;
    }

    iter++;
    if(iter >= burn_in){
      gamma *= t;
    }

    if((iter % CLUSTRVIZ_CHECK_USER_INTERRUPT_RATE) == 0){
      Rcpp::checkUserInterrupt();
    }
  }

  if(iter >= max_iter){
    ClustRVizLogger::warning("CARP ended early -- `max_iter` reached");
  }

  // Now that we are done, we can "drop" unused buffer space before returning to R
  //
  // See explanatory comment in carp_viz.cpp
  UPath.conservativeResize(UPath.rows(), path_iter);
  VPath.conservativeResize(VPath.rows(), path_iter);
  gamma_path.conservativeResize(path_iter);
  v_zeros_path.conservativeResize(v_zeros_path.rows(), path_iter);

  // Wrap up our results and pass them to R
  return Rcpp::List::create(Rcpp::Named("u.path")      = UPath,
                            Rcpp::Named("v.path")      = VPath,
                            Rcpp::Named("v.zero.inds") = v_zeros_path,
                            Rcpp::Named("lambda.path") = gamma_path); // TODO - Change lambda -> gamma in R code
}
