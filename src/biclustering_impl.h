#ifndef CLUSTRVIZ_BICLUSTERING_H
#define CLUSTRVIZ_BICLUSTERING_H 1

#include "clustRviz_base.h"
#include "clustRviz_logging.h"
#include "status.h"

class ConvexBiClustering {
public:
  double gamma; // Current regularization level - need to be able to manipulate this externally

  ConvexBiClustering(const Eigen::MatrixXd& X_,
                     const Eigen::MatrixXd& D_row_,
                     const Eigen::MatrixXd& D_col_,
                     const Eigen::VectorXd& weights_row_,
                     const Eigen::VectorXd& weights_col_,
                     const double rho_,
                     const bool l1_,
                     const bool show_progress_):
    X(X_),
    D_row(D_row_),
    D_col(D_col_),
    weights_row(weights_row_),
    weights_col(weights_col_),
    rho(rho_),
    l1(l1_),
    n(X_.rows()),
    p(X_.cols()),
    num_row_edges(D_row_.rows()),
    num_col_edges(D_col_.cols()),
    sp(show_progress_, D_row_.rows() + D_col_.cols()){

    // Set initial values for optimization variables
    U = X;
    V_row = D_row * U;
    Z_row = V_row;
    V_col = D_col.transpose() * X.transpose();
    Z_col = V_col;

    P = Eigen::MatrixXd::Zero(n, p);
    Q = Eigen::MatrixXd::Zero(n, p);

    v_row_zeros = Eigen::ArrayXi::Zero(num_row_edges);
    v_col_zeros = Eigen::ArrayXi::Zero(num_col_edges);
    gamma = 0;

    sp.set_v_norm_init(V_row.squaredNorm() + V_col.squaredNorm());

    // Initialize storage buffers
    buffer_size = 1.5 * (n + p);
    UPath.resize(n * p, buffer_size);
    V_rowPath.resize(p * num_row_edges, buffer_size);
    V_colPath.resize(n * num_col_edges, buffer_size);
    gamma_path.resize(buffer_size);
    v_row_zeros_path.resize(num_row_edges, buffer_size);
    v_col_zeros_path.resize(num_col_edges, buffer_size);

    // Store initial values
    nzeros_row = 0;
    nzeros_col = 0;
    storage_index = 0;
    store_values();

    // PreCompute chol(I_n + rho D_row^T D_row) and chol(I_p + rho D_col D_col^T) for easy inversions in the ADMM primal update steps
    Eigen::MatrixXd IDTD_row = rho * D_row.transpose() * D_row + Eigen::MatrixXd::Identity(n, n);
    row_primal_solver.compute(IDTD_row);

    Eigen::MatrixXd IDDT_col = rho * D_col * D_col.transpose() + Eigen::MatrixXd::Identity(p, p);
    col_primal_solver.compute(IDDT_col);
  };

  bool is_interesting_iter(){
    // FIXME? A better check would be fusion IDs, not just number
    return (nzeros_row != nzeros_row_old) | (nzeros_col != nzeros_col_old);
  }

  bool multiple_fusions(){
    return (nzeros_row > nzeros_row_old + 1) | (nzeros_col > nzeros_col_old + 1);
  }

  bool is_complete(){
    // Stop when all edges are active (when we have a connected graph, all elements are fused)
    return (nzeros_row == num_row_edges) & (nzeros_col == num_col_edges);
  }

  void full_admm_step(){
    Eigen::MatrixXd T = U + P;
    Eigen::MatrixXd T_prev;
    int k_row = 0;

    do {
      k_row ++;
      T_prev = T;
      /// Row-fusion iterations
      // Primal Update
      T = row_primal_solver.solve(U + P + rho * D_row.transpose() * (V_row - Z_row));

      Eigen::MatrixXd DT = D_row * T;
      ClustRVizLogger::debug("T = ") << T;

      // Copy Update
      Eigen::MatrixXd DTZ = DT + Z_row;
      V_row = MatrixProx(DTZ, gamma / rho, weights_row, l1);
      ClustRVizLogger::debug("V_row = ") << V_row;

      // Dual Update
      Z_row += DT - V_row;
      ClustRVizLogger::debug("Z_row = ") << Z_row;
      /// END Row-fusion iterations

      if(k_row > 150) {
        break; // Avoid infinite looping...
      }
    } while ( (T - T_prev).squaredNorm() > 1e-7);

    // DLPA Updates
    P += U - T;
    Eigen::MatrixXd TQT = T + Q; TQT.transposeInPlace();

    Eigen::MatrixXd S = TQT;
    Eigen::MatrixXd S_prev;
    int k_col = 0;

    do {
      S_prev = S;
      k_col++;

      /// Column-fusion iterations
      // Primal Update
      S = col_primal_solver.solve(TQT + rho * D_col * (V_col - Z_col));
      Eigen::MatrixXd DTS = D_col.transpose() * S;
      ClustRVizLogger::debug("S = ") << S;

      // Copy Update
      Eigen::MatrixXd DTSZ = DTS + Z_col;
      V_col = MatrixProx(DTSZ, gamma / rho, weights_col, l1);
      ClustRVizLogger::debug("V_col = ") << V_col;

      // Dual Update
      Z_col += DTS - V_col;
      ClustRVizLogger::debug("Z_col = ") << Z_col;
      /// END Column-fusion iterations

      if(k_col > 150) {
        break; // Avoid infinite looping...
      }

    } while ( (S - S_prev).squaredNorm() > 1e-7);

    // DLPA Updates + New U
    U = S.transpose();
    Q += T - U;

    // Identify row fusions (rows of V_row which have gone to zero)
    Eigen::VectorXd v_row_norms = V_row.rowwise().squaredNorm();

    for(Eigen::Index i = 0; i < num_row_edges; i++){
      v_row_zeros(i) = v_row_norms(i) == 0;
    }

    nzeros_row = v_row_zeros.sum();

    // Identify column fusions (rows of V_col which have gone to zero)
    // Remember, V_col and Z_col are internal to the "transposed prox" sub-problem
    // so everything is reversed of what we'd expect
    Eigen::VectorXd v_col_norms = V_col.rowwise().squaredNorm();

    for(Eigen::Index i = 0; i < num_col_edges; i++){
      v_col_zeros(i) = v_col_norms(i) == 0;
    }

    nzeros_col = v_col_zeros.sum();

    ClustRVizLogger::debug("Number of row fusions identified ") << nzeros_row;
    ClustRVizLogger::debug("Number of column fusions identified ") << nzeros_col;
  }

  void admm_step(){
    /// Row-fusion iterations
    // Primal Update
    Eigen::MatrixXd T = row_primal_solver.solve(U + P + rho * D_row.transpose() * (V_row - Z_row));

    Eigen::MatrixXd DT = D_row * T;
    ClustRVizLogger::debug("T = ") << T;

    // Copy Update
    Eigen::MatrixXd DTZ = DT + Z_row;
    V_row = MatrixProx(DTZ, gamma / rho, weights_row, l1);
    ClustRVizLogger::debug("V_row = ") << V_row;

    // Dual Update
    Z_row += DT - V_row;
    ClustRVizLogger::debug("Z_row = ") << Z_row;
    /// END Row-fusion iterations

    // DLPA Updates
    P += U - T;
    Eigen::MatrixXd TQT = T + Q; TQT.transposeInPlace();

    /// Column-fusion iterations
    // Primal Update
    Eigen::MatrixXd S = col_primal_solver.solve(TQT + rho * D_col * (V_col - Z_col));
    Eigen::MatrixXd DTS = D_col.transpose() * S;
    ClustRVizLogger::debug("S = ") << S;

    // Copy Update
    Eigen::MatrixXd DTSZ = DTS + Z_col;
    V_col = MatrixProx(DTSZ, gamma / rho, weights_col, l1);
    ClustRVizLogger::debug("V_col = ") << V_col;

    // Dual Update
    Z_col += DTS - V_col;
    ClustRVizLogger::debug("Z_col = ") << Z_col;
    /// END Column-fusion iterations

    // DLPA Updates + New U
    U = S.transpose();
    Q += T - U;

    // Identify row fusions (rows of V_row which have gone to zero)
    Eigen::VectorXd v_row_norms = V_row.rowwise().squaredNorm();

    for(Eigen::Index i = 0; i < num_row_edges; i++){
      v_row_zeros(i) = v_row_norms(i) == 0;
    }

    nzeros_row = v_row_zeros.sum();

    // Identify column fusions (rows of V_col which have gone to zero)
    // Remember, V_col and Z_col are internal to the "transposed prox" sub-problem
    // so everything is reversed of what we'd expect
    Eigen::VectorXd v_col_norms = V_col.rowwise().squaredNorm();

    for(Eigen::Index i = 0; i < num_col_edges; i++){
      v_col_zeros(i) = v_col_norms(i) == 0;
    }

    nzeros_col = v_col_zeros.sum();

    ClustRVizLogger::debug("Number of row fusions identified ") << nzeros_row;
    ClustRVizLogger::debug("Number of column fusions identified ") << nzeros_col;
  }

  void save_fusions(){
    nzeros_row_old = nzeros_row;
    nzeros_col_old = nzeros_col;
  }

  void save_old_values(){
    U_old = U;
    P_old = P;
    Q_old = Q;

    V_row_old = V_row;
    Z_row_old = Z_row;
    v_row_zeros_old = v_row_zeros;

    V_col_old = V_col;
    Z_col_old = Z_col;
    v_col_zeros_old = v_col_zeros;
  }

  void load_old_variables(){
    U = U_old;
    P = P_old;
    Q = Q_old;

    V_row = V_row_old;
    Z_row = Z_row_old;
    V_col = V_col_old;
    Z_col = Z_col_old;
  }

  void load_old_fusions(){
    v_row_zeros = v_row_zeros_old;
    v_col_zeros = v_col_zeros_old;
  }

  bool has_fusions(){
    return (nzeros_row > 0) | (nzeros_col > 0);
  }

  bool admm_converged(){
    return (U - U_old).squaredNorm() < 1e-10;
  }

  void reset_aux(){
    U = U_old = X;
    P.setZero();
    P_old.setZero();
    Q.setZero();
    Q_old.setZero();
  }

  void store_values(){
    if(storage_index >= buffer_size){
      ClustRVizLogger::info("Resizing storage from ") << buffer_size << " to " << 2 * buffer_size << " iterations.";
      buffer_size *= 2; // Double our buffer sizes
      UPath.conservativeResize(UPath.rows(), buffer_size);
      V_rowPath.conservativeResize(V_rowPath.rows(), buffer_size);
      V_colPath.conservativeResize(V_colPath.rows(), buffer_size);
      gamma_path.conservativeResize(buffer_size);
      v_row_zeros_path.conservativeResize(v_row_zeros_path.rows(), buffer_size);
      v_col_zeros_path.conservativeResize(v_col_zeros_path.rows(), buffer_size);
    }

    // Store values
    UPath.col(storage_index)        = Eigen::Map<Eigen::VectorXd>(U.data(), n * p);
    V_rowPath.col(storage_index)    = Eigen::Map<Eigen::VectorXd>(V_row.data(), p * num_row_edges);
    V_colPath.col(storage_index)    = Eigen::Map<Eigen::VectorXd>(V_col.data(), n * num_col_edges);
    gamma_path(storage_index)       = gamma;
    v_row_zeros_path.col(storage_index) = v_row_zeros;
    v_col_zeros_path.col(storage_index) = v_col_zeros;

    storage_index++;
  }

  Rcpp::List build_return_object(){
    // When we are done, we can "drop" unused buffer space before returning to R
    //
    // storage_index is the zero-based index of the next column we would use for storage,
    // but it is also the (one-based) _number_ of columns we want to save so no need
    // to adjust. (NB: conservativeResize takes the target size, not the columns to keep
    // as an argument)
    UPath.conservativeResize(UPath.rows(), storage_index);
    V_rowPath.conservativeResize(V_rowPath.rows(), storage_index);
    V_colPath.conservativeResize(V_colPath.rows(), storage_index);
    gamma_path.conservativeResize(storage_index);
    v_row_zeros_path.conservativeResize(v_row_zeros_path.rows(), storage_index);
    v_col_zeros_path.conservativeResize(v_col_zeros_path.rows(), storage_index);

    return Rcpp::List::create(Rcpp::Named("u_path")          = UPath,
                              Rcpp::Named("v_row_path")      = V_rowPath,
                              Rcpp::Named("v_col_path")      = V_colPath,
                              Rcpp::Named("v_row_zero_inds") = v_row_zeros_path,
                              Rcpp::Named("v_col_zero_inds") = v_col_zeros_path,
                              Rcpp::Named("gamma_path")      = gamma_path);
  }

  void tick(uint iter){
    sp.update(nzeros_row + nzeros_col,
              V_row.squaredNorm() + V_col.squaredNorm(),
              iter,
              gamma);
  }

private:
  // Fixed (non-data-dependent) problem details
  const Eigen::MatrixXd& X; // Data matrix (to be clustered)
  const Eigen::MatrixXd& D_row; // Edge (differencing) matrix
  const Eigen::MatrixXd& D_col;
  const Eigen::VectorXd& weights_row; // Clustering weights
  const Eigen::VectorXd& weights_col;
  const double rho; // ADMM relaxation parameter -- TODO: Factor this out?
  // Theoretically, it's part of the algorithm, not the problem
  // but we need it in the steps...
  bool  l1;         // Is the L1 (true) or L2 (false) norm being used?
  const int n;      // Problem dimensions
  const int p;
  const int num_row_edges;
  const int num_col_edges;
  Eigen::LLT<Eigen::MatrixXd> row_primal_solver; // Cached factorizations for primal updates
  Eigen::LLT<Eigen::MatrixXd> col_primal_solver;

  // Progress printer
  StatusPrinter sp;

  // Current copies of ADMM variables
  Eigen::MatrixXd U;     // Primal Variable
  Eigen::MatrixXd V_row; // Split Variable - row subproblem
  Eigen::MatrixXd Z_row; // Dual Variable - row subproblem
  Eigen::MatrixXd V_col; // Split Variable - column subproblem
  Eigen::MatrixXd Z_col; // Dual Variable - column subproblem
  Eigen::ArrayXi v_row_zeros; // Fusion indicators
  Eigen::ArrayXi v_col_zeros;
  Eigen::Index nzeros_row; // Fusion counts
  Eigen::Index nzeros_col;

  // The DLPA (on which CBASS is based) adds two auxiliary variables -- P & Q --
  // with the same dimensions as the optimization variable, initialized to zero
  Eigen::MatrixXd P;
  Eigen::MatrixXd Q;

  // Old versions (used for back-tracking and fusion counting)
  Eigen::Index nzeros_row_old;
  Eigen::Index nzeros_col_old;
  Eigen::MatrixXd U_old;
  Eigen::MatrixXd P_old;
  Eigen::MatrixXd Q_old;
  Eigen::MatrixXd V_row_old;
  Eigen::MatrixXd Z_row_old;
  Eigen::MatrixXd V_col_old;
  Eigen::MatrixXd Z_col_old;
  Eigen::ArrayXi  v_row_zeros_old;
  Eigen::ArrayXi  v_col_zeros_old;

  // Internal storage buffers
  Eigen::Index buffer_size;
  Eigen::Index storage_index;
  Eigen::MatrixXd UPath;
  Eigen::MatrixXd V_rowPath;
  Eigen::MatrixXd V_colPath;
  Eigen::VectorXd gamma_path;
  Eigen::MatrixXi v_row_zeros_path;
  Eigen::MatrixXi v_col_zeros_path;
};

#endif
