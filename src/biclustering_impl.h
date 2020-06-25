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
    sp(show_progress_, D_row_.rows() + D_col_.cols()),
    DDT_col(D_col_ * D_col_.transpose()),
    DTD_row(D_row_.transpose() * D_row_) {

      // Set initial values for optimization variables
      U = X;
      V_row = D_row * U;
      Z_row = Eigen::MatrixXd::Zero(V_row.rows(), V_row.cols());
      V_col = U * D_col;
      Z_col = Eigen::MatrixXd::Zero(V_col.rows(), V_col.cols());


      v_row_zeros = Eigen::ArrayXi::Zero(num_row_edges);
      v_col_zeros = Eigen::ArrayXi::Zero(num_col_edges);
      gamma = 0;

      //compute alpha
      //TODO: implement a tighter alpha calculation
      double row_max_deg = D_row.cwiseAbs().colwise().sum().maxCoeff();
      double col_max_deg = D_col.cwiseAbs().rowwise().sum().maxCoeff();
      alpha = 2 * (row_max_deg  * col_max_deg);


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

  void admm_step(){
    // U-update
    U = (X + alpha * U + rho * (
        D_row.transpose() * (V_row - Z_row) +
        (V_col - Z_col) * D_col.transpose() -
        DTD_row * U -
        U * DDT_col
      )) / (1 + alpha);

    Eigen::MatrixXd DrowU = D_row * U;
    Eigen::MatrixXd UDcol = U * D_col;
    ClustRVizLogger::debug("U = ") << U;

    // V-updates
    Eigen::MatrixXd DUZ = DrowU + Z_row; //DUZ = D_row * U + Z_row
    V_row = MatrixRowProx(DUZ, gamma / rho, weights_row, l1);
    ClustRVizLogger::debug("V_row = ") << V_row;


    Eigen::MatrixXd UDZ = UDcol + Z_col; //UDZ = (U * D_col + Z_col
    V_col = MatrixRowProx(UDZ, gamma / rho, weights_col, l1);
    ClustRVizLogger::debug("V_col = ") << V_col;


    // Z-updates
    Z_row = Z_row + DrowU - V_row;
    ClustRVizLogger::debug("Z_row = ") << Z_row;

    Z_col = Z_col + UDcol - V_col;
    ClustRVizLogger::debug("Z_col = ") << Z_col;


    // Identify row fusions (rows of V_row which have gone to zero)
    Eigen::VectorXd v_row_norms = V_row.rowwise().squaredNorm();

    for(Eigen::Index i = 0; i < num_row_edges; i++){
      v_row_zeros(i) = v_row_norms(i) == 0;
    }

    nzeros_row = v_row_zeros.sum();

    // Identify column fusions (rows of V_col which have gone to zero)
    Eigen::VectorXd v_col_norms = V_col.colwise().squaredNorm();

    for(Eigen::Index i = 0; i < num_col_edges; i++){
      v_col_zeros(i) = v_col_norms(i) == 0;
    }

    nzeros_col = v_col_zeros.sum();

    double loss;
    if (l1) {
      loss = 0.5 * (X - U).squaredNorm() + gamma * (
        V_row.cwiseAbs().rowwise().sum().dot(weights_row) +
        V_col.cwiseAbs().colwise().sum().dot(weights_col));
    } else {
      loss = 0.5 * (X - U).squaredNorm() + gamma * (
        V_row.rowwise().norm().dot(weights_row) +
        V_col.colwise().norm().dot(weights_col));
    }
    ClustRVizLogger::info("Objective function: ") <<  loss;  


    ClustRVizLogger::debug("Number of row fusions identified ") << nzeros_row;
    ClustRVizLogger::debug("Number of column fusions identified ") << nzeros_col;
  }


  void save_fusions(){
    nzeros_row_old = nzeros_row;
    nzeros_col_old = nzeros_col;
  }

  void save_old_values(){
    U_old = U;

    V_row_old = V_row;
    Z_row_old = Z_row;
    v_row_zeros_old = v_row_zeros;

    V_col_old = V_col;
    Z_col_old = Z_col;
    v_col_zeros_old = v_col_zeros;
  }

  void load_old_variables(){
    U = U_old;

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
    return scaled_squared_norm(U - U_old) < 5e-14 && 
            scaled_squared_norm(Z_row - Z_row_old) +
            scaled_squared_norm(Z_col - Z_col_old) +
            scaled_squared_norm(V_row - V_row_old) +
            scaled_squared_norm(V_col - V_col_old) < 5e-14;
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

  void tick(unsigned int iter){
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
  double alpha; 
  // Theoretically, it's part of the algorithm, not the problem
  // but we need it in the steps...
  bool  l1;         // Is the L1 (true) or L2 (false) norm being used?
  const int n;      // Problem dimensions
  const int p;
  const int num_row_edges;
  const int num_col_edges;

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

  // Precomputed products that are reused in U-update
  const Eigen::MatrixXd DDT_col; // D_col * D_col^T
  const Eigen::MatrixXd DTD_row; // D_row^T * D_row

  // Old versions (used for back-tracking and fusion counting)
  Eigen::Index nzeros_row_old;
  Eigen::Index nzeros_col_old;
  Eigen::MatrixXd U_old;
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
