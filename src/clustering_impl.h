#ifndef CLUSTRVIZ_CLUSTERING_H
#define CLUSTRVIZ_CLUSTERING_H 1

#include "clustRviz_base.h"
#include "clustRviz_logging.h"
#include "status.h"

template <class DataType>
class ConvexClusteringSkeleton {
  using MatrixXt = Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic>;
  using VectorXt = Eigen::Matrix<DataType, Eigen::Dynamic, 1>;

public:
  double gamma; // Current regularization level - need to be able to manipulate this externally

  ConvexClusteringSkeleton(const Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic>& X_,
                           const Eigen::ArrayXXd& M_,
                           const Eigen::MatrixXd& D_,
                           const Eigen::VectorXd& weights_,
                           const double rho_,
                           const bool l1_,
                           const bool show_progress_):
  X(X_),
  M(M_),
  D(D_),
  weights(weights_),
  rho(rho_),
  l1(l1_),
  n(X_.rows()),
  p(X_.cols()),
  num_edges(D_.rows()),
  sp(show_progress_, D_.rows()) {

    // Set initial values for optimization variables
    U = X;
    V = D * U;
    Z = V;
    v_zeros = Eigen::ArrayXi::Zero(num_edges);
    gamma = 0;

    sp.set_v_norm_init(V.squaredNorm());

    // Initialize storage buffers
    buffer_size = 1.5 * n;
    UPath.resize(n * p, buffer_size);
    VPath.resize(p * num_edges, buffer_size);
    gamma_path.resize(buffer_size);
    v_zeros_path.resize(num_edges, buffer_size);

    // Store initial values
    nzeros = 0;
    storage_index = 0;
    store_values();

    // PreCompute chol(I + rho D^TD) for easy inversions in the U update step
    MatrixXt IDTD = rho * D.transpose() * D + MatrixXt::Identity(n, n);
    u_step_solver.compute(IDTD);
  };

  bool is_interesting_iter(){
    // FIXME? A better check would be v_zeros != v_zeros_old for fusion IDs
    // not just number
    return nzeros != nzeros_old;
  }

  bool multiple_fusions(){
    return nzeros > nzeros_old + 1;
  }

  bool is_complete(){
    // Stop when all edges are active (when we have a connected graph, all vertices are fused)
    return nzeros == num_edges;
  }

  inline void u_step(){
    // U-update -- Make this a method so time series clustering can replace
    MatrixXt X_imputed = M * X.array() + (1.0 - M) * U.array();
    U = u_step_solver.solve(X_imputed + rho * D.transpose() * (V - Z));

    ClustRVizLogger::debug("U = ") << U;
  }

  void admm_step(){
    u_step();
    // V-update
    const MatrixXt DUZ = D * U + Z;
    V = MatrixRowProx(DUZ, gamma / rho, weights, l1);
    ClustRVizLogger::debug("V = ") << V;

    // Z-update
    Z = DUZ - V;
    ClustRVizLogger::debug("Z = ") << Z;

    // Identify cluster fusions (rows of V which have gone to zero)
    Eigen::VectorXd v_norms = V.rowwise().squaredNorm();

    for(Eigen::Index i = 0; i < num_edges; i++){
      v_zeros(i) = v_norms(i) == 0;
    }

    nzeros = v_zeros.sum();

    ClustRVizLogger::debug("Number of fusions identified ") << nzeros;
  }

  void save_fusions(){
    nzeros_old = nzeros;
  }

  void save_old_values(){
    U_old = U;
    V_old = V;
    Z_old = Z;
    v_zeros_old = v_zeros;
  }

  void load_old_variables(){
    U = U_old;
    V = V_old;
    Z = Z_old;
  }

  void load_old_fusions(){
    v_zeros = v_zeros_old;
  }

  bool has_fusions(){
    return nzeros > 0;
  }

  bool admm_converged(double thresh = CLUSTRVIZ_DEFAULT_STOP_PRECISION){
    return (scaled_squared_norm(U - U_old) < thresh) &&
           (scaled_squared_norm(V - V_old) < thresh) &&
           (scaled_squared_norm(Z - Z_old) < thresh);
  }

  void store_values(){
    if(storage_index >= buffer_size){
      ClustRVizLogger::info("Resizing storage from ") << buffer_size << " to " << 2 * buffer_size << " iterations.";
      buffer_size *= 2; // Double our buffer sizes
      UPath.conservativeResize(UPath.rows(), buffer_size);
      VPath.conservativeResize(VPath.rows(), buffer_size);
      gamma_path.conservativeResize(buffer_size);
      v_zeros_path.conservativeResize(v_zeros_path.rows(), buffer_size);
    }

    // Store values
    UPath.col(storage_index)        = Eigen::Map<VectorXt>(U.data(), n * p);
    VPath.col(storage_index)        = Eigen::Map<VectorXt>(V.data(), p * num_edges);
    gamma_path(storage_index)       = gamma;
    v_zeros_path.col(storage_index) = v_zeros;

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
    VPath.conservativeResize(VPath.rows(), storage_index);
    gamma_path.conservativeResize(storage_index);
    v_zeros_path.conservativeResize(v_zeros_path.rows(), storage_index);

    return Rcpp::List::create(Rcpp::Named("u_path")      = UPath,
                              Rcpp::Named("v_path")      = VPath,
                              Rcpp::Named("v_zero_inds") = v_zeros_path,
                              Rcpp::Named("gamma_path")  = gamma_path);
  }

  void tick(unsigned int iter){
    sp.update(nzeros, V.squaredNorm(), iter, gamma);
  }

protected:
  // Fixed (non-data-dependent) problem details
  const MatrixXt& X; // Data matrix (to be clustered)
  const Eigen::ArrayXXd& M; // Missing data mask
  const Eigen::MatrixXd& D; // Edge (differencing) matrix
  const Eigen::VectorXd& weights; // Clustering weights
  const double rho; // ADMM relaxation parameter -- TODO: Factor this out?
                    // Theoretically, it's part of the algorithm, not the problem
                    // but we need it in the steps...
  bool  l1;         // Is the L1 (true) or L2 (false) norm being used?
  const int n;      // Problem dimensions
  const int p;
  const int num_edges;
  Eigen::LLT<MatrixXt> u_step_solver; // Cached factorization for u-update

  // Progress printer
  StatusPrinter sp;

  // Current copies of ADMM variables
  MatrixXt U; // Primal variable
  MatrixXt V; // Split variable
  MatrixXt Z; // Dual variable
  Eigen::ArrayXi v_zeros; // Fusion indicators
  Eigen::Index nzeros; // Number of fusions

  // Old versions (used for back-tracking and fusion counting)
  Eigen::Index nzeros_old;
  MatrixXt U_old;
  MatrixXt V_old;
  MatrixXt Z_old;
  Eigen::ArrayXi  v_zeros_old;

  // Internal storage buffers
  Eigen::Index buffer_size;
  Eigen::Index storage_index;
  MatrixXt UPath;
  MatrixXt VPath;
  Eigen::VectorXd gamma_path;
  Eigen::MatrixXi v_zeros_path;
};

#endif
