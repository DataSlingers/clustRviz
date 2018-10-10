# Helper types to hold onto weight selection details
RBFWeights <- function(phi = phi,
                       user_phi = user_phi,
                       dist.method = dist.method,
                       p = p){

  obj <- list(phi         = phi,
              user_phi    = user_phi,
              dist.method = dist.method,
              p           = p,
              name        = "RBF Kernel")

  class(obj) <- c("RBFWeights", "ClusteringWeights")

  obj
}

UserFunction <- function(){
  obj <- list(name = "User-Provided Function")
  class(obj) <- "ClusteringWeights"
  obj
}

UserMatrix <- function(){
  obj <- list(name = "User-Provided Matrix")
  class(obj) <- "ClusteringWeights"
  obj
}

add_sparse_weights <- function(fit_type, user_k, k){
  fit_type$k      <- k
  fit_type$user_k <- user_k
  class(fit_type) <- c(class(fit_type), "SparseClusteringWeights")

  fit_type
}

#' @export
print.RBFWeights <- function(x, indent = 2, ...){
  cat(" - Source: Radial Basis Function Kernel Weights\n")
  cat(" - Distance Metric: ", capitalize_string(x$dist.method), sep = "")
  if (x$dist.method == "minkowski") {
    cat(" (p = ", round(x$p, 3), ")")
  }
  cat("\n")
  cat(" - Scale parameter (phi): ", round(x$phi, 3), " [",
      if (x$user_phi) "User-Supplied" else "Data-Driven", "]\n", sep = "")

  if (inherits(x, "SparseClusteringWeights")) {
    cat(" - Sparsified: ", x$k, " Nearest Neighbors [",
        if (x$user_k) "User-Supplied" else "Data-Driven", "]\n", sep = "")
  }

  cat("\n")
}

#' @export
print.ClusteringWeights <- function(x, indent = 2, ...){
  cat(" - Source: ", x$name, "\n", sep = "")
  if (inherits(x, "SparseClusteringWeights")) {
    cat(" - Sparsified: ", x$k, " Nearest Neighbors [",
        if (x$user_phi) "User-Supplied" else "Data-Driven", "]\n", sep = "")
  }
  cat("\n")
}

#' Factory function - this takes a weight function and wraps it up in some
#' sparsification code. Not intended to be exposed to users, but useful to
#' create sparse analogues of various weight schemes
#' @noRd
#' @importFrom utils modifyList
make_sparse_weights_func <- function(weight_func){

  function(..., k = "auto"){
    dense_weight_func <- weight_func(...)

    function(X){
      dense_fit     <- dense_weight_func(X)
      dense_weights <- dense_fit$weight_mat
      fit_type      <- dense_fit$type

      check_weight_matrix(dense_weights)

      user_k <- (k != "auto")

      if (k == "auto") {
        user_k <- FALSE
        k <- 1

        while (TRUE) {
          adjacency_matrix <- take_k_neighbors(dense_weights, k = k)
          connected <- is_connected_adj_mat(adjacency_matrix)

          if (connected) {
            break
          } else {
            k <- k + 1
          }

          if(k == NCOL(adjacency_matrix)){
            crv_error("Cannot find ", sQuote("k"), " yielding fully connected graph.")
          }
        }
      }

      if (!is_integer_scalar(k)) {
        crv_error("If not `auto,` ", sQuote("k"), " must be an integer scalar (vector of length 1).")
      }

      if (k <= 0) {
        crv_error(sQuote("k"), " must be positive.")
      }

      sparse_weights <- take_k_neighbors(dense_weights, k = k)

      connected <- is_connected_adj_mat(sparse_weights)

      if (!connected) {
        crv_error("k = ", k, " does not give a fully connected graph. Convex (bi)clustering will not converge.")
      }

      list(weight_mat = sparse_weights,
           type = add_sparse_weights(fit_type,
                                     user_k = user_k,
                                     k      = k))
    }
  }
}

#' Construct Clustering Weights Based on the Radial Basis Function (Gaussian/Euclidean) Kernel
#'
#' This is a \emph{factory function} - it returns a \emph{function} which can be
#' used to create a matrix of clustering weights. In particular, it returns a function
#' which takes a n-by-p data matrix \eqn{X} and returns an n-by-n matrix whose
#' \eqn{(i, j)}-th element is given by \eqn{e^{-phi * dist(x_i, x_j)}} where
#' \eqn{x_i}, \eqn{x_j} are the \eqn{i}-th and \eqn{j}-th row of \eqn{X}
#' respectively. The distance metric used is determined by the
#' \code{dist.method} and \code{p} arguments, which are passed to
#' \code{\link[stats]{dist}}.
#'
#' The sparse weights (\code{sparse_rbf_kernel_weights}) are calculated by
#' dropping all but the \code{k} largest weights for each row of the matrix
#' (equivalent to taking the \code{k} nearest neighbors to each point). The
#' weight matrix is symmetrized, so if \emph{a} is a neighbor of \emph{b}, but
#' not vice versa, the edge is still included. If \code{k} is too small, resulting
#' in a non-fully-connected graph, an error is thrown.
#'
#' If \code{phi == "auto"}, a grid of possible \eqn{phi} values are used and
#' the \code{phi} which maximizes the variance of the resulting weights is taken.
#'
#' If \code{k == "auto"}, the smallest \code{k} that still yields a fully connected
#' graph is used.
#'
#' @param phi The scale factor used for the RBF kernel
#' @param dist.method The type of distance used to calculate distances between
#'                    points. See the \code{method} argument to \code{\link[stats]{dist}}.
#' @param p The power of the Minkowski distance (only relevant if \code{method == "minkowski"}).
#'          See the \code{p} argument to \code{\link[stats]{dist}}.
#' @importFrom stats dist var
#' @return A function which, when called, returns a matrix of clustering weights.
#' @examples
#' weight_func <- dense_rbf_kernel_weights()
#' weight_func(presidential_speech)
#'
#' weight_func <- dense_rbf_kernel_weights(phi=0.1, dist.method="canberra")
#' weight_func(presidential_speech)
#'
#' weight_func <- sparse_rbf_kernel_weights()
#' weight_func(presidential_speech)
#' @export
#' @name RBF Kernel Weights
#' @rdname rbf_kernel_weights
dense_rbf_kernel_weights <- function(phi = "auto",
                                     dist.method = c("euclidean",
                                                     "maximum",
                                                     "manhattan",
                                                     "canberra",
                                                     "binary",
                                                     "minkowski"),
                                     p = 2){

  tryCatch(dist.method <- match.arg(dist.method),
           error = function(e){
             crv_error("Unsupported choice of ", sQuote("weight.dist;"),
                       " see the ", sQuote("method"), " argument of ",
                       sQuote("stats::dist"), " for supported distances.")
           })

  if ((p <= 0) || !is_numeric_scalar(p)) {
    crv_error(sQuote("p"),
              " must be a positive scalar; see the ", sQuote("p"),
              " argument of ", sQuote("stats::dist"), " for details.")
  }

  function(X){
    user_phi <- (phi != "auto")

    if (phi == "auto") {
      ## FIXME - This leads to calling dist() one more time than strictly
      ##         necessary...
      phi_range <- 10^(seq(-10, 10, length.out = 21))
      weight_vars <- vapply(phi_range,
                            function(phi) var(exp((-1) * phi * (dist(X, method = dist.method, p = p)[TRUE])^2)),
                            numeric(1))

      phi <- phi_range[which.max(weight_vars)]
    }

    if (!is_numeric_scalar(phi)) {
      crv_error("If not `auto,` ", sQuote("phi"), " must be a numeric scalar (vector of length 1).")
    }

    if (phi <= 0) {
      crv_error(sQuote("phi"), " must be positive.")
    }

    dist_mat <- as.matrix(dist(X, method = dist.method, p = p))
    dist_mat <- exp(-1 * phi * dist_mat^2)

    check_weight_matrix(dist_mat)

    diag(dist_mat) <- 0

    list(weight_mat = dist_mat,
         type = RBFWeights(phi = phi,
                           user_phi = user_phi,
                           dist.method = dist.method,
                           p = p))
  }
}

#' @export
#' @rdname rbf_kernel_weights
#' @param ... Arguments passed through from \code{sparse_rbf_kernel_weights} to
#'            \code{dense_rbf_kernel_weights}
#' @param k The number of neighbors to use
sparse_rbf_kernel_weights <- make_sparse_weights_func(dense_rbf_kernel_weights)

#' Check if an adjacency matrix encodes a connected graph.
#'
#' We re-use our cluster assignment code: if all the edges are "on" and imply
#' that all the points are clustered together, then we have a fully connected graph.
#'
#' @noRd
#' @importFrom Matrix nnzero
#' @importMethodsFrom Matrix t
is_connected_adj_mat <- function(adjacency_matrix){
  adjacency_matrix_ut <- adjacency_matrix * upper.tri(adjacency_matrix);

  edge_list <- which(adjacency_matrix_ut != 0, arr.ind = TRUE)

  get_cluster_assignments(edge_list,
                          matrix(TRUE, ncol = NROW(edge_list), nrow = 1),
                          NROW(adjacency_matrix))[[1]]$no == 1
}

#' @noRd
#' Get the k-th largest element of a vector
#'
#' See https://stackoverflow.com/a/2453619/ and comments thereon
kth_largest_element <- function(x, k){
  -sort(-x, partial = k, decreasing = FALSE)[k]
}

#' @noRd
#' Sparsify a graph (represented by its weighted adjacency matrix), leaving only
#' its k nearest neighbors.
#'
#' We symmetrize the graph, so if a is a nearest-neighbor of b OR if b is
#' a nearest-neighbor of a, the edge is retained
take_k_neighbors <- function(weight_matrix, k){
  if (k == NROW(weight_matrix)) {
    crv_error("k should be less than the size of the graph")
  }
  mask <- apply(weight_matrix, 2, function(x) x >= kth_largest_element(x, k))
  mask <- (mask | t(mask)) # Symmetrize our matrix

  weight_matrix * mask
}

#' @noRd
#' Convert a weight matrix to the vectorized form used by other functions
weight_mat_to_vec <- function(weight_mat){
  t(weight_mat)[lower.tri(weight_mat, diag = FALSE)]
}
