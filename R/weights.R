#' Compute the minimum number of nearest neighbors required to fully
#' cluster all observations
#'
#' @importFrom Matrix which
#' @param X an n.obs x p.vars matrix
#' @param dense.weights a vector vector of dense weights such as returned by
#' \code{DenseWeights}
#' @return k an integer. The smallest number of nearest neighbors required to
#' fully cluster all observations.
#' @keywords internal
#' @examples
#' \dontrun{
#' library(clustRviz)
#' data("presidential_speech")
#' Xdat <- presidential_speech$X
#' Xdat.preprocessed <- scale(Xdat,center=TRUE,scale=FALSE)
#' dense.weights <- DenseWeights(X = Xdat.preprocessed,phi=1e-3)
#' MinKNN(
#' X = Xdat.preprocessed,
#' dense.weights = dense.weights
#' ) -> k.min
#' }
MinKNN <- function(X, dense.weights) {
  n <- nrow(X)
  weight.adj <- WeightAdjacency(dense.weights, n)
  cardE <- sum(weight.adj)
  E <- Matrix::which(weight.adj != 0, arr.ind = TRUE)
  E <- E[order(E[, 1], E[, 2]), ]
  k <- 0
  n.comp <- 2
  while (n.comp != 1) {
    k <- k + 1
    w.sp <- SparseWeights(X = X, dense.weights = dense.weights, k = k)
    CreateAdjacency(E, sp.pattern = as.numeric(w.sp != 0), n = n) -> adj.full
    CreateClusterGraph(adj.full) -> cluster.full
    n.comp <- GetClusters(cluster.full)$no
  }
  k
}

#' Compute dense distanced-based gaussian kernel weights for use in CARP or CBASS
#'
#' @param X an n by p matrix with rows the observations and columns the variables.
#' @param phi a number. Scaling parameter used in the gaussian kernel.
#' @param method a string. Passed to the \code{dist} function. See \code{?dist}
#' @param p The power of Minkowski distance. Passed to the \code{dist} function. See \code{?dist}
#' @return dense.weights a numeric vector of weights
#' @importFrom stats dist
#' @keywords internal
#' @examples
#' \dontrun{
#' library(clustRviz)
#' data("presidential_speech")
#' Xdat <- presidential_speech$X
#' Xdat.preprocessed <- scale(Xdat,center=TRUE,scale=FALSE)
#' dense.weights <- DenseWeights(X = Xdat.preprocessed,phi=1e-3)
#' }
DenseWeights <- function(X,
                         phi = 1,
                         method = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
                         p = 2) {
  method <- match.arg(method)
  exp((-1) * phi * (stats::dist(X, method = method, p = p)[TRUE])^2)
}

#' (slowly) Compute sparse kNN weight vector from dense weights
#'
#' @param X an n by p matrix with rows the observations and columns the variables.
#' @param dense.weights a vector dense weights, such as those returned
#' by \code{DenseWeights}
#' @param k a positive integer less than nobs. The number of nearest neighbors.
#' @return sparse.weights a numeric vector of sparse weights
#' @keywords internal
#' @importFrom FNN get.knn
#' @importFrom utils combn
#' @importFrom dplyr %>%
SparseWeights_LEGACY <- function(X, dense.weights, k) {
  nobs <- nrow(X)
  tmp <- FNN::get.knn(X, k = k, algorithm = "brute")
  lapply(1:nobs, function(ind) {
    vec <- rep(0, times = nobs)
    vec[tmp$nn.index[ind, ]] <- 1
    vec[ind] <- 1
    vec
  }) %>%
    do.call(rbind, .) -> mat.pat

  apply(utils::combn(nobs, 2), 2, function(pair.ind) {
    vec1 <- mat.pat[pair.ind[1], ]
    vec2 <- mat.pat[pair.ind[2], ]
    andvec <- vec1 * vec2
    sum(as.logical(andvec[pair.ind[1]]) | as.logical(andvec[pair.ind[2]]))
  }) -> sp.mask
  dense.weights * sp.mask
}

#' Compute sparse kNN weight vector from dense weights
#'
#' @param X an n by p matrix with rows the observations and columns the variables.
#' @param dense.weights a vector dense weights, such as those returned
#' by \code{DenseWeights}
#' @param k a positive integer less than nobs. The number of nearest neighbors.
#' @return sparse.weights a numeric vector of sparse weights
#' @keywords internal
#' @importFrom FNN get.knn
#' @importFrom utils combn
#' @importFrom dplyr %>%
SparseWeights <- function(X, dense.weights, k) {
  nobs <- nrow(X)
  tmp <- FNN::get.knn(X, k = k, algorithm = "brute")
  lapply(1:nobs, function(ind) {
    vec <- rep(0, times = nobs)
    vec[tmp$nn.index[ind, ]] <- 1
    vec[ind] <- 1
    vec
  }) %>%
    do.call(rbind, .) -> mat.pat


  nnvec1 <- mat.pat
  nnvec1[lower.tri(nnvec1, diag = TRUE)] <- NA
  nnvec1 <- t(nnvec1)
  nnvec1.final <- nnvec1[TRUE][!is.na(nnvec1[TRUE])]
  nnvec1.final

  nnvec2 <- t(mat.pat)
  nnvec2[lower.tri(nnvec2, diag = TRUE)] <- NA
  nnvec2 <- t(nnvec2)
  nnvec2.final <- nnvec2[TRUE][!is.na(nnvec2[TRUE])]
  nnvec2.final

  sp.mask <- as.numeric(as.logical(nnvec1.final) | as.logical(nnvec2.final))
  dense.weights * sp.mask
}

#' Construct adjacency matrix induced by weights
#'
#' @param weights a vector of weights such as returned by \code{SparseWeights}
#' @param nobs the number of observations being clustered
#' @param weighted a logical. If \code{FALSE} created unweighted adjacency matrix
#' determined by support of weight vector. If TRUE create weighted adjacency;
#' default is FALSE.
#' @param upper a logical. If \code{TRUE} only return upper triangular of matrix.
#' If \code{FALSE} return symmetric adjacency
#' @return adj a sparse adjacency matrix
#' @keywords internal
#' @importFrom Matrix Matrix
#' @importFrom Matrix t
WeightAdjacency <- function(weights, nobs, weighted = FALSE, upper = TRUE) {
  adj <- Matrix::Matrix(data = 0, nrow = nobs, ncol = nobs, sparse = TRUE)
  if (!weighted) {
    adj[lower.tri(adj, diag = FALSE)] <- as.numeric(weights != 0)
  } else {
    adj[lower.tri(adj, diag = FALSE)] <- weights
  }
  adj <- Matrix::t(adj)
  if (!upper) {
    adj <- adj + Matrix::t(adj)
  }
  adj
}

#' Distances supported by \code{\link[stats]{dist}}
#' @noRd
SUPPORTED_DISTANCES <- c("euclidean",
                         "maximum",
                         "manhattan",
                         "canberra",
                         "binary",
                         "minkowski")

#' Factory function - this takes a weight function and wraps it up in some
#' sparsification code. Not intended to be exposed to users, but useful to
#' create sparse analogues of various weight schemes
#' @noRd
make_sparse_weights_func <- function(weight_func){

  function(..., k = "auto"){
    dense_weight_func <- weight_func(...)

    function(X){
      dense_weights <- dense_weight_func(X)

      if(k == "auto"){
        k <- 1

        while(TRUE){
          adjacency_matrix <- take_k_neighbors(dense_weights, k=k)
          connected <- is_connected_adj_mat(adjacency_matrix)

          if(connected){
            break
          } else {
            k <- k + 1
          }
        }

      }

      if(!is_integer_scalar(k)){
        stop("If not `auto,` ", sQuote("k"), " must be an integer scalar (vector of length 1).")
      }

      if(k <= 0){
        stop(sQuote("k"), " must be positive.")
      }

      sparse_weights <- take_k_neighbors(dense_weights, k = k)

      connected <- is_connected_adj_mat(sparse_weights)

      if(!connected){
        stop("k = ", k, " does not give a fully connected graph. Convex (bi)clustering will not converge.")
      }

      sparse_weights
    }
  }
}

#' Construct clustering weights based on the Gaussian (RBF) kernel
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
#' The sparse weights (\code{sparse_gaussian_kernel_weights}) are calculated by
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
#' weight_func <- dense_gaussian_kernel_weights()
#' weight_func(presidential_speech)
#'
#' weight_func <- dense_gaussian_kernel_weights(phi=1, dist.method="canberra")
#' weight_func(presidential_speech)
#'
#' weight_func <- sparse_gaussian_kernel_weights()
#' weight_func(presidential_speech)
#' @export
#' @name RBF Kernel Weights
#' @rdname rbf_kernel_weights
dense_gaussian_kernel_weights <- function(phi = "auto",
                                          dist.method = c("euclidean",
                                                          "maximum",
                                                          "manhattan",
                                                          "canberra",
                                                          "binary",
                                                          "minkowski"),
                                          p = 2){

  tryCatch(dist.method <- match.arg(dist.method),
           error=function(e){
             stop("Unsupported choice of ", sQuote("weight.dist;"),
                  " see the ", sQuote("method"), " argument of ",
                  sQuote("stats::dist"), " for supported distances.")
           })

  if ((p <= 0) || !is_numeric_scalar(p)) {
    stop(sQuote("p"),
         " must be a positive scalar; see the ", sQuote("p"),
         " argument of ", sQuote("stats::dist"), " for details.")
  }

  function(X){
    if(phi == "auto"){
      ## FIXME - This leads to calling dist() one more time than strictly
      ##         necessary...
      phi_range <- 10^(seq(-10, 10, length.out=21))
      weight_vars <- vapply(phi_range,
                            function(phi) var(exp((-1) * phi * (dist(X, method = dist.method, p = p)[TRUE])^2)),
                            numeric(1))

      phi <- phi_range[which.max(weight_vars)]
    }

    if(!is_numeric_scalar(phi)){
      stop("If not `auto,` ", sQuote("phi"), " must be a numeric scalar (vector of length 1).")
    }

    if(phi <= 0){
      stop(sQuote("phi"), " must be positive.")
    }

    dist_mat <- as.matrix(dist(X, method = dist.method, p = p))
    dist_mat <- exp(-1 * phi * dist_mat^2)

    diag(dist_mat) <- 0

    dist_mat
  }
}

#' @export
#' @rdname rbf_kernel_weights
#' @param ... Arguments passed through from \code{sparse_gaussian_kernel_weights} to
#'            \code{dense_gaussian_kernel_weights}
#' @param k The number of neighbors to use
sparse_gaussian_kernel_weights <- make_sparse_weights_func(dense_gaussian_kernel_weights)

#' Check if an adjacency matrix encodes a connected graph.
#'
#' We use the fact that that if \eqn{A} is the adjacency matrix of a graph,
#' then \eqn{(A^k)_{ij}} is non-zero if and only if there is a path of exactly length
#' \eqn{k} from $i$ to $j$. If we include the diagonal of \eqn{A}, then
#' \eqn{(A^k)_{ij}} is non-zero if and only if there is a path of length at most
#' \eqn{k} from $i$ to $j$
#'
#' We only deal with symmetric adjacency matrices (undirected graphs)
#'
#' We square \eqn{A} repeatedly (yielding \eqn{A}, \eqn{A^2}, \eqn{A^4}, \eqn{A^8})
#' checking for a fully connected (entirely non-zero) graph, or stopping when
#' we have exceeded the dimension of \eqn{A}.
#'
#' @noRd
#' @importFrom Matrix nnzero
#' @importMethodsFrom Matrix t
is_connected_adj_mat <- function(adjacency_matrix){
  diag(adjacency_matrix) <- 1

  N <- NCOL(adjacency_matrix)

  n <- 1
  num_connected <- nnzero(adjacency_matrix)

  if (num_connected == N^2) {
    return(TRUE)
  }

  while (n < N) {
    adjacency_matrix <- crossprod(adjacency_matrix)
    n <- 2 * n
    num_connected_old <- num_connected
    num_connected <- nnzero(adjacency_matrix)

    if (num_connected_old == num_connected) {
      return(FALSE)
    }

    if (num_connected == N^2) {
      return(TRUE)
    }
  }

  return(nnzero(adjacency_matrix) == N^2)
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
    stop("k should be less than the size of the graph")
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
