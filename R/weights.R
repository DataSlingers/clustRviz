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
