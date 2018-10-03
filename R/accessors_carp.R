## Accessor functions for CARP

#' Get Clustering Results for \code{CARP}
#'
#' \code{get_cluster_labels} returns a factor vector of cluster labels.
#' \code{get_clustered_data} returns a matrix (with the same dimensions and names
#' as the original data), but with the values for each observation replaced by
#' its "estimated" value (\emph{i.e.}, the appropriate cluster centroid).
#' \code{get_cluster_centroids} returns a \code{k}-by-\code{p} matrix of cluster
#' centroids, with the same column names as the original data.
#'
#' @param x An object of class \code{CARP} as produced by \code{\link{CARP}}
#' @param percent A number between 0 and 1, giving the regularization level (as
#'                a fraction of the final regularization level used) at which to
#'                get cluster labels.
#' @param k The desired number of clusters. If no iteration with exactly this
#'          many clusters is found, the first iterate with fewer than \code{k}
#'          clusters is used.
#' @param refit Should "naive" centroids (\code{TRUE}) or the actual centroids
#'              estimated by convex clustering be used? The default (\code{refit = TRUE})
#'              centroids returned are actual centroids (mean) of all elements
#'              assigned to that cluster; if \code{refit = FALSE}, the \eqn{\hat{U}}
#'              from the convex clustering problem is used. Due to the global
#'              shrinkage imposed, these clusters are more "shrunk together" than
#'              the naive clusters.
#' @param ... Additional arguments - if any are provided, an error is signalled.
#' @details \code{get_clustered_data} and \code{get_cluster_centroids} return
#' centroids on the original scale of the data, independent of any pre-processing
#' flags passed to \code{CARP}. Note that exactly one of \code{percent} and
#' \code{k} must be supplied to each function.
#' @examples
#' carp_fit <- CARP(presidential_speech)
#'
#' # Get clustering results from 50% along the path
#' get_cluster_labels(carp_fit, percent = 0.5)
#'
#' # Get labels corresponding to the 3 cluster solution
#' get_cluster_labels(carp_fit, k = 3)
#'
#' # Get 3 cluster centroids
#' get_cluster_centroids(carp_fit, k = 3)
#'
#' # Get the clustered estimates for k = 3 clusters
#' get_clustered_data(carp_fit, k = 3)
#' @export
#' @rdname accessors_carp
get_cluster_labels <- function(x, ...){
  UseMethod("get_cluster_labels")
}

#' @importFrom dplyr select filter summarize pull arrange
#' @importFrom rlang .data
#' @rdname accessors_carp
#' @export
get_cluster_labels.CARP <- function(x, ..., percent, k){

  dots <- list(...)
  if ( length(dots) != 0) {
    if (!is.null(names(dots))) {
      nm <- names(dots)
      nm <- nm[nzchar(nm)]
      crv_error("Unknown argument ", sQuote(nm[1]), " passed to ", sQuote("get_cluster_labels."))
    } else {
      crv_error("Unknown argument passed to ", sQuote("get_cluster_labels."))
    }
  }

  has_percent <- !missing(percent)
  has_k       <- !missing(k)
  n_args      <- has_percent + has_k

  if(n_args != 1){
    crv_error("Exactly one of ", sQuote("percent"), " and ", sQuote("k"), " must be supplied.")
  }

  if(has_k){

    if ( !is_integer_scalar(k) ){
      crv_error(sQuote("k"), " must be an integer scalar (vector of length 1).")
    }

    if( k <= 0 ) {
      crv_error(sQuote("k"), " must be positive.")
    }

    if( k > NROW(x$X) ){
      crv_error(sQuote("k"), " cannot be more than the observations in the original data set (", NROW(x$X), ").")
    }

    percent <- x$carp.cluster.path.vis %>%
                 select(.data$LambdaPercent, .data$NCluster) %>%
                 filter(.data$NCluster <= k) %>%
                 select(.data$LambdaPercent) %>%
                 summarize(percent = min(.data$LambdaPercent)) %>%
                 pull
  }

  if( !is_percent_scalar(percent) ){
    crv_error(sQuote("percent"), " must be a scalar between 0 and 1 (inclusive).")
  }

  cluster_labels_df <- x$carp.cluster.path.vis %>%
                         select(.data$LambdaPercent,
                                .data$ObsLabel,
                                .data$Obs,
                                .data$Cluster,
                                .data$Iter) %>%
                         filter(.data$LambdaPercent >= percent) %>%
                         filter(.data$LambdaPercent == min(.data$LambdaPercent)) %>%
                         filter(.data$Iter == min(.data$Iter)) %>% # In case we have multiple iterations at same lambda
                         arrange(.data$Obs) %>%
                         select(-.data$Obs)

  n_clusters <- num_unique(cluster_labels_df$Cluster)
  labels <- factor(x = cluster_labels_df$Cluster,
                   labels = paste("cluster", seq_len(n_clusters), sep="_"))
  names(labels) <- cluster_labels_df$ObsLabel

  labels
}

#' @rdname accessors_carp
#' @export
get_cluster_centroids <- function(x, ...){
  UseMethod("get_cluster_centroids")
}

#' @rdname accessors_carp
#' @export
get_cluster_centroids.CARP <- function(x, ..., percent, k, refit = TRUE){
  labels <- as.integer(get_cluster_labels(x, ..., percent = percent, k = k))

  K <- num_unique(labels)

  if(refit){
    U <- x$X
  } else {
    U <- get_U(x, ..., percent = percent, k = k)
  }

  p <- NCOL(U)

  centroids <- matrix(NA, ncol = p, nrow = K)
  colnames(centroids) <- colnames(U)

  for(k in seq_len(K)){
    centroids[k, ] <- colMeans(U[labels == k, , drop = FALSE])
  }

  centroids
}

#' @rdname accessors_carp
#' @export
get_clustered_data <- function(x, ...){
  UseMethod("get_clustered_data")
}

#' @rdname accessors_carp
#' @export
get_clustered_data.CARP <- function(x, ..., percent, k, refit = TRUE){
  labels <- as.integer(get_cluster_labels(x, ..., percent = percent, k = k))
  centroids <- get_cluster_centroids(x, ..., percent = percent, k = k, refit = refit)

  clustered_data <- x$X
  N <- NROW(clustered_data)

  for(n in seq_len(N)){
    clustered_data[n, ] <- centroids[labels[n], ]
  }

  clustered_data
}

#' @noRd
get_U <- function(x, ...){
  UseMethod("get_U")
}

#' @noRd
get_U.CARP <- function(x, ..., percent, k){
  dots <- list(...)

  if ( length(dots) != 0) {
    if (!is.null(names(dots))) {
      nm <- names(dots)
      nm <- nm[nzchar(nm)]
      crv_error("Unknown argument ", sQuote(nm[1]), " passed to ", sQuote("get_U."))
    } else {
      crv_error("Unknown argument passed to ", sQuote("get_U."))
    }
  }

  has_percent <- !missing(percent)
  has_k       <- !missing(k)
  n_args      <- has_percent + has_k

  if(n_args != 1){
    crv_error("Exactly one of ", sQuote("percent"), " and ", sQuote("k"), " must be supplied.")
  }

  if(has_k){

    if ( !is_integer_scalar(k) ){
      crv_error(sQuote("k"), " must be an integer scalar (vector of length 1).")
    }

    if( k <= 0 ) {
      crv_error(sQuote("k"), " must be positive.")
    }

    if( k > NROW(x$X) ){
      crv_error(sQuote("k"), " cannot be more than the observations in the original data set (", NROW(x$X), ").")
    }

    percent <- x$carp.cluster.path.vis %>%
      select(.data$LambdaPercent, .data$NCluster) %>%
      filter(.data$NCluster <= k) %>%
      select(.data$LambdaPercent) %>%
      summarize(percent = min(.data$LambdaPercent)) %>%
      pull
  }

  if( !is_percent_scalar(percent) ){
    crv_error(sQuote("percent"), " must be a scalar between 0 and 1 (inclusive).")
  }

  index <- which.min(abs(x$carp.sol.path$lambda.path - percent * max(x$carp.sol.path$lambda.path)))[1]

  raw_u <- matrix(x$carp.sol.path$u.path[, index],
                  nrow = x$n,
                  ncol = x$p,
                  byrow = TRUE) # byrow = TRUE because we get u by vectorizing t(X), not X

  U <- unscale_matrix(raw_u, scale = x$scale_vector, center = x$center_vector)

  ## Add rownames back in
  colnames(U) <- colnames(x$X)
  rownames(U) <- rownames(x$X)

  U
}
