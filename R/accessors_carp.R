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

    percent <- x$cluster_membership %>%
                 select(.data$GammaPercent, .data$NCluster) %>%
                 filter(.data$NCluster <= k) %>%
                 select(.data$GammaPercent) %>%
                 summarize(percent = min(.data$GammaPercent)) %>%
                 pull
  }

  if( !is_percent_scalar(percent) ){
    crv_error(sQuote("percent"), " must be a scalar between 0 and 1 (inclusive).")
  }

  cluster_labels_df <- x$cluster_membership %>%
                         select(.data$GammaPercent,
                                .data$ObsLabel,
                                .data$Obs,
                                .data$Cluster,
                                .data$Iter) %>%
                         filter(.data$GammaPercent >= percent) %>%
                         filter(.data$GammaPercent == min(.data$GammaPercent)) %>%
                         filter(.data$Iter == min(.data$Iter)) %>% # In case we have multiple iterations at same gamma
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

    if ( !is_positive_integer_scalar(k) ){
      crv_error(sQuote("k"), " must be a positive integer scalar (vector of length 1).")
    }

    if( k > NROW(x$X) ){
      crv_error(sQuote("k"), " cannot be more than the observations in the original data set (", NROW(x$X), ").")
    }

    percent <- x$cluster_membership %>%
      select(.data$GammaPercent, .data$NCluster) %>%
      filter(.data$NCluster <= k) %>%
      select(.data$GammaPercent) %>%
      summarize(percent = min(.data$GammaPercent)) %>%
      pull
  }

  if( !is_percent_scalar(percent) ){
    crv_error(sQuote("percent"), " must be a scalar between 0 and 1 (inclusive).")
  }

  gamma_path <- x$cluster_membership %>% select(.data$Iter, .data$GammaPercent) %>%
                                         distinct

  ## Pull out the iter for the closest value of "GammaPercent" to the desired percent
  ## slice(which.min(...)[1]) will pull the "which.min(...)[1]"-th element
  index <- gamma_path %>% slice(which.min(abs(.data$GammaPercent - percent))[1]) %>% pull(.data$Iter)
  raw_u <- x$U[, , index]

  U <- unscale_matrix(raw_u, scale = x$scale_vector, center = x$center_vector)

  ## Add rownames back in
  colnames(U) <- colnames(x$X)
  rownames(U) <- rownames(x$X)

  U
}

#' @noRd
available_features <- function(x, ...){
  UseMethod("available_features")
}

available_features.CARP <- function(x, ...){
  c(paste0("PC", seq_len(NCOL(x$rotation_matrix))), colnames(x$X))
}

#' @noRd
is_raw_feature <- function(x, f, ...){
  UseMethod("is_raw_feature")
}

is_raw_feature.CARP <- function(x, f, ...){
  f %in% colnames(x$X)
}


#' @noRd
is_pc_feature <- function(x, f, ...){
  UseMethod("is_pc_feature")
}

is_pc_feature.CARP <- function(x, f, ...){
  ## First check if the feature name is of the form "PC###" or ".PC###"
  ## If so, check that the implied PC is less than the number of singular vectors we kept
  (grepl(pattern = "[.]?PC[0123456789]+", f)) && (as.integer(gsub("[^0123456789]", "", f)) <= NCOL(x$rotation_matrix))
}

#' @noRd
get_pc_path <- function(x, f, ...){
  UseMethod("get_pc_path")
}

get_pc_path.CARP <- function(x, f, ...){
  pc_num <- as.integer(gsub("[^0123456789]", "", f))

  as.vector(tensor_projection(x$U, x$rotation_matrix[, pc_num, drop = FALSE]))
}

#' @noRd
get_feature_paths <- function(x, features, ...){
  UseMethod("get_feature_paths")
}

#' @noRd
get_feature_paths.CARP <- function(x, features, ...){
  dots <- list(...)
  if (length(dots)) {
    crv_error("Unknown arguments passed to", sQuote("get_feature_paths.CARP."))
  }

  path_info <- x$cluster_membership

  ## Avoid duplicates
  if (anyDuplicated(features)) {
    crv_warning("Some features requested multiple times - omitting duplicates.")
    features <- unique(features)
  }

  for(f in features){
    ## Check that `f` is a valid feature
    if (!is_nonempty_character_scalar(f)) {
      crv_error(sQuote(f), " is not a valid feature name.")
    }

    ## Find f
    if (is_raw_feature(x, f)) {
      ## Get the path for `f` and add it to `path_info`
      path_info[[f]] <- as.vector(x$U[,f,])
    } else if (is_pc_feature(x, f)) {
      ## Get the path for `f` and add it to `path_info`
      path_info[[f]] <- get_pc_path(x, f)
    } else {
      crv_error(sQuote(f), " is not an original feature or principal component.")
    }
  }

  path_info
}
