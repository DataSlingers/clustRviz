## Accessor functions for CBASS

#' Get Clustering Results for \code{CBASS}
#'
#' \code{get_cluster_labels} returns a factor vector of cluster labels.
#' \code{get_cluster_centroids} returns a \code{k1}-by-\code{k2} matrix with the
#' estimated centroid of the \code{k1}-th observation cluster and the \code{k2}-th
#' feature cluster.
#' \code{get_clustered_data} returns a matrix (with the same dimensions and names
#' as the original data), but with the values for each observation replaced by
#' its "estimated" value (\emph{i.e.}, the appropriate cluster centroid).
#'
#' @param x An object of class \code{CARP} as produced by \code{\link{CBASS}}
#' @param percent A number between 0 and 1, giving the regularization level (as
#'                a fraction of the final regularization level used) at which to
#'                get cluster labels.
#' @param k.obs The desired number of observation clusters
#' @param k.var The desired number of variable clusters
#' @param type For \code{get_cluster_labels}, which set of labels to return -
#'             observation (row) or feature (column)
#' @param refit Should "naive" centroids (\code{TRUE}) or the actual centroids
#'              estimated by convex clustering be used? The default (\code{refit = TRUE})
#'              centroids returned are actual centroids (mean) of all elements
#'              assigned to that cluster; if \code{refit = FALSE}, the \eqn{\hat{U}}
#'              from the convex biclustering problem is used. Due to the global
#'              shrinkage imposed, these clusters are more "shrunk together" than
#'              the naive clusters.
#' @param ... Additional arguments - if any are provided, an error is signalled.
#' @details \code{get_clustered_data} returns centroids on the original scale of
#' the data, independent of any pre-processing flags passed to \code{CBASS}.
#' Note that exactly one of \code{percent}, \code{k.obs}, \code{k.var}
#' must be supplied and that that \code{k.obs} (if suppplied) will be
#' used even if \code{type = "var"} and \emph{vice versa}.
#'
#' @examples
#' cbass_fit <- CBASS(presidential_speech)
#'
#' # Get observation clustering results from 50% along the path
#' get_cluster_labels(cbass_fit, percent = 0.5)
#'
#' # Get variable clustering corresponding to the 3 cluster solution
#' get_cluster_labels(cbass_fit, k.var = 3, type = "var")
#'
#' # Get observation clustering corresponding to the 3 variable clusters
#' get_cluster_labels(cbass_fit, k.var = 3, type = "obs")
#'
#' # Get cluster centroids partially down the path
#' get_cluster_centroids(cbass_fit, percent = 0.5)
#'
#' # Get clustered data
#' image(get_clustered_data(cbass_fit, k.obs = 2))
#' @export
#' @rdname accessors_cbass
#' @importFrom dplyr select filter summarize pull arrange
#' @importFrom rlang .data
get_cluster_labels.CBASS <- function(x, ..., percent, k.obs, k.var, type = c("obs", "var")){

  dots <- list(...)
  if ( length(dots) != 0 ) {
    if (!is.null(names(dots))) {
      nm <- names(dots)
      nm <- nm[nzchar(nm)]
      stop("Unknown argument ", sQuote(nm[1]), " passed to ", sQuote("get_cluster_labels."))
    } else {
      stop("Unknown argument passed to ", sQuote("get_cluster_labels."))
    }
  }

  type <- match.arg(type)

  has_percent <- !missing(percent)
  has_k.obs   <- !missing(k.obs)
  has_k.var   <- !missing(k.var)
  n_args      <- has_percent + has_k.obs + has_k.var

  if(n_args != 1){
    stop("Exactly one of ", sQuote("percent,"), " ", sQuote("k.obs"),
         " and ", sQuote("k.var"), " must be supplied.")
  }

  if(has_k.obs){

    if ( !is_integer_scalar(k.obs) ){
      stop(sQuote("k"), " must be an integer scalar (vector of length 1).")
    }

    if( k.obs <= 0 ) {
      stop(sQuote("k.obs"), " must be positive.")
    }

    if( k.obs > NROW(x$X) ){
      stop(sQuote("k.obs"), " cannot be more than the observations in the original data set (", NROW(x$X), ").")
    }

    percent <- x$cbass.cluster.path.vis.obs %>%
      select(.data$LambdaPercent, .data$NCluster) %>%
      filter(.data$NCluster <= k.obs) %>%
      select(.data$LambdaPercent) %>%
      summarize(percent = min(.data$LambdaPercent)) %>%
      pull
  }

  if(has_k.var){

    if ( !is_integer_scalar(k.var) ){
      stop(sQuote("k"), " must be an integer scalar (vector of length 1).")
    }

    if( k.var <= 0 ) {
      stop(sQuote("k.var"), " must be positive.")
    }

    if( k.var > NCOL(x$X) ){
      stop(sQuote("k.var"), " cannot be more than the features in the original data set (", NCOL(x$X), ").")
    }

    percent <- x$cbass.cluster.path.vis.var %>%
      select(.data$LambdaPercent, .data$NCluster) %>%
      filter(.data$NCluster <= k.var) %>%
      select(.data$LambdaPercent) %>%
      summarize(percent = min(.data$LambdaPercent)) %>%
      pull
  }

  if( !is_percent_scalar(percent) ){
    stop(sQuote("percent"), " must be a scalar between 0 and 1 (inclusive).")
  }

  if(type == "obs"){
    cluster_labels_df <- x$cbass.cluster.path.vis.obs %>%
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
  } else {
    cluster_labels_df <- x$cbass.cluster.path.vis.var %>%
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
  }

  labels
}

#' @export
#' @rdname accessors_cbass
get_cluster_centroids.CBASS <- function(x, ..., percent, k.var, k.obs, refit = TRUE){
  obs_labels <- as.integer(get_cluster_labels(x, ...,
                                              percent = percent,
                                              k.var = k.var,
                                              k.obs = k.obs,
                                              type  = "obs"))

  var_labels <- as.integer(get_cluster_labels(x, ...,
                                              percent = percent,
                                              k.var = k.var,
                                              k.obs = k.obs,
                                              type  = "var"))

  centroids <- matrix(NA, nrow = num_unique(obs_labels), ncol = num_unique(var_labels))

  if(refit){
    U <- x$X
  } else {
    U <- get_U(x, ..., percent = percent, k.var = k.var, k.obs = k.obs)
  }

  for(o in unique(obs_labels)){
    for(v in unique(var_labels)){
      centroids[o, v] <- mean(U[obs_labels == o, var_labels == v])
    }
  }

  centroids
}

#' @export
#' @rdname accessors_cbass
get_clustered_data.CBASS <- function(x, ..., percent, k.var, k.obs, refit = TRUE){
  obs_labels <- as.integer(get_cluster_labels(x, ...,
                                              percent = percent,
                                              k.var = k.var,
                                              k.obs = k.obs,
                                              type  = "obs"))

  var_labels <- as.integer(get_cluster_labels(x, ...,
                                              percent = percent,
                                              k.var = k.var,
                                              k.obs = k.obs,
                                              type  = "var"))

  centroids <- get_cluster_centroids(x, ...,
                                     percent = percent,
                                     k.var = k.var,
                                     k.obs = k.obs,
                                     refit = refit)

  X <- x$X
  clustered_data <- X * NA

  for(o in unique(obs_labels)){
    for(v in unique(var_labels)){
      clustered_data[obs_labels == o, var_labels == v] <- centroids[o, v]
    }
  }

  clustered_data
}

#' @noRd
get_U.CBASS <- function(x, ..., percent, k.var, k.obs){
  dots <- list(...)

  if ( length(dots) != 0) {
    if (!is.null(names(dots))) {
      nm <- names(dots)
      nm <- nm[nzchar(nm)]
      stop("Unknown argument ", sQuote(nm[1]), " passed to ", sQuote("get_U."))
    } else {
      stop("Unknown argument passed to ", sQuote("get_U."))
    }
  }

  has_percent <- !missing(percent)
  has_k.obs   <- !missing(k.obs)
  has_k.var   <- !missing(k.var)
  n_args      <- has_percent + has_k.obs + has_k.var

  if(n_args != 1){
    stop("Exactly one of ", sQuote("percent,"), " ", sQuote("k.obs"),
         " and ", sQuote("k.var"), " must be supplied.")
  }

  if(has_k.obs){

    if ( !is_integer_scalar(k.obs) ){
      stop(sQuote("k"), " must be an integer scalar (vector of length 1).")
    }

    if( k.obs <= 0 ) {
      stop(sQuote("k.obs"), " must be positive.")
    }

    if( k.obs > NROW(x$X) ){
      stop(sQuote("k.obs"), " cannot be more than the observations in the original data set (", NROW(x$X), ").")
    }

    percent <- x$cbass.cluster.path.vis.obs %>%
      select(.data$LambdaPercent, .data$NCluster) %>%
      filter(.data$NCluster <= k.obs) %>%
      select(.data$LambdaPercent) %>%
      summarize(percent = min(.data$LambdaPercent)) %>%
      pull
  }

  if(has_k.var){

    if ( !is_integer_scalar(k.var) ){
      stop(sQuote("k"), " must be an integer scalar (vector of length 1).")
    }

    if( k.var <= 0 ) {
      stop(sQuote("k.var"), " must be positive.")
    }

    if( k.var > NCOL(x$X) ){
      stop(sQuote("k.var"), " cannot be more than the features in the original data set (", NCOL(x$X), ").")
    }

    percent <- x$cbass.cluster.path.vis.var %>%
      select(.data$LambdaPercent, .data$NCluster) %>%
      filter(.data$NCluster <= k.var) %>%
      select(.data$LambdaPercent) %>%
      summarize(percent = min(.data$LambdaPercent)) %>%
      pull
  }

  if( !is_percent_scalar(percent) ){
    stop(sQuote("percent"), " must be a scalar between 0 and 1 (inclusive).")
  }

  index <- which.min(abs(x$cbass.sol.path$lambda.path - percent * max(x$cbass.sol.path$lambda.path)))[1]

  raw_u <- matrix(x$cbass.sol.path$u.path[, index],
                  nrow = x$n.obs,
                  ncol = x$p.var,
                  byrow = TRUE) # byrow = TRUE because we get u by vectorizing t(X), not X

  U <- raw_u + x$mean_adjust

  ## Add rownames back in
  colnames(U) <- colnames(x$X)
  rownames(U) <- rownames(x$X)

  U
}
