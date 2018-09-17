## Accessor functions for CBASS

#' Get Clustering Results for \code{CBASS}
#'
#' \code{get_cluster_labels} returns a factor vector of cluster labels.
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
#' @param ... Additional arguments - if any are provided, an error is signalled.
#' @details \code{get_clustered_data} returns centroids on the original scale of
#' the data, independent of any pre-processing flags passed to \code{CBASS}.
#' These centroids are "naive" in the sense that they are means of the data points
#' in that cluster and are not shrunk. That is to say, they are not
#' calculated using the solution of the convex biclustering problem.
#'
#' Note that exactly one of \code{percent}, \code{k.obs}, \code{k.var}
#' must be supplied and that that \code{k.obs} (if suppplied) will be
#' used even if \code{type = "var"} and \emph{vice versa}.
#'
#' @examples
#' cbass_fit <- CBASS(presidential_speech)
#'
#' # Get observation clustering results from 50% along the path
#' get_cluster_labels(carp_fit, percent = 0.5)
#'
#' # Get variable clustering corresponding to the 3 cluster solution
#' get_cluster_labels(cbass_fit, k.var = 3, type = "var")
#'
#' # Get observation clustering corresponding to the 3 variable clusters
#' get_cluster_labels(cbass_fit, k.var = 3, type = "obs")
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
             .data$Cluster) %>%
      filter(.data$LambdaPercent >= percent) %>%
      filter(.data$LambdaPercent == min(.data$LambdaPercent)) %>%
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
             .data$Cluster) %>%
      filter(.data$LambdaPercent >= percent) %>%
      filter(.data$LambdaPercent == min(.data$LambdaPercent)) %>%
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
get_cluster_centroids.CBASS <- function(x, ..., percent, k.var, k.obs)

#' @export
#' @rdname accessors_cbass
get_clustered_data.CBASS <- function(x, ..., percent, k.var, k.obs){
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

  X <- x$X
  clustered_data <- X * NA

  for(o in unique(obs_labels)){
    for(v in unique(var_labels)){
      clustered_data[obs_labels == o, var_labels == v] <- mean(X[obs_labels == o, var_labels == v])
    }
  }

  clustered_data
}


