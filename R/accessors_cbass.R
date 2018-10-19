## Accessor functions for CBASS

#' Get Clustering Results for \code{CBASS}
#'
#' \code{get_cluster_labels} returns a factor vector of cluster labels.
#' \code{get_cluster_centroids} returns a \code{k1}-by-\code{k2} matrix with the
#' estimated centroid of the \code{k1}-th row cluster and the \code{k2}-th
#' column cluster.
#' \code{get_clustered_data} returns a matrix (with the same dimensions and names
#' as the original data), but with the values for each row replaced by
#' its "estimated" value (\emph{i.e.}, the appropriate cluster centroid).
#'
#' @param x An object of class \code{CARP} as produced by \code{\link{CBASS}}
#' @param percent A number between 0 and 1, giving the regularization level (as
#'                a fraction of the final regularization level used) at which to
#'                get cluster labels.
#' @param k.row The desired number of row clusters
#' @param k.col The desired number of column clusters
#' @param type For \code{get_cluster_labels}, which set of labels to return -
#'             row (observation) or column (feature)
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
#' Note that exactly one of \code{percent}, \code{k.row}, \code{k.col}
#' must be supplied and that that \code{k.row} (if suppplied) will be
#' used even if \code{type = "col"} and \emph{vice versa}.
#'
#' @examples
#' cbass_fit <- CBASS(presidential_speech)
#'
#' # Get row clustering results from 50% along the path
#' get_cluster_labels(cbass_fit, percent = 0.5)
#'
#' # Get column clustering corresponding to the 3 cluster solution
#' get_cluster_labels(cbass_fit, k.col = 3, type = "col")
#'
#' # Get row clustering corresponding to the 3 column clusters
#' get_cluster_labels(cbass_fit, k.col = 3, type = "row")
#'
#' # Get cluster centroids partially down the path
#' get_cluster_centroids(cbass_fit, percent = 0.5)
#'
#' # Get clustered data
#' image(get_clustered_data(cbass_fit, k.row = 2))
#' @export
#' @rdname accessors_cbass
#' @importFrom dplyr select filter summarize pull arrange
#' @importFrom rlang .data
get_cluster_labels.CBASS <- function(x, ..., percent, k.row, k.col, type = c("row", "col")){

  dots <- list(...)
  if ( length(dots) != 0 ) {
    if (!is.null(names(dots))) {
      nm <- names(dots)
      nm <- nm[nzchar(nm)]
      crv_error("Unknown argument ", sQuote(nm[1]), " passed to ", sQuote("get_cluster_labels."))
    } else {
      crv_error("Unknown argument passed to ", sQuote("get_cluster_labels."))
    }
  }

  type <- match.arg(type)

  has_percent <- !missing(percent)
  has_k.row   <- !missing(k.row)
  has_k.col   <- !missing(k.col)
  n_args      <- has_percent + has_k.row + has_k.col

  if(n_args != 1){
    crv_error("Exactly one of ", sQuote("percent,"), " ", sQuote("k.row"),
             " and ", sQuote("k.col"), " must be supplied.")
  }

  if(has_k.row){

    if ( !is_positive_integer_scalar(k.row) ){
      crv_error(sQuote("k"), " must be a positive integer scalar (vector of length 1).")
    }

    if( k.row > NROW(x$X) ){
      crv_error(sQuote("k.row"), " cannot be more than the rows in the original data set (", NROW(x$X), ").")
    }

    percent <- x$row_fusions$cluster_membership %>%
      select(.data$GammaPercent, .data$NCluster) %>%
      filter(.data$NCluster <= k.row) %>%
      select(.data$GammaPercent) %>%
      summarize(percent = min(.data$GammaPercent)) %>%
      pull
  }

  if(has_k.col){

    if ( !is_positive_integer_scalar(k.col) ){
      crv_error(sQuote("k"), " must be a positive integer scalar (vector of length 1).")
    }

    if( k.col > NCOL(x$X) ){
      crv_error(sQuote("k.col"), " cannot be more than the columns in the original data set (", NCOL(x$X), ").")
    }

    percent <- x$col_fusions$cluster_membership %>%
      select(.data$GammaPercent, .data$NCluster) %>%
      filter(.data$NCluster <= k.col) %>%
      select(.data$GammaPercent) %>%
      summarize(percent = min(.data$GammaPercent)) %>%
      pull
  }

  if( !is_percent_scalar(percent) ){
    crv_error(sQuote("percent"), " must be a scalar between 0 and 1 (inclusive).")
  }

  if(type == "row"){
    cluster_labels_df <- x$row_fusions$cluster_membership %>%
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
  } else {
    cluster_labels_df <- x$col_fusions$cluster_membership %>%
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
  }

  labels
}

#' @export
#' @rdname accessors_cbass
get_cluster_centroids.CBASS <- function(x, ..., percent, k.row, k.col, refit = TRUE){
  row_labels <- as.integer(get_cluster_labels(x, ...,
                                              percent = percent,
                                              k.col = k.col,
                                              k.row = k.row,
                                              type  = "row"))

  col_labels <- as.integer(get_cluster_labels(x, ...,
                                              percent = percent,
                                              k.col = k.col,
                                              k.row = k.row,
                                              type  = "col"))

  centroids <- matrix(NA, nrow = num_unique(row_labels), ncol = num_unique(col_labels))

  if(refit){
    U <- x$X
  } else {
    U <- get_U(x, ..., percent = percent, k.col = k.col, k.row = k.row)
  }

  for(r in unique(row_labels)){
    for(c in unique(col_labels)){
      centroids[r, c] <- mean(U[row_labels == r, col_labels == c])
    }
  }

  centroids
}

#' @export
#' @rdname accessors_cbass
get_clustered_data.CBASS <- function(x, ..., percent, k.row, k.col, refit = TRUE){
  row_labels <- as.integer(get_cluster_labels(x, ...,
                                              percent = percent,
                                              k.row = k.row,
                                              k.col = k.col,
                                              type  = "row"))

  col_labels <- as.integer(get_cluster_labels(x, ...,
                                              percent = percent,
                                              k.row = k.row,
                                              k.col = k.col,
                                              type  = "col"))

  centroids <- get_cluster_centroids(x, ...,
                                     percent = percent,
                                     k.row = k.row,
                                     k.col = k.col,
                                     refit = refit)

  X <- x$X
  clustered_data <- X * NA

  for(r in unique(row_labels)){
    for(c in unique(col_labels)){
      clustered_data[row_labels == r, col_labels == c] <- centroids[r, c]
    }
  }

  clustered_data
}

#' @noRd
get_U.CBASS <- function(x, ..., percent, k.row, k.col){
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
  has_k.row   <- !missing(k.row)
  has_k.col   <- !missing(k.col)
  n_args      <- has_percent + has_k.row + has_k.col

  if(n_args != 1){
    crv_error("Exactly one of ", sQuote("percent,"), " ", sQuote("k.row"),
              " and ", sQuote("k.col"), " must be supplied.")
  }

  if(has_k.row){

    if ( !is_positive_integer_scalar(k.row) ){
      crv_error(sQuote("k"), " must be a positive integer scalar (vector of length 1).")
    }

    if( k.row > NROW(x$X) ){
      crv_error(sQuote("k.row"), " cannot be more than the rows in the original data set (", NROW(x$X), ").")
    }

    percent <- x$row_fusions$cluster_membership %>%
      select(.data$GammaPercent, .data$NCluster) %>%
      filter(.data$NCluster <= k.row) %>%
      select(.data$GammaPercent) %>%
      summarize(percent = min(.data$GammaPercent)) %>%
      pull
  }

  if(has_k.col){

    if ( !is_positive_integer_scalar(k.col) ){
      crv_error(sQuote("k"), " must be a positive integer scalar (vector of length 1).")
    }

    if( k.col > NCOL(x$X) ){
      crv_error(sQuote("k.col"), " cannot be more than the columns in the original data set (", NCOL(x$X), ").")
    }

    percent <- x$col_fusions$cluster_membership %>%
      select(.data$GammaPercent, .data$NCluster) %>%
      filter(.data$NCluster <= k.col) %>%
      select(.data$GammaPercent) %>%
      summarize(percent = min(.data$GammaPercent)) %>%
      pull
  }

  if( !is_percent_scalar(percent) ){
    crv_error(sQuote("percent"), " must be a scalar between 0 and 1 (inclusive).")
  }

  ## Pull out the iter for the closest value of "GammaPercent" to the desired percent
  ## slice(which.min(...)[1]) will pull the "which.min(...)[1]"-th element
  index <- x$row_fusions$cluster_membership %>%
              slice(which.min(abs(.data$GammaPercent - percent))[1]) %>% pull(.data$Iter)
  raw_u <- x$row_fusions$U[, , index] ## The choice of row_fusions here is a bit arbitrary
                                      ## It would be better to improve post-processing to
                                      ## do simultaneous interpolation for rows and columns

  U <- raw_u + x$mean_adjust

  ## Add rownames back in
  colnames(U) <- colnames(x$X)
  rownames(U) <- rownames(x$X)

  U
}


is_raw_feature.CBASS <- function(x, f, type = c("row", "col"), ...){
  type <- match.arg(type)

  ## This is a bit confusing, but if we are clustering rows, our "feature" are
  ## columns (the usual case), while if we are clustering columns, our features are now
  ## the rows
  if (type == "row") {
    f %in% colnames(x$X)
  } else {
    f %in% rownames(x$X)
  }
}

is_pc_feature.CBASS <- function(x, f, type = c("row", "col"), ...){
  type <- match.arg(type)

  ## First check if the feature name is of the form "PC###" or ".PC###"
  ## If so, check that the implied PC is less than the number of singular vectors we kept

  if (type == "row") {
    (grepl(pattern = "[.]?PC[0123456789]+", f)) && (as.integer(gsub("[^0123456789]", "", f)) <= NCOL(x$row_fusions$rotation_matrix))
  } else {
    (grepl(pattern = "[.]?PC[0123456789]+", f)) && (as.integer(gsub("[^0123456789]", "", f)) <= NCOL(x$col_fusions$rotation_matrix))
  }
}

get_pc_path.CBASS <- function(x, f, type = c("row", "col"), ...){
  type <- match.arg(type)
  pc_num <- as.integer(gsub("[^0123456789]", "", f))

  if (type == "row") {
    as.vector(tensor_projection(x$row_fusions$U, x$row_fusions$rotation_matrix[, pc_num, drop = FALSE]))
  } else {
    as.vector(tensor_projection(x$col_fusions$U, x$col_fusions$rotation_matrix[, pc_num, drop = FALSE]))
  }
}

get_feature_paths.CBASS <- function(x, features, type = c("row", "col"), ...){

  dots <- list(...)
  if (length(dots)) {
    crv_error("Unknown arguments passed to", sQuote("get_feature_paths.CBASS."))
  }

  type <- match.arg(type)

  path_info <- x[[if(type == "row") "row_fusions" else "col_fusions"]]$cluster_membership

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
    if (is_raw_feature(x, f, type = type)) {
      ## Get the path for `f` and add it to `path_info`
      if(type == "row"){
        path_info[[f]] <- as.vector(x$row_fusions$U[,f,])
      } else {
        path_info[[f]] <- as.vector(x$col_fusions$U[,f,])
      }
    } else if (is_pc_feature(x, f)) {
      ## Get the path for `f` and add it to `path_info`
      path_info[[f]] <- get_pc_path(x, f, type = type)
    } else {
      crv_error(sQuote(f), " is not an original feature or principal component.")
    }
  }
  path_info
}
