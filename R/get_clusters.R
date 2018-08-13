#' Method for returning CARP and CBASS clustering solutions
#'
#' See \code{clustering.CARP} and \code{clustering.CBASS} for details
#'
#' @param x a CARP or CBASS object
#' @param ... additional arguements passed to \code{clustering.CARP} or
#' \code{clustering.CBASS}
#' @return CARP clustering solutions or CBASS biclustering solutions
#' @export
clustering <- function(x, ...) {
  UseMethod("clustering", x)
}

#' Get clustering solution from a CARP object
#'
#' Returns cluster labels and cluster means at a point along the CARP path,
#' or the entire cluster sequence of cluster labels and means
#'
#' Passing either the desired number of clusters (\code{k}) or the percent
#' regularization (\code{percent}) returns the clustering assignment
#' and cluster means at the specific point along the CARP path. If neither
#' \code{k} nor \code{percent} are specified, all clusteirng assignments and
#' mean matricies are returned.
#'
#' @param x A CARP object returned by \code{CARP}
#' @param k An interger between 1 and \code{n.obs}. The number of unique
#' clusters
#' @param percent A number between 0 and 1. The percent of regularization at
#' which to cut the path.
#' @param ... Unused additional generic arguements
#' @return A list with elements
#' \describe{
#' \item{\code{clustering.assignment}}{
#' In the case where either \code{k} or \code{percent} is specified, a vector
#' of cluster labels of length \code{n.obs}.
#' In the case where neither \code{k} nor \code{percent} is specified, a
#' matrix of size \code{n.obs} by \code{n.obs}, each row specifying a unique
#' cluster assignment along the CARP path.
#' }
#' \item{\code{cluster.means}}{
#' In the case where either \code{k} or \code{percent} is specified, a matrix
#' of dimension \code{p.var} by \code{length(unique(clustering.assignment))}
#' with each column a cluster mean.
#' In the case where neither \code{k} nor \code{percent} is specificed, a
#' list of matricies of length \code{n.obs}, representing the cluster means
#' for each cluster assignment along the CARP path.
#' }
#' }
#' @importFrom stats cutree
#' @export
#' @examples
#' library(clustRviz)
#' data("presidential_speech")
#' Xdat <- presidential_speech[1:10,1:4]
#' carp.fit <- CARP(X=Xdat)
#' # Return the CARP iterate with k=5 clusters
#' carp.clustering <- clustering(carp.fit,k=5)
#' # Examine the cluster labels
#' carp.clustering$clustering.assignment
#' # Examine the cluster means
#' head(carp.clustering$cluster.means)
#' # Return the whole sequence of solutions
#' carp.clustering.full <- clustering(carp.fit)
#' # Examine the k=5 solution again
#' carp.clustering.full$clustering.assignment[5,]
#' # Examine the k=5 means again
#' head(carp.clustering.full$cluster.means[[5]])
clustering.CARP <- function(x, k = NULL, percent = NULL, ...) {
  if (!is.null(k)) {
    clust.assign <- stats::cutree(x$carp.dend, k = k)
    lapply(unique(clust.assign), function(cl.lab) {
      apply(
        matrix(t(x$X)[, clust.assign == cl.lab], nrow = x$p.var),
        1,
        mean
      )
    }) %>%
      do.call(cbind, .) -> clust.means
    clust.assign <- paste("cl", clust.assign, sep = "")
    colnames(clust.means) <- unique(clust.assign)
  } else if (!is.null(percent)) {
    clust.assign <- stats::cutree(x$carp.dend, h = percent)
    lapply(unique(clust.assign), function(cl.lab) {
      apply(
        matrix(t(x$X)[, clust.assign == cl.lab], nrow = x$p.var),
        1,
        mean
      )
    }) %>%
      do.call(cbind, .) -> clust.means
    clust.assign <- paste("cl", clust.assign, sep = "")
    colnames(clust.means) <- unique(clust.assign)
  } else {
    lapply(1:x$n.obs, function(k) {
      stats::cutree(x$carp.dend, k)
    }) %>%
      do.call(rbind, .) -> clust.assign
    apply(clust.assign, 1, function(cl.ass) {
      lapply(unique(cl.ass), function(cl.lab) {
        apply(
          matrix(t(x$X)[, cl.ass == cl.lab], nrow = x$p.var),
          1,
          mean
        )
      }) %>%
        do.call(cbind, .) -> tmp.means
      colnames(tmp.means) <- paste("cl", 1:length(unique(cl.ass)), sep = "")
      tmp.means
    }) -> clust.means
    clust.assign <- matrix(paste("cl", clust.assign, sep = ""), nrow = x$n.obs)
  }
  list(
    clustering.assignment = clust.assign,
    cluster.means = clust.means
  )
}


#' Get biclustering solution from a CBASS object
#'
#' Returns (observation and variable) cluster labels and
#' cluster mean matrix at a point along the CBASS path.
#'
#' Passing either the desired number of clusters (\code{k}) or the percent
#' regularization (\code{percent}) returns the clustering assignment
#' and cluster mean matrix at the specific point along the CBASS path.
#'
#' @param x A CBASS object returned by \code{CBASS}
#' @param k.obs An interger between 1 and \code{n.obs}. The number of unique
#' observation clusters.
#' @param k.var An interger between 1 and \code{p.var}. The number of unique
#' variable clusters.
#' @param percent A number between 0 and 1. The percent of regularization at
#' which to cut the path.
#' @param ... Unused additional generic arguements
#' @return A list with elements
#' \describe{
#' \item{\code{clustering.assignment.obs}}{
#' A vector of observation cluster labels of length \code{n.obs}.
#' }
#' \item{\code{clustering.assignment.var}}{
#' A vector of variable cluster labels of length \code{p.var}.
#' }
#' \item{\code{cluster.mean.matrix}}{
#' The cluster mean matrix of size \code{n.obs} by \code{p.var}.
#' }
#' }
#' @importFrom stats cutree
#' @importFrom purrr map_dfr
#' @importFrom dplyr tibble
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr slice
#' @importFrom dplyr select
#' @importFrom dplyr n
#' @export
#' @examples
#' \dontrun{
#' library(clustRviz)
#' data("presidential_speech")
#' Xdat <- presidential_speech
#' cbass.fit <- CBASS(X=Xdat)
#' cbass.clustering <- clustering(cbass.fit,percent = .8)
#' }
clustering.CBASS <- function(x, k.obs = NULL, k.var = NULL, percent = NULL, ...) {
  Lambda <- NObsCl <- NVarCl <- Percent <- NULL

  n.not.null <- sum(
    c(
      !is.null(k.obs),
      !is.null(k.var),
      !is.null(percent)
    )
  )
  if (n.not.null != 1) {
    stop("Select exactly one of k.obs, k.var, or percent")
  }
  lam.vec <- x$cbass.sol.path$lambda.path %>% as.vector()
  max.lam <- max(lam.vec)
  lam.vec %>%
    purrr::map_dfr(.f = function(cur.lam) {
      # find lambda closest in column path
      cur.col.lam.ind <- which.min(abs(x$cbass.cluster.path.obs$lambda.path.inter - cur.lam))
      # find clustering solution in column path
      cur.col.clust.assignment <- x$cbass.cluster.path.obs$clust.path[[cur.col.lam.ind]]$membership
      cur.col.clust.labels <- unique(cur.col.clust.assignment)
      cur.col.nclust <- length(cur.col.clust.labels)
      # find lambda closest in rowumn path
      cur.row.lam.ind <- which.min(abs(x$cbass.cluster.path.var$lambda.path.inter - cur.lam))
      # find clustering solution in rowumn path
      cur.row.clust.assignment <- x$cbass.cluster.path.var$clust.path[[cur.row.lam.ind]]$membership
      cur.row.clust.labels <- unique(cur.row.clust.assignment)
      cur.row.nclust <- length(cur.row.clust.labels)
      dplyr::tibble(
        Lambda = cur.lam,
        NObsCl = cur.col.nclust,
        NVarCl = cur.row.nclust
      )
    }) %>%
    dplyr::mutate(
      Percent = Lambda / max.lam
    ) -> cut.table

  if (!is.null(k.obs)) {
    cut.table %>%
      dplyr::filter(NObsCl <= k.obs) %>%
      dplyr::slice(1) %>%
      dplyr::select(Lambda) %>%
      unlist() %>%
      unname() -> cur.lam
  } else if (!is.null(k.var)) {
    cut.table %>%
      dplyr::filter(NVarCl <= k.var) %>%
      dplyr::slice(1) %>%
      dplyr::select(Lambda) %>%
      unlist() %>%
      unname() -> cur.lam
  } else if (!is.null(percent)) {
    cut.table %>%
      dplyr::filter(Percent >= percent) %>%
      dplyr::slice(1) %>%
      dplyr::select(Lambda) %>%
      unlist() %>%
      unname() -> cur.lam
  } else {
    stop("Select exactly one of k.obs, k.var, or percent")
  }
  # find lambda closest in column path
  cur.col.lam.ind <- which.min(abs(x$cbass.cluster.path.obs$lambda.path.inter - cur.lam))
  # find clustering solution in column path
  cur.col.clust.assignment <- x$cbass.cluster.path.obs$clust.path[[cur.col.lam.ind]]$membership
  cur.col.clust.labels <- unique(cur.col.clust.assignment)
  cur.col.nclust <- length(cur.col.clust.labels)
  # find lambda closest in row path
  cur.row.lam.ind <- which.min(abs(x$cbass.cluster.path.var$lambda.path.inter - cur.lam))
  # find clustering solution in row path
  cur.row.clust.assignment <- x$cbass.cluster.path.var$clust.path[[cur.row.lam.ind]]$membership
  cur.row.clust.labels <- unique(cur.row.clust.assignment)
  cur.row.nclust <- length(cur.row.clust.labels)

  if (x$X.center.global) {
    X.heat <- x$X
    X.heat <- X.heat - mean(X.heat)
    X.heat <- t(X.heat)
    X <- x$X
    X <- X - mean(X)
    X <- t(X)
  } else {
    X.heat <- t(x$X)
    X <- t(x$X)
  }
  colnames(X.heat) <- x$obs.labels
  rownames(X.heat) <- x$var.labels
  for (col.label.ind in seq_along(cur.col.clust.labels)) {
    cur.col.label <- cur.col.clust.labels[col.label.ind]
    col.inds <- which(cur.col.clust.assignment == cur.col.label)
    for (row.label.ind in seq_along(cur.row.clust.labels)) {
      cur.row.label <- cur.row.clust.labels[row.label.ind]
      row.inds <- which(cur.row.clust.assignment == cur.row.label)
      mean.value <- mean(X[row.inds, col.inds])
      X.heat[row.inds, col.inds] <- mean.value
    }
  }


  clust.assign.obs <- paste("cl", cur.col.clust.assignment, sep = "")
  clust.assign.var <- paste("cl", cur.row.clust.assignment, sep = "")
  list(
    clustering.assignment.obs = clust.assign.obs,
    clustering.assignment.var = clust.assign.var,
    cluster.mean.matrix = X.heat
  )
}
