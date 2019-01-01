#' Compute \code{CBASS} (Convex BiClustering) Solution Path
#'
#' \code{CBASS} returns a fast approximation to the Convex BiClustering
#' solution path along with visualizations such as dendrograms and
#' heatmaps. \code{CBASS} solves the Convex Biclustering problem using an efficient
#' Algorithmic Regularization scheme.
#'
#' @param X The data matrix (\eqn{X \in R^{n \times p}}{X})
#' @param row_weights One of the following: \itemize{
#'                    \item A function which, when called with argument \code{X},
#'                          returns a n-by-n matrix of fusion weights.
#'                    \item A matrix of size n-by-n containing fusion weights
#'                    }
#'                    Note that the weights will be renormalized to sum to
#'                    \eqn{1/\sqrt{n}} internally.
#' @param col_weights One of the following: \itemize{
#'                    \item A function which, when called with argument \code{t(X)},
#'                          returns a p-by-p matrix of fusion weights. (Note the
#'                          transpose.)
#'                    \item A matrix of size p-by-p containing fusion weights
#'                    }
#'                    Note that the weights will be renormalized to sum to
#'                    \eqn{1/\sqrt{p}} internally.
#' @param row_labels A character vector of length \eqn{n}: row (observation) labels
#' @param col_labels A character vector of length \eqn{p}: column (variable) labels
#' @param X.center.global A logical: Should \code{X} be centered globally?
#'                        \emph{I.e.}, should the global mean of \code{X} be subtracted?
#' @param back_track A logical: Should back-tracking be used to exactly identify fusions?
#'                   By default, back-tracking is not used.
#' @param exact A logical: Should the exact solution be computed using an iterative algorithm?
#'              By default, algorithmic regularization is applied and the exact solution
#'              is not computed. Setting \code{exact = TRUE} often significantly increases
#'              computation time.
#' @param norm Which norm to use in the fusion penalty? Currently only \code{1}
#'             and \code{2} (default) are supported.
#' @param t A number greater than 1: the size of the multiplicative update to
#'          the cluster fusion regularization parameter (not used by
#'          back-tracking variants). Typically on the scale of \code{1.005} to \code{1.1}.
#' @param npcs An integer >= 2. The number of principal components to compute
#'             for path visualization.
#' @param dendrogram.scale A character string denoting how the scale of dendrogram
#'                         regularization proportions should be visualized.
#'                         Choices are \code{'original'} or \code{'log'}; if not
#'                         provided, a data-driven heuristic choice is used.
#' @param ... Unused arguements. An error will be thrown if any unrecognized
#'            arguments as given.
#' @param status Should a status message be printed to the console?
#' @return An object of class \code{CBASS} containing the following elements (among others):
#'         \itemize{
#'         \item \code{X}: the original data matrix
#'         \item \code{n}: the number of observations (rows of \code{X})
#'         \item \code{p}: the number of variables (columns of \code{X})
#'         \item \code{alg.type}: the \code{CBASS} variant used
#'         \item \code{row_fusions}: A record of row fusions - see the documentation
#'                                   of \code{\link{CARP}} for details of what this
#'                                   may include.
#'         \item \code{col_fusions}: A record of column fusions - see the documentation
#'                                   of \code{\link{CARP}} for details of what this
#'                                   may include.
#'         }
#' @export
#' @examples
#' \dontrun{
#' cbass_fit <- CBASS(presidential_speech)
#' print(cbass_fit)
#' plot(cbass_fit)
#' }
CBASS <- function(X,
                  ...,
                  row_weights = sparse_rbf_kernel_weights(k = "auto",
                                                          phi = "auto",
                                                          dist.method = "euclidean",
                                                          p = 2),
                  col_weights = sparse_rbf_kernel_weights(k = "auto",
                                                          phi = "auto",
                                                          dist.method = "euclidean",
                                                          p = 2),
                  row_labels = rownames(X),
                  col_labels = colnames(X),
                  X.center.global = TRUE,
                  t = 1.01,
                  back_track = FALSE,
                  exact = FALSE,
                  norm = 2,
                  npcs = min(4L, NCOL(X), NROW(X)),
                  dendrogram.scale = NULL,
                  status = (interactive() && (clustRviz_logger_level() %in% c("MESSAGE", "WARNING", "ERROR")))) {

  tic <- Sys.time()

  ####################
  ##
  ## Input validation
  ##
  ####################

  dots <- list(...)

  if (length(dots) != 0L) {
    if (!is.null(names(dots))) {
      crv_error("Unknown argument ", sQuote(names(dots)[1L]), " passed to ", sQuote("CBASS."))
    } else {
      crv_error("Unknown ", sQuote("..."), " arguments passed to ", sQuote("CBASS."))
    }
  }

  if (!is.matrix(X)) {
    crv_warning(sQuote("X"), " should be a matrix, not a " , class(X)[1],
                ". Converting with as.matrix().")
    X <- as.matrix(X)
  }

  if (!is.numeric(X)) {
    crv_error(sQuote("X"), " must be numeric.")
  }

  if (anyNA(X)) {
    crv_error(sQuote("CBASS"), " cannot handle missing data.")
  }

  if (!all(is.finite(X))) {
    crv_error("All elements of ", sQuote("X"), " must be finite.")
  }

  if (!is.logical(X.center.global) || is.na(X.center.global) || (length(X.center.global) != 1L)) {
    crv_error(sQuote("X.center.global"), "must be either ", sQuote("TRUE"), " or ", sQuote("FALSE."))
  }

  if (!is.null(dendrogram.scale)) {
    if (dendrogram.scale %not.in% c("original", "log")) {
      crv_error("If not NULL, ", sQuote("dendrogram.scale"), " must be either ", sQuote("original"), " or ", sQuote("log."))
    }
  }

  if ( (!is_integer_scalar(npcs)) || (npcs < 2) || (npcs > NCOL(X)) || (npcs > NROW(X)) ){
    crv_error(sQuote("npcs"), " must be an integer scalar between 2 and ", sQuote("min(dim(X))."))
  }

  if (!is_logical_scalar(back_track)) {
    crv_error(sQuote("back_track"), "must be either ", sQuote("TRUE"), " or ", sQuote("FALSE."))
  }

  if (!is_logical_scalar(exact)) {
    crv_error(sQuote("exact"), "must be either ", sQuote("TRUE"), " or ", sQuote("FALSE."))
  }

  if (norm %not.in% c(1, 2)){
    crv_error(sQuote("norm"), " must be either 1 or 2.")
  }

  l1 <- (norm == 1)

  if ( (!is_numeric_scalar(t)) || (t <= 1) ) {
    crv_error(sQuote("t"), " must be a scalar greater than 1.")
  }

  ## Get row (observation) labels
  if (is.null(row_labels)) {
    row_labels <- paste0("Row", seq_len(NROW(X)))
  }

  if ( length(row_labels) != NROW(X) ){
    crv_error(sQuote("row_labels"), " must be of length ", sQuote("NROW(X)."))
  }

  rownames(X) <- row_labels <- make.unique(as.character(row_labels), sep = "_")

  ## Get column (variable) labels
  if (is.null(col_labels)) {
    col_labels <- paste0("Col", seq_len(NCOL(X)))
  }

  if ( length(col_labels) != NCOL(X) ){
    crv_error(sQuote("col_labels"), " must be of length ", sQuote("NCOL(X)."))
  }

  colnames(X) <- col_labels <- make.unique(as.character(col_labels), sep = "_")

  n <- NROW(X)
  p <- NCOL(X)

  # Preprocess X
  X.orig <- X
  if (X.center.global) {
    mean_adjust <- mean(X)
    X <- X - mean_adjust
  } else {
    mean_adjust <- 0
  }

  crv_message("Pre-computing column weights and edge sets")
  # Calculate column (variable/feature)-clustering weights
  if (is.function(col_weights)) { # Usual case, `col_weights` is a function which calculates the weight matrix
    col_weight_result <- col_weights(t(X))

    if (is.matrix(col_weight_result)) {
      col_weight_matrix <- col_weight_result
      col_weight_type   <- UserFunction()
    } else {
      col_weight_matrix <- col_weight_result$weight_mat
      col_weight_type   <- col_weight_result$type
    }
  } else if (is.matrix(col_weights)) {

    if (!is_square(col_weights)) {
      crv_error(sQuote("col_weights"), " must be a square matrix.")
    }

    if (NROW(col_weights) != NCOL(X)) {
      crv_error(sQuote("NROW(col_weights)"), " must be equal to ", sQuote("NCOL(X)."))
    }

    col_weight_matrix <- col_weights
    col_weight_type   <- UserMatrix()
  } else {
    crv_error(sQuote("CBASS"), " does not know how to handle ", sQuote("col_weights"),
              " of class ", class(col_weights)[1], ".")
  }

  if (any(col_weight_matrix < 0) || anyNA(col_weight_matrix)) {
    crv_error("All column fusion weights must be positive or zero.")
  }

  if (!is_connected_adj_mat(col_weight_matrix != 0)) {
    crv_error("Weights for columns do not imply a connected graph. Biclustering will not succeed.")
  }

  crv_message("Pre-computing row weights and edge sets")
  # Calculate row (observation)-clustering weights
  if (is.function(row_weights)) { # Usual case, `row_weights` is a function which calculates the weight matrix
    row_weight_result <- row_weights(X)

    if (is.matrix(row_weight_result)) {
      row_weight_matrix <- row_weight_result
      row_weight_type   <- UserFunction()
    } else {
      row_weight_matrix <- row_weight_result$weight_mat
      row_weight_type   <- row_weight_result$type
    }
  } else if (is.matrix(row_weights)) {

    if (!is_square(row_weights)) {
      crv_error(sQuote("row_weights"), " must be a square matrix.")
    }

    if (NROW(row_weights) != NROW(X)) {
      crv_error(sQuote("NROW(row_weights)"), " must be equal to ", sQuote("NROW(X)."))
    }

    row_weight_matrix <- row_weights
    row_weight_type   <- UserMatrix()
  } else {
    crv_error(sQuote("CBASS"), " does not know how to handle ", sQuote("row_weights"),
              " of class ", class(row_weights)[1], ".")
  }

  if (any(row_weight_matrix < 0) || anyNA(row_weight_matrix)) {
    crv_error("All row fusion weights must be positive or zero.")
  }

  if (!is_connected_adj_mat(row_weight_matrix != 0)) {
    crv_error("Weights for rows do not imply a connected graph. Biclustering will not succeed.")
  }

  row_weights <- weight_mat_to_vec(row_weight_matrix)
  col_weights <- weight_mat_to_vec(col_weight_matrix)

  ## Rescale to ensure coordinated fusions
  ##
  ## It is important that the observation (col) weights sum to 1/sqrt(p)
  ## and the feature/variable (row) weights sum to 1/sqrt(n).
  row_weights <- row_weights / (sum(row_weights) * sqrt(n))
  col_weights <- col_weights / (sum(col_weights) * sqrt(p))

  row_weight_matrix_ut <- row_weight_matrix * upper.tri(row_weight_matrix);

  row_edge_list <- which(row_weight_matrix_ut != 0, arr.ind = TRUE)
  row_edge_list <- row_edge_list[order(row_edge_list[, 1], row_edge_list[, 2]), ]
  num_edge_rows <- NROW(row_edge_list)
  D_row <- matrix(0, ncol = n, nrow = num_edge_rows)
  D_row[cbind(seq_len(num_edge_rows), row_edge_list[,1])] <-  1
  D_row[cbind(seq_len(num_edge_rows), row_edge_list[,2])] <- -1

  col_weight_matrix_ut <- col_weight_matrix * upper.tri(col_weight_matrix);

  col_edge_list <- which(col_weight_matrix_ut != 0, arr.ind = TRUE)
  col_edge_list <- col_edge_list[order(col_edge_list[, 1], col_edge_list[, 2]), ]
  num_edge_cols <- NROW(col_edge_list)
  D_col <- matrix(0, ncol = num_edge_cols, nrow = p)
  D_col[cbind(col_edge_list[,1], seq_len(num_edge_cols))] <-  1
  D_col[cbind(col_edge_list[,2], seq_len(num_edge_cols))] <- -1

  crv_message("Computing Convex Bi-Clustering [CBASS] Path")
  tic_inner <- Sys.time()

  cbass.sol.path <- CBASScpp(X,
                             D_row,
                             D_col,
                             t = t,
                             epsilon = .clustRvizOptionsEnv[["epsilon"]],
                             weights_row = row_weights[row_weights != 0],
                             weights_col = col_weights[col_weights != 0],
                             rho = .clustRvizOptionsEnv[["rho"]],
                             max_iter = .clustRvizOptionsEnv[["max_iter"]],
                             burn_in = .clustRvizOptionsEnv[["burn_in"]],
                             viz_max_inner_iter = .clustRvizOptionsEnv[["viz_max_inner_iter"]],
                             viz_initial_step = .clustRvizOptionsEnv[["viz_initial_step"]],
                             viz_small_step = .clustRvizOptionsEnv[["viz_small_step"]],
                             keep = .clustRvizOptionsEnv[["keep"]],
                             l1 = l1,
                             show_progress = status,
                             back_track = back_track,
                             exact = exact)

  toc_inner <- Sys.time()

  ## FIXME - Convert gamma.path to a single column matrix instead of a vector
  ##         RcppArmadillo returns a arma::vec as a n-by-1 matrix
  ##         RcppEigen returns an Eigen::VectorXd as a n-length vector
  ##         Something downstream cares about the difference, so just change
  ##         the type here for now
  cbass.sol.path$gamma_path <- matrix(cbass.sol.path$gamma_path, ncol=1)

  crv_message("Post-processing rows")

  post_processing_results_row <- ConvexClusteringPostProcess(X = X,
                                                             edge_matrix      = row_edge_list,
                                                             gamma_path       = cbass.sol.path$gamma_path,
                                                             u_path           = cbass.sol.path$u_path,
                                                             v_path           = cbass.sol.path$v_row_path,
                                                             v_zero_indices   = cbass.sol.path$v_row_zero_inds,
                                                             labels           = row_labels,
                                                             dendrogram_scale = dendrogram.scale,
                                                             npcs             = npcs)

  crv_message("Post-processing columns")

  post_processing_results_col <- ConvexClusteringPostProcess(X = t(X),
                                                             edge_matrix      = col_edge_list,
                                                             gamma_path       = cbass.sol.path$gamma_path,
                                                             u_path           = cbass.sol.path$u_path,
                                                             v_path           = cbass.sol.path$v_col_path,
                                                             v_zero_indices   = cbass.sol.path$v_col_zero_inds,
                                                             labels           = col_labels,
                                                             dendrogram_scale = dendrogram.scale,
                                                             npcs             = npcs,
                                                             internal_transpose = TRUE)

  cbass.fit <- list(
    X = X.orig,
    n = n,
    p = p,
    row_fusions = list(
      labels = row_labels,
      weight_type = row_weight_type,
      weights = row_weight_matrix,
      U = post_processing_results_row$U,
      D = D_row,
      dendrogram = post_processing_results_row$dendrogram,
      rotation_matrix = post_processing_results_row$rotation_matrix,
      cluster_membership = post_processing_results_row$membership_info
    ),
    col_fusions = list(
      labels = col_labels,
      weight_type = col_weight_type,
      weights = col_weight_matrix,
      U = post_processing_results_col$U,
      D = D_col,
      dendrogram = post_processing_results_col$dendrogram,
      rotation_matrix = post_processing_results_col$rotation_matrix,
      cluster_membership = post_processing_results_col$membership_info
    ),
    # General flags
    back_track = back_track,
    exact = exact,
    norm = norm,
    t = t,
    X.center.global = X.center.global,
    mean_adjust = mean_adjust,
    time = Sys.time() - tic,
    fit_time = toc_inner - tic_inner
  )

  if (.clustRvizOptionsEnv[["keep_debug_info"]]) {
    cbass.fit[["debug"]] <- list(path = cbass.sol.path,
                                 col  = post_processing_results_col[["debug"]],
                                 row  = post_processing_results_row[["debug"]])
  }

  class(cbass.fit) <- "CBASS"

  return(cbass.fit)
}

#' Print \code{CBASS} Results
#'
#' Prints a brief descripton of a fitted \code{CBASS} object.
#'
#' @details The \code{as.dendrogram} and \code{as.hclust} methods convert the
#'          \code{CBASS} output to an object of class \code{dendrogram} or \code{hclust}
#'          respectively.
#'
#' @param x an object of class \code{CARP} as returned by \code{\link{CARP}}
#' @param object an object of class \code{CARP} as returned by \code{\link{CARP}}
#' @param ... Additional unused arguments
#' @export
#' @rdname print_cbass
#' @examples
#' cbass_fit <- CBASS(X=presidential_speech)
#' print(cbass_fit)
print.CBASS <- function(x, ...) {
  if(x$exact){
    if(x$back_track){
      alg_string = "DLPA+ADMM-VIZ [Exact Solver + Back-Tracking Fusion Search]"
    } else {
      alg_string = paste0("DLPA+ADMM (t = ", round(x$t, 3), ") [Exact Solver]")
    }
  } else {
    if(x$back_track){
      alg_string = "CBASS-VIZ [Back-Tracking Fusion Search]"
    } else {
      alg_string = paste0("CBASS (t = ", round(x$t, 3), ")")
    }
  }

  if(x$norm == 1){
    alg_string <- paste(alg_string, "[L1]")
  }

  if(x$norm == 1){
    alg_string <- paste(alg_string, "[L1]")
  }

  cat("CBASS Fit Summary\n")
  cat("====================\n\n")
  cat("Algorithm:", alg_string, "\n")
  cat("Fit Time:", sprintf("%2.3f %s", x$fit_time, attr(x$fit_time, "units")), "\n")
  cat("Total Time:", sprintf("%2.3f %s", x$time, attr(x$time, "units")), "\n\n")

  cat("Number of Rows:", x$n, "\n")
  cat("Number of Columns:", x$p, "\n\n")

  cat("Pre-processing options:\n")
  cat(" - Global centering:", x$X.center.global, "\n\n")

  cat("Row Weights:\n")
  print(x$row_fusions$weight_type)

  cat("Column Weights:\n")
  print(x$col_fusions$weight_type)

  invisible(x)
}


#' @export
#' @importFrom stats as.dendrogram
#' @param type Either \code{"row"} or \code{"col"}, indicating whether
#'             the row or column dendrogram should be returned.
#' @rdname print_cbass
as.dendrogram.CBASS <- function(object, ..., type = c("row", "col")){
  type <- match.arg(type)

  if(type == "row"){
    as.dendrogram(object$row_fusions$dendrogram)
  } else {
    as.dendrogram(object$col_fusions$dendrogram)
  }
}

#' @importFrom stats as.hclust
#' @rdname print_cbass
as.hclust.CBASS <- function(x, ..., type = c("row", "col")){
  type <- match.arg(type)

  if(type == "row"){
    x$row_fusions$dendrogram
  } else {
    x$col_fusions$dendrogram
  }
}
