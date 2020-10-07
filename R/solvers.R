#' Compute Convex Clustering Solution Path on a User-Specified Grid
#'
#' \code{convex_clustering} calculates the convex clustering solution path
#' at a user-specified grid of lambda values (or just a single value). It is,
#' in general, difficult to know a useful set of lambda values \emph{a priori},
#' so this function is more useful for timing comparisons and methodological
#' research than applied work.
#'
#' Compared to the \code{\link{CARP}} function, the returned object
#' is much more "bare-bones," containing only the estimated \eqn{U} matrices,
#' and no information used for dendrogram or path visualizations.
#'
#' @param X The data matrix (\eqn{X \in R^{n \times p}}{X}): rows correspond to
#'          the observations (to be clustered) and columns to the variables (which
#'          will not be clustered). If \code{X} has missing values - \code{NA} or
#'          \code{NaN} values - they will be automatically imputed.
#' @param lambda_grid A user-supplied set of \eqn{\lambda}{lambda} values at which
#'                    to solve the convex clustering problem. These must be strictly
#'                    positive values and will be automatically sorted internally.
#' @param X.center A logical: Should \code{X} be centered columnwise?
#' @param X.scale A logical: Should \code{X} be scaled columnwise?
#' @param norm Which norm to use in the fusion penalty? Currently only \code{1}
#'             and \code{2} (default) are supported.
#' @param ... Unused arguements. An error will be thrown if any unrecognized
#'            arguments as given. All arguments other than \code{X} must be given
#'            by name.
#' @param weights One of the following: \itemize{
#'                \item A function which, when called with argument \code{X},
#'                      returns an b-by-n matrix of fusion weights.
#'                \item A matrix of size n-by-n containing fusion weights
#'                }
#' @param impute_func A function used to impute missing data in \code{X}. By default,
#'                    the \code{\link[missForest]{missForest}} function from the
#'                    package of the same name is used. This provides a flexible
#'                    potentially non-linear imputation function. This function
#'                    has to return a data matrix with no \code{NA} values.
#'                    Note that, consistent with base \code{R}, both \code{NaN}
#'                    and \code{NA} are treaded as "missing values" for imputation.
#' @param status Should a status message be printed to the console?
#' @return An object of class \code{convex_clustering} containing the following elements (among others):
#'         \itemize{
#'         \item \code{X}: the original data matrix
#'         \item \code{n}: the number of observations (rows of \code{X})
#'         \item \code{p}: the number of variables (columns of \code{X})
#'         \item \code{X.center}: a logical indicating whether \code{X} was centered
#'                                column-wise before clustering
#'         \item \code{X.scale}: a logical indicating whether \code{X} was scaled
#'                               column-wise before centering
#'         \item \code{weight_type}: a record of the scheme used to create
#'                                   fusion weights
#'         \item \code{U}: a tensor (3-array) of clustering solutions
#'         }
#' @importFrom utils data
#' @importFrom dplyr %>% mutate group_by ungroup as_tibble n_distinct
#' @importFrom rlang %||%
#' @importFrom stats var
#' @importFrom missForest missForest
#' @export
#' @examples
#' clustering_fit <- convex_clustering(presidential_speech[1:10,1:4], lambda_grid = 1:100)
#' print(clustering_fit)
convex_clustering <- function(X,
                              ...,
                              lambda_grid,
                              weights = sparse_rbf_kernel_weights(k = "auto",
                                                                  phi = "auto",
                                                                  dist.method = "euclidean",
                                                                  p = 2),
                              X.center = TRUE,
                              X.scale = FALSE,
                              norm = 2,
                              impute_func = function(X) {if(anyNA(X)) missForest(X)$ximp else X},
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
      crv_error("Unknown argument ", sQuote(names(dots)[1L]), " passed to ", sQuote("convex_clustering."))
    } else {
      crv_error("Unknown ", sQuote("..."), " arguments passed to ", sQuote("convex_clustering."))
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

  # Missing data mask: M_{ij} = 1 means we see X_{ij};
  M <- 1 - is.na(X)

  # Impute missing values in X
  # By default, we use the "Missing Forest" function from the missForest package
  # though other imputers can be supplied by the user.
  X.orig <- X

  if(anyNA(X)) {
    X <- impute_func(X)
  }

  ## Check that imputation was successful.
  if (anyNA(X)) {
    crv_error("Imputation failed. Missing values found in ", sQuote("X"), " even after imputation.")
  }

  if (!all(is.finite(X))) {
    crv_error("All elements of ", sQuote("X"), " must be finite.")
  }

  if (!is_logical_scalar(X.center)) {
    crv_error(sQuote("X.center"), "must be either ", sQuote("TRUE"), " or ", sQuote("FALSE."))
  }

  if (!is_logical_scalar(X.scale)) {
    crv_error(sQuote("X.scale"), "must be either ", sQuote("TRUE"), " or ", sQuote("FALSE."))
  }

  if (norm %not.in% c(1, 2)){
    crv_error(sQuote("norm"), " must be either 1 or 2.")
  }

  if(missing(lambda_grid)){
    crv_error(sQuote("lambda_grid"), " must be supplied.")
  }

  if(!all(is.finite(lambda_grid))){
    crv_error(sQuote("lambda_grid"), " containts non-finite values.")
  }

  if(any(lambda_grid < 0)){
    crv_error(sQuote("lambda_grid"), " must be strictly positive.")
  }

  if(any(lambda_grid == 0)){
    crv_error(sQuote("lambda_grid"), " must be strictly positive - 0 will be automatically added.")
  }

  if(is.unsorted(lambda_grid)){
    crv_warning(sQuote("lambda_grid"), " is unsorted: sorting for maximum performance.")
    lambda_grid <- sort(lambda_grid)
  }

  l1 <- (norm == 1)

  n <- NROW(X)
  p <- NCOL(X)

  # Center and scale X
  if (X.center | X.scale) {
    X <- scale(X, center = X.center, scale = X.scale)
  }

  scale_vector  <- attr(X, "scaled:scale", exact=TRUE)  %||% rep(1, p)
  center_vector <- attr(X, "scaled:center", exact=TRUE) %||% rep(0, p)

  crv_message("Pre-computing weights and edge sets")

  # Calculate clustering weights
  if (is.function(weights)) { # Usual case, `weights` is a function which calculates the weight matrix
    weight_result <- weights(X)

    if (is.matrix(weight_result)) {
      weight_matrix <- weight_result
      weight_type   <- UserFunction()
    } else {
      weight_matrix <- weight_result$weight_mat
      weight_type   <- weight_result$type
    }
  } else if (is.matrix(weights)) {

    if (!is_square(weights)) {
      crv_error(sQuote("weights"), " must be a square matrix.")
    }

    if (NROW(weights) != NROW(X)) {
      crv_error(sQuote("NROW(weights)"), " must be equal to ", sQuote("NROW(X)."))
    }

    weight_matrix <- weights
    weight_type   <- UserMatrix()
  } else {
    crv_error(sQuote("convex_clustering"), " does not know how to handle ", sQuote("weights"),
              " of class ", class(weights)[1], ".")
  }

  if (any(weight_matrix < 0) || anyNA(weight_matrix)) {
    crv_error("All fusion weights must be positive or zero.")
  }

  weight_matrix_ut <- weight_matrix * upper.tri(weight_matrix);

  edge_list <- which(weight_matrix_ut != 0, arr.ind = TRUE)
  edge_list <- edge_list[order(edge_list[, 1], edge_list[, 2]), ]
  cardE <- NROW(edge_list)
  D <- matrix(0, ncol = n, nrow = cardE)
  D[cbind(seq_len(cardE), edge_list[,1])] <-  1
  D[cbind(seq_len(cardE), edge_list[,2])] <- -1

  weight_vec <- weight_mat_to_vec(weight_matrix)

  crv_message("Computing Convex Clustering Solutions")
  tic_inner <- Sys.time()

  clustering_sol <- ConvexClusteringCPP(X = X,
                                        M = M,
                                        D = D,
                                        lambda_grid = lambda_grid,
                                        weights = weight_vec[weight_vec != 0],
                                        rho = .clustRvizOptionsEnv[["rho"]],
                                        thresh = .clustRvizOptionsEnv[["stopping_threshold"]],
                                        max_iter = .clustRvizOptionsEnv[["max_iter"]],
                                        max_inner_iter = .clustRvizOptionsEnv[["max_inner_iter"]],
                                        l1 = l1,
                                        show_progress = status)

  toc_inner <- Sys.time()

  crv_message("Post-processing")

  lambda_grid <- clustering_sol$gamma_path
  U_raw <- array(clustering_sol$u_path,
                 dim = c(n, p, length(lambda_grid)),
                 dimnames = list(rownames(X.orig),
                                 colnames(X.orig),
                                 paste0("Lambda_", seq_along(lambda_grid) - 1)))

  center_tensor <- aperm(array(center_vector, dim = c(p, n, length(lambda_grid))), c(2, 1, 3))
  scale_tensor  <- aperm(array(scale_vector, dim = c(p, n, length(lambda_grid))), c(2, 1, 3))

  U <- U_raw * scale_tensor + center_tensor

  convex_clustering_fit <- list(
    X = X.orig,
    M = M,
    D = D,
    U = U,
    n = n,
    p = p,
    norm = norm,
    lambda_grid = lambda_grid,
    weights = weight_matrix,
    weight_type = weight_type,
    X.center = X.center,
    center_vector = center_vector,
    X.scale = X.scale,
    scale_vector = scale_vector,
    time = Sys.time() - tic,
    fit_time = toc_inner - tic_inner
  )

  class(convex_clustering_fit) <- "ConvexClustering"

  return(convex_clustering_fit)
}

#' Print \code{\link{convex_clustering}} Results
#'
#' Prints a brief descripton of a fitted \code{convex_clustering} object.
#'
#' Reports number of observations and variables of dataset, any preprocessing
#' done by the \code{\link{convex_clustering}} function, and regularization
#' details.
#'
#' @param x an object of class \code{convex_clustering} as returned by
#'          \code{\link{convex_clustering}}
#' @param ... Additional unused arguments
#' @export
#' @examples
#' cluster_fit <- convex_clustering(presidential_speech, lambda_grid = 1:5)
#' print(cluster_fit)
print.ConvexClustering <- function(x, ...) {
  alg_string <- paste("ADMM", if(x$norm == 1) "[L1]" else "[L2]")

  cat("Convex Clustering Fit Summary\n")
  cat("=============================\n\n")
  cat("Algorithm:", alg_string, "\n")
  cat("Grid:", length(x$lambda_grid), "values of lambda. \n")
  cat("Fit Time:", sprintf("%2.3f %s", x$fit_time, attr(x$fit_time, "units")), "\n")
  cat("Total Time:", sprintf("%2.3f %s", x$time, attr(x$time, "units")), "\n\n")
  cat("Number of Observations:", x$n, "\n")
  cat("Number of Variables:   ", x$p, "\n\n")

  cat("Pre-processing options:\n")
  cat(" - Columnwise centering:", x$X.center, "\n")
  cat(" - Columnwise scaling:  ", x$X.scale, "\n\n")

  cat("Weights:\n")
  print(x$weight_type)

  invisible(x)
}

#' Compute Convex BiClustering Solution Path on a User-Specified Grid
#'
#' \code{convex_biclustering} calculates the convex biclustering solution path
#' at a user-specified grid of lambda values (or just a single value). It is,
#' in general, difficult to know a useful set of lambda values \emph{a priori},
#' so this function is more useful for timing comparisons and methodological
#' research than applied work.
#'
#' Compared to the \code{\link{CBASS}} function, the returned object
#' is much more "bare-bones," containing only the estimated \eqn{U} matrices,
#' and no information used for dendrogram or path visualizations.
#'
#' @param X The data matrix (\eqn{X \in R^{n \times p}}{X}).
#'          If \code{X} has missing values - \code{NA} or
#'          \code{NaN} values - they will be automatically imputed.
#' @param lambda_grid A user-supplied set of \eqn{\lambda}{lambda} values at which
#'                    to solve the convex biclustering problem. These must be strictly
#'                    positive values and will be automatically sorted internally.
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
#' @param X.center.global A logical: Should \code{X} be centered globally?
#'                        \emph{I.e.}, should the global mean of \code{X} be subtracted?
#' @param norm Which norm to use in the fusion penalty? Currently only \code{1}
#'             and \code{2} (default) are supported.
#' @param ... Unused arguements. An error will be thrown if any unrecognized
#'            arguments as given.
#' @param status Should a status message be printed to the console?
#' @return An object of class \code{convex_biclustering} containing the
#'         following elements (among others):
#'         \itemize{
#'         \item \code{X}: the original data matrix
#'         \item \code{n}: the number of observations (rows of \code{X})
#'         \item \code{p}: the number of variables (columns of \code{X})
#'         \item \code{U}: a tensor (3-array) of clustering solutions
#'         }
#' @export
#' @examples
#' \dontrun{
#' biclustering_fit <- convex_biclustering(presidential_speech, lambda_grid = 1:100)
#' print(biclustering_fit)
#' }
convex_biclustering <- function(X,
                                ...,
                                lambda_grid,
                                row_weights = sparse_rbf_kernel_weights(k = "auto",
                                                                        phi = "auto",
                                                                        dist.method = "euclidean",
                                                                        p = 2),
                                col_weights = sparse_rbf_kernel_weights(k = "auto",
                                                                        phi = "auto",
                                                                        dist.method = "euclidean",
                                                                        p = 2),
                                X.center.global = TRUE,
                                norm = 2,
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
      crv_error("Unknown argument ", sQuote(names(dots)[1L]), " passed to ", sQuote("convex_biclustering."))
    } else {
      crv_error("Unknown ", sQuote("..."), " arguments passed to ", sQuote("convex_biclustering."))
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

  # Missing data mask: M_{ij} = 1 means we see X_{ij};
  M <- 1 - is.na(X)

  # Impute missing values in X via the global mean
  X.orig <- X

  if(anyNA(X)) {
    X[is.na(X)] <- mean(X, na.rm = TRUE)
  }

  ## Check that imputation was successful.
  if (anyNA(X)) {
    crv_error("Imputation failed. Missing values found in ", sQuote("X"), " even after imputation.")
  }

  if (!all(is.finite(X))) {
    crv_error("All elements of ", sQuote("X"), " must be finite.")
  }

  if (!is.logical(X.center.global) || anyNA(X.center.global) || (length(X.center.global) != 1L)) {
    crv_error(sQuote("X.center.global"), "must be either ", sQuote("TRUE"), " or ", sQuote("FALSE."))
  }

  if (norm %not.in% c(1, 2)){
    crv_error(sQuote("norm"), " must be either 1 or 2.")
  }

  l1 <- (norm == 1)

  if(missing(lambda_grid)){
    crv_error(sQuote("lambda_grid"), " must be supplied.")
  }

  if(!all(is.finite(lambda_grid))){
    crv_error(sQuote("lambda_grid"), " containts non-finite values.")
  }

  if(any(lambda_grid < 0)){
    crv_error(sQuote("lambda_grid"), " must be strictly positive.")
  }

  if(any(lambda_grid == 0)){
    crv_error(sQuote("lambda_grid"), " must be strictly positive - 0 will be automatically added.")
  }

  if(is.unsorted(lambda_grid)){
    crv_warning(sQuote("lambda_grid"), " is unsorted: sorting for maximum performance.")
    lambda_grid <- sort(lambda_grid)
  }

  n <- NROW(X)
  p <- NCOL(X)

  # Preprocess X
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
    crv_error(sQuote("convex_biclustering"), " does not know how to handle ", sQuote("col_weights"),
              " of class ", class(col_weights)[1], ".")
  }

  if (any(col_weight_matrix < 0) || anyNA(col_weight_matrix)) {
    crv_error("All column fusion weights must be positive or zero.")
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
    crv_error(sQuote("convex_biclustering"), " does not know how to handle ", sQuote("row_weights"),
              " of class ", class(row_weights)[1], ".")
  }

  if (any(row_weight_matrix < 0) || anyNA(row_weight_matrix)) {
    crv_error("All row fusion weights must be positive or zero.")
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

  crv_message("Computing Convex Bi-Clustering Path")
  tic_inner <- Sys.time()

  biclustering_sol <- ConvexBiClusteringCPP(X = X,
                                            M = M,
                                            D_row = D_row,
                                            D_col = D_col,
                                            lambda_grid = lambda_grid,
                                            weights_row = row_weights[row_weights != 0],
                                            weights_col = col_weights[col_weights != 0],
                                            rho = .clustRvizOptionsEnv[["rho"]],
                                            thresh = .clustRvizOptionsEnv[["stopping_threshold"]],
                                            max_iter = .clustRvizOptionsEnv[["max_iter"]],
                                            max_inner_iter = .clustRvizOptionsEnv[["max_inner_iter"]],
                                            l1 = l1,
                                            show_progress = status)

  toc_inner <- Sys.time()

  crv_message("Post-processing")

  lambda_grid <- biclustering_sol$gamma_path
  U <- array(biclustering_sol$u_path,
             dim = c(n, p, length(lambda_grid)),
             dimnames = list(rownames(X.orig),
                             colnames(X.orig),
                             paste0("Lambda_", seq_along(lambda_grid) - 1))) +
           mean_adjust

  convex_biclustering_fit <- list(
    X = X.orig,
    M = M,
    D_row = D_row,
    D_col = D_col,
    U = U,
    n = n,
    p = p,
    norm = norm,
    lambda_grid = lambda_grid,
    row_weights = row_weight_matrix,
    row_weight_type = row_weight_type,
    col_weights = col_weight_matrix,
    col_weight_type = col_weight_type,
    X.center.global = X.center.global,
    mean_adjust = mean_adjust,
    time = Sys.time() - tic,
    fit_time = toc_inner - tic_inner
  )

  class(convex_biclustering_fit) <- "ConvexBiClustering"

  return(convex_biclustering_fit)
}

#' Print \code{\link{convex_biclustering}} Results
#'
#' Prints a brief descripton of a fitted \code{convex_biclustering} object.
#'
#' Reports number of observations and variables of dataset, any preprocessing
#' done by the \code{\link{convex_biclustering}} function, and regularization
#' details.
#'
#' @param x an object of class \code{convex_biclustering} as returned by
#'          \code{\link{convex_biclustering}}
#' @param ... Additional unused arguments
#' @export
#' @examples
#' bicluster_fit <- convex_biclustering(presidential_speech, lambda_grid = 1:5)
#' print(bicluster_fit)
print.ConvexBiClustering <- function(x, ...) {
  alg_string <- paste("ADMM", if(x$norm == 1) "[L1]" else "[L2]")

  cat("Convex Clustering Fit Summary\n")
  cat("=============================\n\n")
  cat("Algorithm:", alg_string, "\n")
  cat("Grid:", length(x$lambda_grid), "values of lambda. \n")
  cat("Fit Time:", sprintf("%2.3f %s", x$fit_time, attr(x$fit_time, "units")), "\n")
  cat("Total Time:", sprintf("%2.3f %s", x$time, attr(x$time, "units")), "\n\n")

  cat("Number of Rows:", x$n, "\n")
  cat("Number of Columns:", x$p, "\n\n")

  cat("Pre-processing options:\n")
  cat(" - Global centering:", x$X.center.global, "\n\n")

  cat("Row Weights:\n")
  print(x$row_weight_type)

  cat("Column Weights:\n")
  print(x$col_weight_type)

  invisible(x)
}

#' Compute Time Series (Phase Aligned Complex) Convex Clustering Solution Path on a User-Specified Grid
#'
#' \code{ts_convex_clustering} calculates the phase-aligned complex convex clustering solution path
#' at a user-specified grid of lambda values (or just a single value). It is,
#' in general, difficult to know a useful set of lambda values \emph{a priori},
#' so this function is more useful for timing comparisons and methodological
#' research than applied work.
#'
#' Compared to the \code{\link{TROUT}} function, the returned object
#' is much more "bare-bones," containing only the estimated \eqn{U} matrices,
#' and no information used for dendrogram or path visualizations.
#'
#' @param X The data matrix (\eqn{X \in C^{n \times p}}{X}): rows correspond to
#'          the observations (to be clustered) and columns to the variables (which
#'          will not be clustered).
#' @param lambda_grid A user-supplied set of \eqn{\lambda}{lambda} values at which
#'                    to solve the convex clustering problem. These must be strictly
#'                    positive values and will be automatically sorted internally.
#' @param X.center A logical: Should \code{X} be centered columnwise?
#' @param X.scale A logical: Should \code{X} be scaled columnwise?
#' @param norm Which norm to use in the fusion penalty? Currently only \code{1}
#'             and \code{2} (default) are supported.
#' @param ... Unused arguements. An error will be thrown if any unrecognized
#'            arguments as given. All arguments other than \code{X} must be given
#'            by name.
#' @param weights One of the following: \itemize{
#'                \item A function which, when called with argument \code{X},
#'                      returns an b-by-n matrix of fusion weights.
#'                \item A matrix of size n-by-n containing fusion weights
#'                }
#' @param status Should a status message be printed to the console?
#' @return An object of class \code{convex_clustering} containing the following elements (among others):
#'         \itemize{
#'         \item \code{X}: the original data matrix
#'         \item \code{n}: the number of observations (rows of \code{X})
#'         \item \code{p}: the number of variables (columns of \code{X})
#'         \item \code{X.center}: a logical indicating whether \code{X} was centered
#'                                column-wise before clustering
#'         \item \code{X.scale}: a logical indicating whether \code{X} was scaled
#'                               column-wise before centering
#'         \item \code{weight_type}: a record of the scheme used to create
#'                                   fusion weights
#'         \item \code{U}: a tensor (3-array) of clustering solutions
#'         }
#' @importFrom utils data
#' @importFrom dplyr %>% mutate group_by ungroup as_tibble n_distinct
#' @importFrom rlang %||%
#' @importFrom stats var
#' @export
ts_convex_clustering <- function(X,
                                 ...,
                                 lambda_grid,
                                 weights = sparse_rbf_kernel_weights(k = "auto",
                                                                     phi = "auto",
                                                                     dist.method = "trout",
                                                                     p = 2),
                                 X.center = TRUE,
                                 X.scale = FALSE,
                                 norm = 2,
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
      crv_error("Unknown argument ", sQuote(names(dots)[1L]), " passed to ", sQuote("ts_convex_clustering."))
    } else {
      crv_error("Unknown ", sQuote("..."), " arguments passed to ", sQuote("ts_convex_clustering."))
    }
  }

  if (!is.matrix(X)) {
    crv_warning(sQuote("X"), " should be a matrix, not a " , class(X)[1],
                ". Converting with as.matrix().")
    X <- as.matrix(X)
  }

  if (!is.complex(X)) {
    crv_error(sQuote("X"), " must be complex.")
  }

  # Missing data mask: M_{ij} = 1 means we see X_{ij};
  # Currently we don't support missing values in X so this better be all ones
  M <- 1 - is.na(X)

  ## Check that imputation was successful.
  if (anyNA(X)) {
    crv_error("ts_convex_clustering() does not support missing values in ", sQuote("X."))
  }

  if (!all(is.finite(X))) {
    crv_error("All elements of ", sQuote("X"), " must be finite.")
  }

  if (!is_logical_scalar(X.center)) {
    crv_error(sQuote("X.center"), "must be either ", sQuote("TRUE"), " or ", sQuote("FALSE."))
  }

  if (!is_logical_scalar(X.scale)) {
    crv_error(sQuote("X.scale"), "must be either ", sQuote("TRUE"), " or ", sQuote("FALSE."))
  }

  if (norm %not.in% c(1, 2)){
    crv_error(sQuote("norm"), " must be either 1 or 2.")
  }

  if(missing(lambda_grid)){
    crv_error(sQuote("lambda_grid"), " must be supplied.")
  }

  if(!all(is.finite(lambda_grid))){
    crv_error(sQuote("lambda_grid"), " containts non-finite values.")
  }

  if(any(lambda_grid < 0)){
    crv_error(sQuote("lambda_grid"), " must be strictly positive.")
  }

  if(any(lambda_grid == 0)){
    crv_error(sQuote("lambda_grid"), " must be strictly positive - 0 will be automatically added.")
  }

  if(is.unsorted(lambda_grid)){
    crv_warning(sQuote("lambda_grid"), " is unsorted: sorting for maximum performance.")
    lambda_grid <- sort(lambda_grid)
  }

  l1 <- (norm == 1)

  n <- NROW(X)
  p <- NCOL(X)

  X.orig <- X
  # Center and scale X
  if (X.center | X.scale) {
    X <- scale(X, center = X.center, scale = X.scale)
  }

  scale_vector  <- attr(X, "scaled:scale", exact=TRUE)  %||% rep(1, p)
  center_vector <- attr(X, "scaled:center", exact=TRUE) %||% rep(0, p)

  crv_message("Pre-computing weights and edge sets")

  # Calculate clustering weights
  if (is.function(weights)) { # Usual case, `weights` is a function which calculates the weight matrix
    weight_result <- weights(X)

    if (is.matrix(weight_result)) {
      weight_matrix <- weight_result
      weight_type   <- UserFunction()
    } else {
      weight_matrix <- weight_result$weight_mat
      weight_type   <- weight_result$type
    }
  } else if (is.matrix(weights)) {

    if (!is_square(weights)) {
      crv_error(sQuote("weights"), " must be a square matrix.")
    }

    if (NROW(weights) != NROW(X)) {
      crv_error(sQuote("NROW(weights)"), " must be equal to ", sQuote("NROW(X)."))
    }

    weight_matrix <- weights
    weight_type   <- UserMatrix()
  } else {
    crv_error(sQuote("ts_convex_clustering"), " does not know how to handle ", sQuote("weights"),
              " of class ", class(weights)[1], ".")
  }

  if (any(weight_matrix < 0) || anyNA(weight_matrix)) {
    crv_error("All fusion weights must be positive or zero.")
  }

  weight_matrix_ut <- weight_matrix * upper.tri(weight_matrix);

  edge_list <- which(weight_matrix_ut != 0, arr.ind = TRUE)
  edge_list <- edge_list[order(edge_list[, 1], edge_list[, 2]), ]
  cardE <- NROW(edge_list)
  D <- matrix(0, ncol = n, nrow = cardE)
  D[cbind(seq_len(cardE), edge_list[,1])] <-  1
  D[cbind(seq_len(cardE), edge_list[,2])] <- -1

  weight_vec <- weight_mat_to_vec(weight_matrix)

  crv_message("Computing Time-Series Convex Clustering Solutions")
  tic_inner <- Sys.time()

  clustering_sol <- TroutClusteringCPP(X = X,
                                       M = M,
                                       D = D,
                                       lambda_grid = lambda_grid,
                                       weights = weight_vec[weight_vec != 0],
                                       rho = .clustRvizOptionsEnv[["rho"]],
                                       thresh = .clustRvizOptionsEnv[["stopping_threshold"]],
                                       max_iter = .clustRvizOptionsEnv[["max_iter"]],
                                       max_inner_iter = .clustRvizOptionsEnv[["max_inner_iter"]],
                                       l1 = l1,
                                       show_progress = status)

  toc_inner <- Sys.time()

  crv_message("Post-processing")

  lambda_grid <- clustering_sol$gamma_path
  U_raw <- array(clustering_sol$u_path,
                 dim = c(n, p, length(lambda_grid)),
                 dimnames = list(rownames(X.orig),
                                 colnames(X.orig),
                                 paste0("Lambda_", seq_along(lambda_grid) - 1)))

  center_tensor <- aperm(array(center_vector, dim = c(p, n, length(lambda_grid))), c(2, 1, 3))
  scale_tensor  <- aperm(array(scale_vector, dim = c(p, n, length(lambda_grid))), c(2, 1, 3))

  U <- U_raw * scale_tensor + center_tensor

  trout_clustering_fit <- list(
    X = X.orig,
    M = M,
    D = D,
    U = U,
    n = n,
    p = p,
    norm = norm,
    lambda_grid = lambda_grid,
    weights = weight_matrix,
    weight_type = weight_type,
    X.center = X.center,
    center_vector = center_vector,
    X.scale = X.scale,
    scale_vector = scale_vector,
    time = Sys.time() - tic,
    fit_time = toc_inner - tic_inner
  )

  class(trout_clustering_fit) <- "TimeSeriesConvexClustering"

  return(trout_clustering_fit)
}

#' Print \code{\link{ts_convex_clustering}} Results
#'
#' Prints a brief descripton of a fitted \code{ts_convex_clustering} object.
#'
#' Reports number of observations and variables of dataset, any preprocessing
#' done by the \code{\link{ts_convex_clustering}} function, and regularization
#' details.
#'
#' @param x an object of class \code{ts_convex_clustering} as returned by
#'          \code{\link{ts_convex_clustering}}
#' @param ... Additional unused arguments
#' @export
print.TimeSeriesConvexClustering <- function(x, ...) {
  alg_string <- paste("ADMM", if(x$norm == 1) "[L1]" else "[L2]")

  cat("Time Series [Phase-Aligned Convex] Convex Clustering Fit Summary\n")
  cat("================================================================\n\n")
  cat("Algorithm:", alg_string, "\n")
  cat("Grid:", length(x$lambda_grid), "values of lambda. \n")
  cat("Fit Time:", sprintf("%2.3f %s", x$fit_time, attr(x$fit_time, "units")), "\n")
  cat("Total Time:", sprintf("%2.3f %s", x$time, attr(x$time, "units")), "\n\n")
  cat("Number of Observations:", x$n, "\n")
  cat("Number of Variables:   ", x$p, "\n\n")

  cat("Pre-processing options:\n")
  cat(" - Columnwise centering:", x$X.center, "\n")
  cat(" - Columnwise scaling:  ", x$X.scale, "\n\n")

  cat("Weights:\n")
  print(x$weight_type)

  invisible(x)
}
