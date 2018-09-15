#' Compute \code{CBASS} (Convex BiClustering) Solution Path
#'
#' \code{CBASS} returns a fast approximation to the Convex BiClustering
#' solution path along with visualizations such as dendrograms and
#' heatmaps. \code{CBASS} solves the Convex Biclustering problem via
#' Algorithmic Regularization Paths. A seqeunce of biclustering
#' solutions is returned along with several visualizations.
#'
#' @param X The data matrix (\eqn{X \in R^{n \times p}}{X}): rows correspond to
#'          the observations (to be clustered) and columns to the variables (which
#'          will not be clustered).
#' @param verbose Any of the values \code{0}, \code{1}, or \code{2}. Higher values
#'                correspond to more verbose output while running.
#' @param obs_weights One of the following: \itemize{
#'                    \item A function which, when called with argument \code{X},
#'                          returns a n-by-n matrix of fusion weights.
#'                    \item A matrix of size n-by-ncontaining fusion weights
#'                    }
#'                    Note that the weights will be renormalized to sum to
#'                    \eqn{1/\sqrt{n}} internally.
#' @param var_weights One of the following: \itemize{
#'                    \item A function which, when called with argument \code{t(X)},
#'                          returns a p-by-p matrix of fusion weights. (Note the
#'                          transpose.)
#'                    \item A matrix of size p-by-p containing fusion weights
#'                    }
#'                    Note that the weights will be renormalized to sum to
#'                    \eqn{1/\sqrt{p}} internally.
#' @param obs_labels A character vector of length \eqn{n}: observations (row) labels
#' @param var_labels A character vector of length \eqn{p}: variable (column) labels
#' @param X.center.global A logical: Should \code{X} be centered globally?
#'                        \emph{I.e.}, should the global mean of \code{X} be subtracted?
#' @param rho For advanced users only (not advisable to change): the penalty
#'            parameter used for the augmented Lagrangian.
#' @param max.iter An integer: the maximum number of \code{CBASS} iterations to perform.
#' @param burn.in An integer: the number of initial iterations at a fixed
#'                (small) value of \eqn{\lambda}
#' @param alg.type Which \code{CBASS} variant to use. Allowed values are \itemize{
#'        \item \code{"cbass"} - The standard \code{CBASS} algorithm with \eqn{L2} penalty;
#'        \item \code{"cbassviz"} - The back-tracking \code{CBASS} algorithm with \eqn{L2} penalty;
#'        \item \code{"cbassl1"} - The standard \code{CBASS} algorithm with \eqn{L1} penalty; and
#'        \item \code{"cbassvizl1"} - The back-tracking \code{CBASS} algorithm with \eqn{L1} penalty.}
#' @param t A number greater than 1: the size of the multiplicative update to
#'          the cluster fusion regularization parameter (not used by
#'          back-tracking variants). Typically on the scale of \code{1.005} to \code{1.1}.
#' @param npcs An integer >= 2. The number of principal components to compute
#'             for path visualization.
#' @param ... Unused arguements. An error will be thrown if any unrecognized
#'            arguments as given.
#' @return An object of class \code{CBASS} containing the following elements (among others):
#'         \itemize{
#'         \item \code{X}: the original data matrix
#'         \item \code{n.obs}: the number of observations (rows of \code{X})
#'         \item \code{p.var}: the number of variables (columns of \code{X})
#'         \item \code{alg.type}: the \code{CBASS} variant used
#'         \item \code{X.center.global}: a logical indicating whether \code{X}
#'                                       was globally centered before clustering
#'         \item \code{var_weight_type}: a record of the scheme used to create
#'                                       fusion weights for the variables
#'         \item \code{obs_weight_type}: a record of the scheme used to create
#'                                       fusion weights for the observations
#'         \item \code{burn.in}: an integer indicating the number of "burn-in"
#'                               iterations performed
#'         \item \code{carp.dend.obs}: a dendrogram (object of class
#'                                     \code{\link[stats]{hclust}}) containing
#'                                     the clustering solution path for the observations
#'         \item \code{carp.dend.var}: a dendrogram (object of class
#'                                     \code{\link[stats]{hclust}}) containing
#'                                     the clustering solution path for the variables
#'         \item \code{cbass.cluster.path.obs}: The \code{CBASS} solution path for the observations
#'         \item \code{cbass.cluster.path.var}: The \code{CBASS} solution path for the variables
#'         \item \code{obs.labels}: a character vector of length \code{n.obs} containing
#'                                  observation (row) labels
#'         \item \code{var.labels}: a character vector of length \code{p.var} containing
#'                                  variable (column) labels
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
                  verbose = 1L,
                  var_weights = sparse_rbf_kernel_weights(k = "auto",
                                                          phi = "auto",
                                                          dist.method = "euclidean",
                                                          p = 2),
                  obs_weights = sparse_rbf_kernel_weights(k = "auto",
                                                          phi = "auto",
                                                          dist.method = "euclidean",
                                                          p = 2),
                  obs_labels = rownames(X),
                  var_labels = colnames(X),
                  X.center.global = TRUE,
                  rho = 1.0,
                  t = 1.01,
                  max.iter = 1000000L,
                  burn.in = 50L,
                  alg.type = c("cbassviz", "cbassvizl1", "cbass", "cbassl1"),
                  npcs = min(4L, NCOL(X))) {

  ####################
  ##
  ## Input validation
  ##
  ####################

  dots <- list(...)

  if (length(dots) != 0L) {
    if (!is.null(names(dots))) {
      stop("Unknown argument ", sQuote(names(dots)[1L]), " passed to ", sQuote("CBASS."))
    } else {
      stop("Unknown ", sQuote("..."), " arguments passed to ", sQuote("CBASS."))
    }
  }

  if (!is.matrix(X)) {
    warning(sQuote("X"), " should be a matrix, not a " , class(X)[1],
            ". Converting with as.matrix().")
    X <- as.matrix(X)
  }

  if (!is.numeric(X)) {
    stop(sQuote("X"), " must be numeric.")
  }

  if (anyNA(X)) {
    stop(sQuote("CBASS"), " cannot handle missing data.")
  }

  if (!all(is.finite(X))) {
    stop("All elements of ", sQuote("X"), " must be finite.")
  }

  if (!is.logical(X.center.global) || is.na(X.center.global) || (length(X.center.global) != 1L)) {
    stop(sQuote("X.center.global"), "must be either ", sQuote("TRUE"), " or ", sQuote("FALSE."))
  }

  if ( (!is_numeric_scalar(rho)) || (rho <= 0)) {
    stop(sQuote("rho"), "must be a positive scalar (vector of length 1).")
  }

  if ( (!is_integer_scalar(npcs)) || (npcs < 2) || (npcs > NCOL(X)) ){
    stop(sQuote("npcs"), " must be an integer scalar between 2 and ", sQuote("NCOL(X)."))
  }

  if ( (!is_integer_scalar(max.iter)) || (max.iter <= 1L) ) {
    stop(sQuote("max.iter"), " must be a positive integer scalar and at least 2.")
  }

  if ( (!is_integer_scalar(burn.in)) || (burn.in <= 0L) || (burn.in >= max.iter) ) {
    stop(sQuote("burn.in"), " must be a positive integer less than ", sQuote("max.iter."))
  }

  alg.type <- match.arg(alg.type)

  if ( (!is_numeric_scalar(t)) || (t <= 1) ) {
    stop(sQuote("t"), " must be a scalar greater than 1.")
  }

  ## Get row (observation) labels
  if (is.null(obs_labels)) {
    obs_labels <- paste0("Obs", seq_len(NROW(X)))
  }

  if ( length(obs_labels) != NROW(X) ){
    stop(sQuote("obs_labels"), " must be of length ", sQuote("NROW(X)."))
  }

  rownames(X) <- obs_labels <- make.unique(as.character(obs_labels), sep = "_")

  ## Get column (variable) labels
  if (is.null(var_labels)) {
    var_labels <- paste0("Var", seq_len(NCOL(X)))
  }

  if ( length(var_labels) != NCOL(X) ){
    stop(sQuote("var_labels"), " must be of length ", sQuote("NCOL(X)."))
  }

  colnames(X) <- var_labels <- make.unique(as.character(var_labels), sep = "_")

  if (is.logical(verbose)) {
    verbose.basic <- TRUE
    verbose.deep <- FALSE
  } else if (verbose == 1) {
    verbose.basic <- TRUE
    verbose.deep <- FALSE
  } else if (verbose == 2) {
    verbose.basic <- TRUE
    verbose.deep <- TRUE
  } else {
    verbose.basic <- FALSE
    verbose.deep <- FALSE
  }

  n.obs <- NROW(X)
  p.var <- NCOL(X)

  # Preprocess X
  X.orig <- X
  if (X.center.global) {
    X <- X - mean(X)
  }

  ## Transform to a form suitable for down-stream computation
  X <- t(X) ## TODO: Ask JN why we keep this. (It matches the convention in Chi, Allen, Baraniuk)

  # Calculate variable/feature (row)-clustering weights
  if (is.function(var_weights)) { # Usual case, `var_weights` is a function which calculates the weight matrix
    var_weight_result <- var_weights(X)

    if (is.matrix(var_weight_result)) {
      var_weight_matrix <- var_weight_result
      var_weight_type   <- UserFunction()
    } else {
      var_weight_matrix <- var_weight_result$weight_mat
      var_weight_type   <- var_weight_result$type
    }
  } else if (is.matrix(var_weights)) {

    if (!is_square(var_weights)) {
      stop(sQuote("var_weights"), " must be a square matrix.")
    }

    if (NROW(var_weights) != NROW(X)) {
      stop(sQuote("NROW(var_weights)"), " must be equal to ", sQuote("NROW(X)."))
    }

    var_weight_matrix <- var_weights
    var_weight_type   <- UserMatrix()
  } else {
    stop(sQuote("CBASS"), " does not know how to handle ", sQuote("var_weights"),
         " of class ", class(var_weights)[1], ".")
  }

  if (any(var_weight_matrix < 0) || anyNA(var_weight_matrix)) {
    stop("All fusion weights for variables must be positive or zero.")
  }

  if (!is_connected_adj_mat(var_weight_matrix != 0)) {
    stop("Weights for variables do not imply a connected graph. Biclustering will not succeed.")
  }

  # Calculate observation (column)-clustering weights
  if (is.function(obs_weights)) { # Usual case, `obs_weights` is a function which calculates the weight matrix
    obs_weight_result <- obs_weights(t(X))

    if (is.matrix(obs_weight_result)) {
      obs_weight_matrix <- obs_weight_result
      obs_weight_type   <- UserFunction()
    } else {
      obs_weight_matrix <- obs_weight_result$weight_mat
      obs_weight_type   <- obs_weight_result$type
    }
  } else if (is.matrix(obs_weights)) {

    if (!is_square(obs_weights)) {
      stop(sQuote("obs_weights"), " must be a square matrix.")
    }

    if (NROW(obs_weights) != NCOL(X)) {
      stop(sQuote("NROW(obs_weights)"), " must be equal to ", sQuote("NCOL(X)."))
    }

    obs_weight_matrix <- obs_weights
    obs_weight_type   <- UserMatrix()
  } else {
    stop(sQuote("CBASS"), " does not know how to handle ", sQuote("obs_weights"),
         " of class ", class(obs_weights)[1], ".")
  }

  if (any(obs_weight_matrix < 0) || anyNA(obs_weight_matrix)) {
    stop("All fusion weights for observations must be positive or zero.")
  }

  if (!is_connected_adj_mat(obs_weight_matrix != 0)) {
    stop("Weights for observations do not imply a connected graph. Clustering will not succeed.")
  }

  ## NB: We are following Chi, Allen, and Baraniuk so "row" here refers to
  ##     features (variables) instead of the more typical observations
  ##     and vice versa for columns
  row_weights <- weight_mat_to_vec(var_weight_matrix)
  col_weights <- weight_mat_to_vec(obs_weight_matrix)

  ## Rescale to ensure coordinated fusions
  ##
  ## It is important that the observation (col) weights sum to 1/sqrt(p)
  ## and the feature/variable (row) weights sum to 1/sqrt(n).
  row_weights <- row_weights / (sum(row_weights) * sqrt(n.obs))
  col_weights <- col_weights / (sum(col_weights) * sqrt(p.var))

  PreCompList.row <- ConvexClusteringPreCompute(X = t(X),
                                                weights = row_weights,
                                                rho = rho,
                                                verbose = verbose.deep)
  cardE.row <- NROW(PreCompList.row$E)

  PreCompList.col <- ConvexClusteringPreCompute(X = X,
                                                weights = col_weights,
                                                rho = rho,
                                                verbose = verbose.deep)

  cardE.col <- NROW(PreCompList.col$E)

  if (verbose.basic) message("Computing CBASS Path")

  if (alg.type %in% c("cbassviz", "cbassvizl1")) {
    cbass.sol.path <- CBASS_VIZcpp(x = X[TRUE],
                                   n = as.integer(n.obs),
                                   p = as.integer(p.var),
                                   lambda_init = 1e-6,
                                   weights_row = row_weights[row_weights != 0],
                                   weights_col = col_weights[col_weights != 0],
                                   uinit_row = as.matrix(PreCompList.row$uinit),
                                   uinit_col = as.matrix(PreCompList.col$uinit),
                                   vinit_row = as.matrix(PreCompList.row$vinit),
                                   vinit_col = as.matrix(PreCompList.col$vinit),
                                   premat_row = PreCompList.row$PreMat,
                                   premat_col = PreCompList.col$PreMat,
                                   IndMat_row = PreCompList.row$ind.mat,
                                   IndMat_col = PreCompList.col$ind.mat,
                                   EOneIndMat_row = PreCompList.row$E1.ind.mat,
                                   EOneIndMat_col = PreCompList.col$E1.ind.mat,
                                   ETwoIndMat_row = PreCompList.row$E2.ind.mat,
                                   ETwoIndMat_col = PreCompList.col$E2.ind.mat,
                                   rho = rho,
                                   max_iter = as.integer(max.iter),
                                   burn_in = burn.in,
                                   verbose = verbose.deep,
                                   ti = 10,
                                   t_switch = 1.01,
                                   keep = 10,
                                   l1 = (alg.type == "cbassvizl1"))
  } else {
    cbass.sol.path <- CBASScpp(x = X[TRUE],
                               n = as.integer(n.obs),
                               p = as.integer(p.var),
                               lambda_init = 1e-6,
                               t = t,
                               weights_row = row_weights[row_weights != 0],
                               weights_col = col_weights[col_weights != 0],
                               uinit_row = as.matrix(PreCompList.row$uinit),
                               uinit_col = as.matrix(PreCompList.col$uinit),
                               vinit_row = as.matrix(PreCompList.row$vinit),
                               vinit_col = as.matrix(PreCompList.col$vinit),
                               premat_row = PreCompList.row$PreMat,
                               premat_col = PreCompList.col$PreMat,
                               IndMat_row = PreCompList.row$ind.mat,
                               IndMat_col = PreCompList.col$ind.mat,
                               EOneIndMat_row = PreCompList.row$E1.ind.mat,
                               EOneIndMat_col = PreCompList.col$E1.ind.mat,
                               ETwoIndMat_row = PreCompList.row$E2.ind.mat,
                               ETwoIndMat_col = PreCompList.col$E2.ind.mat,
                               rho = rho,
                               max_iter = as.integer(max.iter),
                               burn_in = burn.in,
                               verbose = verbose.deep,
                               keep = 10,
                               l1 = (alg.type == "cbassl1"))
  }

  ## FIXME - Convert lambda.path to a single column matrix instead of a vector
  ##         RcppArmadillo returns a arma::vec as a n-by-1 matrix
  ##         RcppEigen returns an Eigen::VectorXd as a n-length vector
  ##         Something downstream cares about the difference, so just change
  ##         the type here for now
  cbass.sol.path$lambda.path <- matrix(cbass.sol.path$lambda.path, ncol=1)

  if (verbose.basic) message("Post-processing")

  ISP(
    sp.path = cbass.sol.path$v.row.zero.inds %>% t(),
    v.path = cbass.sol.path$v.row.path,
    u.path = cbass.sol.path$u.path,
    lambda.path = cbass.sol.path$lambda.path,
    cardE = sum(row_weights != 0)
  ) -> cbass.cluster.path.row

  clust.path.row <- get_cluster_assignments(PreCompList.row$E, cbass.cluster.path.row$sp.path.inter, p.var)
  clust.path.dups.row <- duplicated(clust.path.row, fromLast = FALSE)

  cbass.cluster.path.row[["clust.path"]] <- clust.path.row
  cbass.cluster.path.row[["clust.path.dups"]] <- clust.path.dups.row

  cbass.dend.row <- CreateDendrogram(cbass.cluster.path.row, var_labels)

  ISP(
    sp.path = cbass.sol.path$v.col.zero.inds %>% t(),
    v.path = cbass.sol.path$v.col.path,
    u.path = cbass.sol.path$u.path,
    lambda.path = cbass.sol.path$lambda.path,
    cardE = sum(col_weights != 0)
  ) -> cbass.cluster.path.col

  clust.path.col <- get_cluster_assignments(PreCompList.col$E, cbass.cluster.path.col$sp.path.inter, n.obs)
  clust.path.dups.col <- duplicated(clust.path.col, fromLast = FALSE)

  cbass.cluster.path.col[["clust.path"]] <- clust.path.col
  cbass.cluster.path.col[["clust.path.dups"]] <- clust.path.dups.col

  cbass.dend.col <- CreateDendrogram(cbass.cluster.path.col, obs_labels)

  cbass.fit <- list(
    X = X.orig,
    cbass.sol.path = cbass.sol.path,
    cbass.cluster.path.obs = cbass.cluster.path.col,
    cbass.cluster.path.var = cbass.cluster.path.row,
    cbass.dend.var = cbass.dend.row,
    cbass.dend.obs = cbass.dend.col,
    n.obs = n.obs,
    p.var = p.var,
    var_weight_type = var_weight_type,
    obs_weight_type = obs_weight_type,
    burn.in = burn.in,
    alg.type = alg.type,
    t = t,
    X.center.global = X.center.global,
    obs.labels = obs_labels,
    var.labels = var_labels
  )

  class(cbass.fit) <- "CBASS"

  return(cbass.fit)
}

#' Print \code{CBASS} Results
#'
#' Prints a brief descripton of a fitted \code{CBASS} object.
#'
#' Reports number of observations and variables of dataset, any preprocessing
#' done by the \code{\link{CBASS}} function, regularization weight information,
#' and the variant of \code{CBASS}.
#'
#' @param x an object of class \code{CBASS} as returned by \code{\link{CBASS}}
#' @param ... Additional unused arguments
#' @export
#' @examples
#' cbass_fit <- CBASS(X=presidential_speech[1:10,1:4])
#' print(cbass_fit)
print.CBASS <- function(x, ...) {
  alg_string <- switch(x$alg.type,
                       cbass      = paste0("CBASS (t = ", round(x$t, 3), ")"),
                       cbassl1    = paste0("CBASS (t = ", round(x$t, 3), ") [L1]"),
                       cbassviz   = "CBASS-VIZ",
                       cbassvizl1 = "CBASS-VIZ [L1]")

  cat("CBASS Fit Summary\n")
  cat("====================\n\n")
  cat("Algorithm: ", alg_string, "\n\n")

  cat("Number of Observations: ", x$n.obs, "\n")
  cat("Number of Variables:    ", x$p.var, "\n\n")

  cat("Pre-processing options:\n")
  cat(" - Global centering: ", x$X.center.global, "\n\n")

  cat("Observation Weights:\n")
  print(x$obs_weight_type)

  cat("Feature Weights:\n")
  print(x$var_weight_type)

  cat("Raw Data:\n")
  print(x$X[1:min(5, x$n.obs), 1:min(5, x$p.var)])

  invisible(x)
}
