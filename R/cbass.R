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
#' @param control A list containing advanced parameters for the \code{CBASS} algorithm,
#'                typically created by \code{\link{cbass.control}}.
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
#' @param ... Additional arguments used to control the behavior of \code{CBASS}; see
#'            \code{\link{cbass.control}} for details.
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
#' @importFrom stats var
#' @export
#' @examples
#' \dontrun{
#' cbass_fit <- CBASS(presidential_speech)
#' print(cbass_fit)
#' plot(cbass_fit)
#' }
CBASS <- function(X,
                  verbose = 1L,
                  var_weights = sparse_gaussian_kernel_weights(k = "auto",
                                                               phi = "auto",
                                                               dist.method = "euclidean",
                                                               p = 2),
                  obs_weights = sparse_gaussian_kernel_weights(k = "auto",
                                                               phi = "auto",
                                                               dist.method = "euclidean",
                                                               p = 2),
                  ...,
                  control = NULL) {

  if (!is.matrix(X)) {
    warning(sQuote("X"), " should be a matrix, not a " , class(X)[1],
            ". Converting with as.matrix().")
    X <- as.matrix(X)
  }

  if (!is.numeric(X)) {
    stop(sQuote("X"), " must be numeric.")
  }

  n.obs <- NROW(X)
  p.var <- NCOL(X)

  if (anyNA(X)) {
    stop(sQuote("CBASS"), " cannot handle missing data.")
  }

  if (!all(is.finite(X))) {
    stop("All elements of ", sQuote("X"), " must be finite.")
  }

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

  internal.control <- cbass.control(...)
  if (!is.null(control)) {
    internal.control <- modifyList(internal.control, control)
  }

  obs.labels <- internal.control$obs.labels
  var.labels <- internal.control$var.labels
  X.center.global <- internal.control$X.center.global
  rho <- internal.control$rho
  max.iter <- internal.control$max.iter
  burn.in <- internal.control$burn.in
  alg.type <- internal.control$alg.type
  t <- internal.control$t
  npcs <- internal.control$npcs

  # get labels
  if (is.null(obs.labels)) {
    if (!is.null(rownames(X))) {
      n.labels <- rownames(X)
    } else {
      n.labels <- 1:NROW(X)
    }
  } else {
    if (length(obs.labels) == n.obs) {
      n.labels <- obs.labels
    } else {
      stop("obs.labels should be length NROW(X)")
    }
  }
  if (is.null(var.labels)) {
    if (!is.null(colnames(X))) {
      p.labels <- colnames(X)
    } else {
      p.labels <- 1:NCOL(X)
    }
  } else {
    if (length(var.labels) == p.var) {
      p.labels <- var.labels
    } else {
      stop("var.labels should be length NCOL(X)")
    }
  }

  # center and scale
  colnames(X) <- var.labels
  rownames(X) <- obs.labels
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
    stop(sQuote("CARP"), " does not know how to handle ", sQuote("var_weights"),
         " of class ", class(var_weights)[1], ".")
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
    stop(sQuote("CARP"), " does not know how to handle ", sQuote("obs_weights"),
         " of class ", class(obs_weights)[1], ".")
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

  PreCompList.row <- suppressMessages(
    ConvexClusteringPreCompute(
      X = t(X),
      weights = row_weights,
      rho = rho
    )
  )
  cardE.row <- NROW(PreCompList.row$E)

  PreCompList.col <- suppressMessages(
    ConvexClusteringPreCompute(
      X = X,
      weights = col_weights,
      rho = rho
    )
  )
  cardE.col <- NROW(PreCompList.col$E)

  if (verbose.basic) message("Computing CBASS Path")

  if (alg.type %in% c("cbassviz", "cbassvizl1")) {
    bicarp.sol.path <- CBASS_VIZcpp(x = X[TRUE],
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
    bicarp.sol.path <- CBASScpp(x = X[TRUE],
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
  bicarp.sol.path$lambda.path <- matrix(bicarp.sol.path$lambda.path, ncol=1)

  if (verbose.basic) message("Post-processing")

  ISP(
    sp.path = bicarp.sol.path$v.row.zero.inds %>% t(),
    v.path = bicarp.sol.path$v.row.path,
    u.path = bicarp.sol.path$u.path,
    lambda.path = bicarp.sol.path$lambda.path,
    cardE = sum(row_weights != 0)
  ) -> bicarp.cluster.path.row

  clust.path.row <- get_cluster_assignments(PreCompList.row$E, bicarp.cluster.path.row$sp.path.inter, p.var)
  clust.path.dups.row <- duplicated(clust.path.row, fromLast = FALSE)

  bicarp.cluster.path.row[["clust.path"]] <- clust.path.row
  bicarp.cluster.path.row[["clust.path.dups"]] <- clust.path.dups.row

  bicarp.dend.row <- CreateDendrogram(bicarp.cluster.path.row, p.labels)

  ISP(
    sp.path = bicarp.sol.path$v.col.zero.inds %>% t(),
    v.path = bicarp.sol.path$v.col.path,
    u.path = bicarp.sol.path$u.path,
    lambda.path = bicarp.sol.path$lambda.path,
    cardE = sum(col_weights != 0)
  ) -> bicarp.cluster.path.col

  clust.path.col <- get_cluster_assignments(PreCompList.col$E, bicarp.cluster.path.col$sp.path.inter, n.obs)
  clust.path.dups.col <- duplicated(clust.path.col, fromLast = FALSE)

  bicarp.cluster.path.col[["clust.path"]] <- clust.path.col
  bicarp.cluster.path.col[["clust.path.dups"]] <- clust.path.dups.col

  bicarp.dend.col <- CreateDendrogram(bicarp.cluster.path.col, n.labels)

  cbass.fit <- list(
    X = X.orig,
    cbass.sol.path = bicarp.sol.path,
    cbass.cluster.path.obs = bicarp.cluster.path.col,
    cbass.cluster.path.var = bicarp.cluster.path.row,
    cbass.dend.var = bicarp.dend.row,
    cbass.dend.obs = bicarp.dend.col,
    n.obs = n.obs,
    p.var = p.var,
    var_weight_type = var_weight_type,
    obs_weight_type = obs_weight_type,
    burn.in = burn.in,
    alg.type = alg.type,
    t = t,
    X.center.global = X.center.global,
    obs.labels = n.labels,
    var.labels = p.labels
  )

  class(cbass.fit) <- "CBASS"

  return(cbass.fit)
}

#' Control for \code{CBASS} fits
#'
#' Set \code{CBASS} algorithm parameters
#'
#' This function constructs a list containing additional arguments to control
#' the behavior of the \code{CBASS} algorithm. It is typically only used internally
#' by \code{\link{CBASS}}, but may be useful to advanced users who wish to
#' construct the \code{control} argument directly.
#'
#' @param obs.labels A character vector of length \eqn{n}: observations (row) labels
#' @param var.labels A character vector of length \eqn{p}: variable (column) labels
#' @param X.center.global A logical: Should \code{X} be centered globally?
#'                        \emph{I.e.}, should the global mean of \code{X} be subtracted?
#' @param rho For advanced users only (not advisable to change): the penalty
#'            parameter used for the augmented Lagrangian.
#' @param max.iter An integer: the maximum number of CARP iterations.
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
#' @return A list containing the \code{CBASS} algorithm parameters.
#' @export
cbass.control <- function(obs.labels = NULL,
                          var.labels = NULL,
                          X.center.global = TRUE,
                          rho = 1,
                          t = 1.01,
                          max.iter = as.integer(1e6),
                          burn.in = as.integer(50),
                          alg.type = "cbassviz",
                          npcs = as.integer(4),
                          ...) {

  dots <- list(...)

  if (length(dots) != 0L) {
    if (!is.null(names(dots))) {
      stop("Unknown argument ", sQuote(names(dots)[1L]), " passed to ", sQuote("CARP."))
    } else {
      stop("Unknown ", sQuote("..."), " arguments passed to ", sQuote("CARP."))
    }
  }

  if (!is.logical(X.center.global) || is.na(X.center.global) || (length(X.center.global) != 1L)) {
    stop(sQuote("X.center.global"), "must be either ", sQuote("TRUE"), " or ", sQuote("FALSE."))
  }

  if ((rho < 0) || is.na(rho) || (length(rho) != 1L)) {
    stop(sQuote("rho"), "must a be non-negative scalar.")
  }

  if (!is.null(npcs)) {
    if (!is.integer(npcs) || npcs <= 1L) {
      stop(sQuote("npcs"), " must be at least 2.")
    }
  }

  if (!is.integer(max.iter) || (max.iter <= 0) || (length(max.iter) != 1L)) {
    stop(sQuote("max.iter"), " must be a positive integer.")
  }

  if (!is.integer(burn.in) || (burn.in <= 0) || (burn.in >= max.iter)) {
    stop(sQuote("burn.in"), " must be a positive integer less than ", sQuote("max.iter."))
  }

  if (alg.type %not.in% c("cbassviz", "cbass", "cbassl1", "cbassvizl1")) {
    stop("Unrecognized value of ", sQuote("alg.type;"), " see help for allowed values.")
  }

  if ((t <= 1) || is.na(t) || (length(t) != 1L)) {
    stop(sQuote("t"), " must be a scalar greater than 1.")
  }

  list(
    obs.labels = obs.labels,
    var.labels = var.labels,
    X.center.global = X.center.global,
    rho = rho,
    t = t,
    max.iter = max.iter,
    burn.in = burn.in,
    alg.type = alg.type
  )
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
