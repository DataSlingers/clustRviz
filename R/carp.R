#' Compute \code{CARP} (Convex Clustering) Solution Path
#'
#' \code{CARP} returns a fast approximation to the Convex Clustering
#' solution path along with visualizations such as dendrograms and
#' cluster paths. \code{CARP} solves the Convex Clustering problem via
#' Algorithmic Regularization Paths. A seqeunce of clustering
#' solutions is returned along with several visualizations.
#'
#' @param X The data matrix (\eqn{X \in R^{n \times p}}{X}): rows correspond to
#'          the observations (to be clustered) and columns to the variables (which
#'          will not be clustered).
#' @param labels A character vector of length \eqn{n}: observations (row) labels
#' @param X.center A logical: Should \code{X} be centered columnwise?
#' @param X.scale A logical: Should \code{X} be scaled columnwise?
#' @param rho For advanced users only (not advisable to change): the penalty
#'            parameter used for the augmented Lagrangian.
#' @param max.iter An integer: the maximum number of \code{CARP} iterations.
#' @param burn.in An integer: the number of initial iterations at a fixed
#'                (small) value of \eqn{\lambda}
#' @param alg.type Which \code{CARP} variant to use. Allowed values are \itemize{
#'        \item \code{"carp"} - The standard \code{CARP} algorithm with \eqn{L2} penalty;
#'        \item \code{"carpviz"} - The back-tracking \code{CARP} algorithm with \eqn{L2} penalty;
#'        \item \code{"carpl1"} - The standard \code{CARP} algorithm with \eqn{L1} penalty; and
#'        \item \code{"carpvizl1"} - The back-tracking \code{CARP} algorithm with \eqn{L1} penalty.}
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
#'            arguments as given. All arguments other than \code{X} must be given
#'            by name.
#' @param weights One of the following: \itemize{
#'                \item A function which, when called with argument \code{X},
#'                      returns an b-by-n matrix of fusion weights.
#'                \item A matrix of size n-by-n containing fusion weights
#'                }
#' @return An object of class \code{CARP} containing the following elements (among others):
#'         \itemize{
#'         \item \code{X}: the original data matrix
#'         \item \code{n.obs}: the number of observations (rows of \code{X})
#'         \item \code{p.var}: the number of variables (columns of \code{X})
#'         \item \code{alg.type}: the \code{CARP} variant used
#'         \item \code{X.center}: a logical indicating whether \code{X} was centered
#'                                column-wise before clustering
#'         \item \code{X.scale}: a logical indicating whether \code{X} was scaled
#'                               column-wise before centering
#'         \item \code{burn.in}: an integer indicating the number of "burn-in"
#'                               iterations performed
#'         \item \code{weight_type}: a record of the scheme used to create
#'                                   fusion weights
#'         \item \code{carp.dend}: a dendrogram (object of class
#'                                 \code{\link[stats]{hclust}}) containing
#'                                 the clustering solution path
#'         \item \code{carp.cluster.path.vis}: The \code{CARP} solution path
#'         }
#' @importFrom utils data
#' @importFrom dplyr %>% mutate group_by ungroup as_tibble n_distinct
#' @importFrom rlang %||%
#' @importFrom stats var
#' @export
#' @examples
#' carp_fit <- CARP(presidential_speech[1:10,1:4])
#' print(carp_fit)
#' plot(carp_fit)
CARP <- function(X,
                 ...,
                 weights = sparse_rbf_kernel_weights(k = "auto",
                                                     phi = "auto",
                                                     dist.method = "euclidean",
                                                     p = 2),
                 labels = rownames(X),
                 X.center = TRUE,
                 X.scale = FALSE,
                 rho = 1.0,
                 max.iter = 1000000L,
                 burn.in = 50L,
                 alg.type = c("carpviz", "carpvizl1", "carp", "carpl1"),
                 t = 1.05,
                 npcs = min(4L, NCOL(X), NROW(X)),
                 dendrogram.scale = NULL) {

  tic <- Sys.time()

  ####################
  ##
  ## Input validation
  ##
  ####################

  dots <- list(...)

  if (length(dots) != 0L) {
    if (!is.null(names(dots))) {
      crv_error("Unknown argument ", sQuote(names(dots)[1L]), " passed to ", sQuote("CARP."))
    } else {
      crv_error("Unknown ", sQuote("..."), " arguments passed to ", sQuote("CARP."))
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
    crv_error(sQuote("CARP"), " cannot handle missing data.")
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

  if ( (!is_numeric_scalar(rho)) || (rho <= 0)) {
    crv_error(sQuote("rho"), "must be a positive scalar (vector of length 1).")
  }

  if ( (!is_integer_scalar(max.iter)) || (max.iter <= 1L) ) {
    crv_error(sQuote("max.iter"), " must be a positive integer scalar and at least 2.")
  }

  if ( (!is_integer_scalar(burn.in)) || (burn.in <= 0L) || (burn.in >= max.iter) ) {
    crv_error(sQuote("burn.in"), " must be a positive integer less than ", sQuote("max.iter."))
  }

  alg.type <- match.arg(alg.type)

  if ( (!is_numeric_scalar(t)) || (t <= 1) ) {
    crv_error(sQuote("t"), " must be a scalar greater than 1.")
  }

  if (!is.null(dendrogram.scale)) {
    if (dendrogram.scale %not.in% c("original", "log")) {
      crv_error("If not NULL, ", sQuote("dendrogram.scale"), " must be either ", sQuote("original"), " or ", sQuote("log."))
    }
  }

  if ( (!is_integer_scalar(npcs)) || (npcs < 2) || (npcs > NCOL(X)) || (npcs > NROW(X)) ){
    crv_error(sQuote("npcs"), " must be an integer scalar between 2 and ", sQuote("min(dim(X))."))
  }

  ## Get row (observation) labels
  if (is.null(labels)) {
    labels <- paste0("Obs", seq_len(NROW(X)))
  }

  if ( length(labels) != NROW(X) ){
    crv_error(sQuote("labels"), " must be of length ", sQuote("NROW(X)."))
  }

  rownames(X) <- labels <- make.unique(as.character(labels), sep="_")

  n.obs <- NROW(X)
  p.var <- NCOL(X)

  # Center and scale X
  X.orig <- X
  if (X.center | X.scale) {
    X <- scale(X, center = X.center, scale = X.scale)
  }

  scale_vector  <- attr(X, "scaled:scale", exact=TRUE)  %||% rep(1, p.var)
  center_vector <- attr(X, "scaled:center", exact=TRUE) %||% rep(0, p.var)

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
    crv_error(sQuote("CARP"), " does not know how to handle ", sQuote("weights"),
              " of class ", class(weights)[1], ".")
  }

  if (any(weight_matrix < 0) || anyNA(weight_matrix)) {
    crv_error("All fusion weights must be positive or zero.")
  }

  if (!is_connected_adj_mat(weight_matrix != 0)) {
    crv_error("Weights do not imply a connected graph. Clustering will not succeed.")
  }

  weight_matrix_ut <- weight_matrix * upper.tri(weight_matrix);

  edge_list <- which(weight_matrix_ut != 0, arr.ind = TRUE)
  edge_list <- edge_list[order(edge_list[, 1], edge_list[, 2]), ]
  cardE <- NROW(edge_list)
  D <- matrix(0, ncol = n.obs, nrow = cardE)
  D[cbind(seq_len(cardE), edge_list[,1])] <-  1
  D[cbind(seq_len(cardE), edge_list[,2])] <- -1

  weight_vec <- weight_mat_to_vec(weight_matrix)

  crv_message("Computing CARP Path")

  if (alg.type %in% c("carpvizl1", "carpviz")) {
      carp.sol.path <- CARP_VIZcpp(X,
                                   D,
                                   lambda_init = 1e-8,
                                   weights = weight_vec[weight_vec != 0],
                                   rho = rho,
                                   max_iter = as.integer(max.iter),
                                   burn_in = as.integer(burn.in),
                                   ti = 10,
                                   t_switch = 1.01,
                                   keep = 1,
                                   l1 = (alg.type == "carpvizl1"))
  } else {
      carp.sol.path <- CARPcpp(X,
                               D,
                               lambda_init = 1e-8,
                               t = t,
                               weights = weight_vec[weight_vec != 0],
                               rho = rho,
                               max_iter = as.integer(max.iter),
                               burn_in = as.integer(burn.in),
                               keep = 1,
                               l1 = (alg.type == "carpl1"))
  }

  ## FIXME - Convert lambda.path to a single column matrix instead of a vector
  ##         RcppArmadillo returns a arma::vec as a n-by-1 matrix
  ##         RcppEigen returns an Eigen::VectorXd as a n-length vector
  ##         Something downstream cares about the difference, so just change
  ##         the type here for now
  carp.sol.path$lambda.path <- matrix(carp.sol.path$lambda.path, ncol=1)

  crv_message("Post-processing")

  post_processing_results <- ConvexClusteringPostProcess(X = X,
                                                         edge_matrix      = edge_list,
                                                         lambda_path      = carp.sol.path$lambda.path,
                                                         u_path           = carp.sol.path$u.path,
                                                         v_path           = carp.sol.path$v.path,
                                                         v_zero_indices   = carp.sol.path$v.zero.inds,
                                                         labels           = labels,
                                                         dendrogram_scale = dendrogram.scale,
                                                         npcs             = npcs)

  carp.fit <- list(
    X = X.orig,
    carp.dend = post_processing_results$dendrogram,
    carp.cluster.path.vis = post_processing_results$paths,
    carp.sol.path = carp.sol.path,
    cardE = cardE,
    n.obs = n.obs,
    p.var = p.var,
    weight_type = weight_type,
    burn.in = burn.in,
    alg.type = alg.type,
    t = t,
    X.center = X.center,
    center_vector = center_vector,
    X.scale = X.scale,
    scale_vector = scale_vector,
    time = Sys.time() - tic
  )

  class(carp.fit) <- "CARP"

  return(carp.fit)
}

#' Print \code{CARP} Results
#'
#' Prints a brief descripton of a fitted \code{CARP} object.
#'
#' Reports number of observations and variables of dataset, any preprocessing
#' done by the \code{\link{CARP}} function, regularization weight information,
#' and the variant of \code{CARP} used.
#'
#' @param x an object of class \code{CARP} as returned by \code{\link{CARP}}
#' @param ... Additional unused arguments
#' @export
#' @examples
#' carp_fit <- CARP(presidential_speech[1:10,1:4])
#' print(carp_fit)
print.CARP <- function(x, ...) {
  alg_string <- switch(x$alg.type,
                       carp      = paste0("CARP (t = ", round(x$t, 3), ")"),
                       carpl1    = paste0("CARP (t = ", round(x$t, 3), ") [L1]"),
                       carpviz   = "CARP-VIZ",
                       carpvizl1 = "CARP-VIZ [L1]")

  cat("CARP Fit Summary\n")
  cat("====================\n\n")
  cat("Algorithm:", alg_string, "\n")
  cat("Time:", sprintf("%2.3f %s", x$time, attr(x$time, "units")), "\n\n")

  cat("Number of Observations:", x$n.obs, "\n")
  cat("Number of Variables:   ", x$p.var, "\n\n")

  cat("Pre-processing options:\n")
  cat(" - Columnwise centering:", x$X.center, "\n")
  cat(" - Columnwise scaling:  ", x$X.scale, "\n\n")

  cat("Weights:\n")
  print(x$weight_type)

  cat("Raw Data:\n")
  print(x$X[1:min(5, x$n.obs), 1:min(5, x$p.var)])

  invisible(x)
}
