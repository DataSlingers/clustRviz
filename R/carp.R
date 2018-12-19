#' Compute \code{CARP} (Convex Clustering) Solution Path
#'
#' \code{CARP} returns a fast approximation to the Convex Clustering
#' solution path along with visualizations such as dendrograms and
#' cluster paths. \code{CARP} solves the Convex Clustering problem via an efficient
#' Algorithmic Regularization scheme.
#'
#' @param X The data matrix (\eqn{X \in R^{n \times p}}{X}): rows correspond to
#'          the observations (to be clustered) and columns to the variables (which
#'          will not be clustered).
#' @param labels A character vector of length \eqn{n}: observations (row) labels
#' @param X.center A logical: Should \code{X} be centered columnwise?
#' @param X.scale A logical: Should \code{X} be scaled columnwise?
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
#'            arguments as given. All arguments other than \code{X} must be given
#'            by name.
#' @param weights One of the following: \itemize{
#'                \item A function which, when called with argument \code{X},
#'                      returns an b-by-n matrix of fusion weights.
#'                \item A matrix of size n-by-n containing fusion weights
#'                }
#' @param status Should a status message be printed to the console?
#' @return An object of class \code{CARP} containing the following elements (among others):
#'         \itemize{
#'         \item \code{X}: the original data matrix
#'         \item \code{n}: the number of observations (rows of \code{X})
#'         \item \code{p}: the number of variables (columns of \code{X})
#'         \item \code{alg.type}: the \code{CARP} variant used
#'         \item \code{X.center}: a logical indicating whether \code{X} was centered
#'                                column-wise before clustering
#'         \item \code{X.scale}: a logical indicating whether \code{X} was scaled
#'                               column-wise before centering
#'         \item \code{weight_type}: a record of the scheme used to create
#'                                   fusion weights
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
                 back_track = FALSE,
                 exact = FALSE,
                 norm = 2,
                 t = 1.05,
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

  n <- NROW(X)
  p <- NCOL(X)

  # Center and scale X
  X.orig <- X
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
  D <- matrix(0, ncol = n, nrow = cardE)
  D[cbind(seq_len(cardE), edge_list[,1])] <-  1
  D[cbind(seq_len(cardE), edge_list[,2])] <- -1

  weight_vec <- weight_mat_to_vec(weight_matrix)

  crv_message("Computing Convex Clustering [CARP] Path")
  tic_inner <- Sys.time()

  carp.sol.path <- CARPcpp(X,
                           D,
                           t = t,
                           epsilon = .clustRvizOptionsEnv[["epsilon"]],
                           weights = weight_vec[weight_vec != 0],
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
  carp.sol.path$gamma_path <- matrix(carp.sol.path$gamma_path, ncol=1)

  crv_message("Post-processing")

  post_processing_results <- ConvexClusteringPostProcess(X = X,
                                                         edge_matrix      = edge_list,
                                                         gamma_path       = carp.sol.path$gamma_path,
                                                         u_path           = carp.sol.path$u_path,
                                                         v_path           = carp.sol.path$v_path,
                                                         v_zero_indices   = carp.sol.path$v_zero_inds,
                                                         labels           = labels,
                                                         dendrogram_scale = dendrogram.scale,
                                                         npcs             = npcs,
                                                         smooth_U         = TRUE)

  carp.fit <- list(
    X = X.orig,
    D = D,
    U = post_processing_results$U,
    dendrogram = post_processing_results$dendrogram,
    rotation_matrix = post_processing_results$rotation_matrix,
    cluster_membership = post_processing_results$membership_info,
    n = n,
    p = p,
    weights = weight_matrix,
    weight_type = weight_type,
    back_track = back_track,
    exact = exact,
    norm = norm,
    t = t,
    X.center = X.center,
    center_vector = center_vector,
    X.scale = X.scale,
    scale_vector = scale_vector,
    time = Sys.time() - tic,
    fit_time = toc_inner - tic_inner
  )

  if (.clustRvizOptionsEnv[["keep_debug_info"]]) {
    carp.fit[["debug"]] <- list(path = carp.sol.path,
                                row  = post_processing_results$debug)
  }

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
#' @details The \code{as.dendrogram} and \code{as.hclust} methods convert the
#'          \code{CBASS} output to an object of class \code{dendrogram} or \code{hclust}
#'          respectively.
#'
#' @param x an object of class \code{CARP} as returned by \code{\link{CARP}}
#' @param object an object of class \code{CARP} as returned by \code{\link{CARP}}
#' @param ... Additional unused arguments
#' @export
#' @rdname print_carp
#' @examples
#' carp_fit <- CARP(presidential_speech)
#' print(carp_fit)
print.CARP <- function(x, ...) {
  if(x$exact){
    if(x$back_track){
      alg_string = "ADMM-VIZ [Exact Solver + Back-Tracking Fusion Search]"
    } else {
      alg_string = paste0("ADMM (t = ", round(x$t, 3), ") [Exact Solver]")
    }
  } else {
    if(x$back_track){
      alg_string = "CARP-VIZ [Back-Tracking Fusion Search]"
    } else {
      alg_string = paste0("CARP (t = ", round(x$t, 3), ")")
    }
  }

  if(x$norm == 1){
    alg_string <- paste(alg_string, "[L1]")
  }

  cat("CARP Fit Summary\n")
  cat("====================\n\n")
  cat("Algorithm:", alg_string, "\n")
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

#' @export
#' @importFrom stats as.dendrogram
#' @rdname print_carp
as.dendrogram.CARP <- function(object, ...){
  as.dendrogram(object$dendrogram)
}

#' @export
#' @importFrom stats as.hclust
#' @rdname print_carp
as.hclust.CARP <- function(x, ...){
  x$dendrogram
}
