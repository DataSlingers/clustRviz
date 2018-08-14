#' Compute \code{CBASS} (Convex BiClustering) Solution Path
#'
#' \code{CBASS} returns a fast approximation to the Convex BiClustering
#' solution path along with visualizations such as dendrograms and
#' heatmaps. Visualizations may be static, interactive, or both.
#'
#' \code{CBASS} solves the Convex Biclustering problem via
#' Algorithmic Regularization Paths. A seqeunce of biclustering
#' solutions is returned along with several visualizations.
#'
#' @param X The data matrix (\eqn{X \in R^{n \times p}}{X}): rows correspond to
#'          the observations (to be clustered) and columns to the variables (which
#'          will not be clustered).
#' @param verbose Any of the values \code{0}, \code{1}, or \code{2}. Higher values
#'                correspond to more verbose output while running.
#' @param interactive A logical. Should interactive paths and dendrograms be
#'                    returned?
#' @param static A logical. Should static paths and dendrograms be returned?
#' @param control A list containing advanced parameters for the \code{CBASS} algorithm,
#'                typically created by \code{\link{cbass.control}}.
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
#'         \item \code{X.scale}: a logical indicating whether \code{X} was scaled
#'                               column-wise before centering
#'         \item \code{burn.in}: an integer indicating the number of "burn-in"
#'                               iterations performed
#'         \item \code{k.obs}: the number of neighbors used to create sparse
#'                             clustering weights for the observations
#'         \item \code{k.var}: the number of neighbors used to create sparse
#'                             clustering weights for the variables
#'         \item \code{phi.obs}: the scale factor of the RBF kernel used to calculate
#'                               clustering weights for the observations
#'         \item \code{phi.var}: the scale factor of the RBF kernel used to calculate
#'                               clustering weights for the variables
#'         \item \code{carp.dend.obs}: If \code{static=TRUE}, an dendrogram (object of
#'                                     class \code{\link[stats]{hclust}}) containing
#'                                     the clustering solution path for the observations
#'         \item \code{carp.dend.var}: If \code{static=TRUE}, an dendrogram (object of
#'                                     class \code{\link[stats]{hclust}}) containing
#'                                     the clustering solution path for the variables
#'         \item \code{cbass.cluster.path.obs}: The \code{CBASS} solution path for the observations
#'         \item \code{cbass.cluster.path.var}: The \code{CBASS} solution path for the variables
#'         \item \code{static}: a logical indicating whether static visualizations are available for this \code{CBASS} object
#'         \item \code{interactive}: a logical indicating whether interactive visualizations are available for this \code{CBASS} object
#'         \item \code{obs.labels}: a character vector of length \code{n.obs} containing
#'                                  observation (row) labels
#'         \item \code{var.labels}: a character vector of length \code{p.var} containing
#'                                  variable (column) labels
#'         }
#' @return obs.labels a vector of length n.obs containing observations (row) labels
#' @return var.labels a vector of length p.var containing variable (column) labels
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
                  interactive = TRUE,
                  static = TRUE,
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
  k.obs <- internal.control$k.obs
  k.var <- internal.control$k.var
  phi <- internal.control$phi
  rho <- internal.control$rho
  weights.obs <- internal.control$weights.obs
  weights.var <- internal.control$weights.var
  obs.weight.dist <- internal.control$obs.weight.dist
  var.weight.dist <- internal.control$var.weight.dist
  obs.weight.dist.p <- internal.control$obs.weight.dist.p
  var.weight.dist.p <- internal.control$var.weight.dist.p
  ncores <- internal.control$ncores
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
  if (!is.null(phi)) {
    if (phi <= 0) {
      stop("phi should be positive.")
    }
  }


  # center and scale
  colnames(X) <- var.labels
  rownames(X) <- obs.labels
  X.orig <- X
  if (X.center.global) {
    X <- X - mean(X)
    X <- t(X)
  } else {
    X <- t(X)
  }

  if (!is.null(weights.var)) {
    if (length(weights.var) == choose(p.var, 2)) {
      weights.row <- weights.var
    } else {
      stop("weights.var should be length choose(NCOL(X),2)")
    }
  } else {
    if (is.null(phi)) {
      phi.vec <- 10^(-10:10)
      sapply(phi.vec, function(phi) {
        stats::var(DenseWeights(X = X, phi = phi, method = var.weight.dist, p = var.weight.dist.p))
      }) %>%
        which.max() %>%
        phi.vec[.] -> phi
    }
    phi.row <- phi / n.obs
    weights.row <- DenseWeights(X = X, phi = phi.row, method = var.weight.dist, p = var.weight.dist.p)
    if (is.null(k.var)) {
      k.row <- MinKNN(X = X, dense.weights = weights.row)
    } else {
      k.row <- k.var
    }
    weights.row <- SparseWeights(X = X, dense.weights = weights.row, k = k.row)
  }
  weights.row <- weights.row / sum(weights.row)
  weights.row <- weights.row / sqrt(n.obs)

  PreCompList.row <- suppressMessages(
    ConvexClusteringPreCompute(
      X = t(X),
      weights = weights.row,
      ncores = ncores, rho = rho
    )
  )
  cardE.row <- NROW(PreCompList.row$E)

  if (!is.null(weights.obs)) {
    if (length(weights.obs) == choose(n.obs, 2)) {
      weights.cols <- weights.obs
    } else {
      stop("weights.obs should have length choose(NROW(X),2)")
    }
  } else {
    phi.col <- phi / p.var
    weights.col <- DenseWeights(X = t(X), phi = phi.col, method = obs.weight.dist, p = obs.weight.dist.p)
    if (is.null(k.obs)) {
      k.col <- MinKNN(X = t(X), dense.weights = weights.col)
    } else {
      k.col <- k.obs
    }
    weights.col <- SparseWeights(X = t(X), dense.weights = weights.col, k = k.col)
  }
  weights.col <- weights.col / sum(weights.col)
  weights.col <- weights.col / sqrt(p.var)

  PreCompList.col <- suppressMessages(
    ConvexClusteringPreCompute(
      X = X,
      weights = weights.col,
      ncores = ncores, rho = rho
    )
  )
  cardE.col <- NROW(PreCompList.col$E)



  if (verbose.basic) message("Computing CBASS Path\n")
  switch(
    alg.type,
    cbassviz = {
      BICARPL2_VIS(
        x = X[TRUE],
        n = as.integer(n.obs),
        p = as.integer(p.var),
        lambda_init = 1e-6,
        weights_row = weights.row[weights.row != 0],
        weights_col = weights.col[weights.col != 0],
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
        verbose_inner = verbose.deep,
        try_tol = 1e-3,
        ti = 10,
        t_switch = 1.01
      ) -> bicarp.sol.path
    },
    cbassvizl1 = {
      BICARPL1_VIS(
        x = X[TRUE],
        n = as.integer(n.obs),
        p = as.integer(p.var),
        lambda_init = 1e-6,
        weights_row = weights.row[weights.row != 0],
        weights_col = weights.col[weights.col != 0],
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
        verbose_inner = verbose.deep,
        try_tol = 1e-3,
        ti = 10,
        t_switch = 1.01
      ) -> bicarp.sol.path
    },
    cbass = {
      BICARPL2_NF_FRAC(
        x = X[TRUE],
        n = as.integer(n.obs),
        p = as.integer(p.var),
        lambda_init = 1e-6,
        t = t,
        weights_row = weights.row[weights.row != 0],
        weights_col = weights.col[weights.col != 0],
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
        keep = 10
      ) -> bicarp.sol.path
    },
    cbassl1 = {
      BICARPL1_NF_FRAC(
        x = X[TRUE],
        n = as.integer(n.obs),
        p = as.integer(p.var),
        lambda_init = 1e-6,
        t = t,
        weights_row = weights.row[weights.row != 0],
        weights_col = weights.col[weights.col != 0],
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
        keep = 10
      ) -> bicarp.sol.path
    }
  )

  if (verbose.basic) message("Post-processing\n")
  ISP(
    sp.path = bicarp.sol.path$v.row.zero.inds %>% t(),
    v.path = bicarp.sol.path$v.row.path,
    u.path = bicarp.sol.path$u.path,
    lambda.path = bicarp.sol.path$lambda.path,
    cardE = sum(weights.row != 0)
  ) -> bicarp.cluster.path.row
  bicarp.cluster.path.row$sp.path.inter %>% duplicated(fromLast = FALSE) -> sp.path.dups.row
  adj.path.row <- CreateAdjacencyPath(PreCompList.row$E, sp.path = bicarp.cluster.path.row$sp.path.inter, p.var)
  clust.graph.path.row <- CreateClusterGraphPath(adj.path.row)
  clust.path.row <- GetClustersPath(clust.graph.path.row)
  clust.path.dups.row <- duplicated(clust.path.row, fromLast = FALSE)

  bicarp.cluster.path.row[["sp.path.dups"]] <- sp.path.dups.row
  bicarp.cluster.path.row[["adj.path"]] <- adj.path.row
  bicarp.cluster.path.row[["clust.graph.path"]] <- clust.graph.path.row
  bicarp.cluster.path.row[["clust.path"]] <- clust.path.row
  bicarp.cluster.path.row[["clust.path.dups"]] <- clust.path.dups.row

  bicarp.dend.row <- CreateDendrogram(bicarp.cluster.path.row, p.labels)



  ISP(
    sp.path = bicarp.sol.path$v.col.zero.inds %>% t(),
    v.path = bicarp.sol.path$v.col.path,
    u.path = bicarp.sol.path$u.path,
    lambda.path = bicarp.sol.path$lambda.path,
    cardE = sum(weights.col != 0)
  ) -> bicarp.cluster.path.col
  bicarp.cluster.path.col$sp.path.inter %>% duplicated(fromLast = FALSE) -> sp.path.dups.col
  adj.path.col <- CreateAdjacencyPath(PreCompList.col$E, sp.path = bicarp.cluster.path.col$sp.path.inter, n.obs)
  clust.graph.path.col <- CreateClusterGraphPath(adj.path.col)
  clust.path.col <- GetClustersPath(clust.graph.path.col)
  clust.path.dups.col <- duplicated(clust.path.col, fromLast = FALSE)

  bicarp.cluster.path.col[["sp.path.dups"]] <- sp.path.dups.col
  bicarp.cluster.path.col[["adj.path"]] <- adj.path.col
  bicarp.cluster.path.col[["clust.graph.path"]] <- clust.graph.path.col
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
    phi.var = phi.row,
    phi.obs = phi.col,
    k.obs = k.col,
    k.var = k.row,
    burn.in = burn.in,
    alg.type = alg.type,
    t = t,
    X.center.global = X.center.global,
    static = static,
    interactive = interactive,
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
#' @param phi A positive real number: the scale factor used in the RBF kernel

#' @param weights.obs A vector of positive number of length \code{choose(n,2)}.
#' @param k.obs An positive integer: the number of neighbors used to create sparse weights
#' @param obs.weight.dist A string indicating the distance metric used to calculate weights.
#'                        See \code{\link[stats]{distance}} for details.
#' @param obs.weight.dist.p The exponent used to calculate the Minkowski distance if
#'                          \code{weight.dist = "minkowski"}.
#'                          See \code{\link[stats]{distance}} for details.
#' @param weights.var A vector of positive number of length \code{choose(n,2)}.
#' @param k.var An positive integer: the number of neighbors used to create sparse weights
#' @param var.weight.dist A string indicating the distance metric used to calculate weights.
#'                        See \code{\link[stats]{distance}} for details.
#' @param var.weight.dist.p The exponent used to calculate the Minkowski distance if
#'                          \code{weight.dist = "minkowski"}.
#'                          See \code{\link[stats]{distance}} for details.
#' @param ncores An positive integer: the number of cores to use.
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
                          phi = NULL,
                          weights.obs = NULL,
                          weights.var = NULL,
                          obs.weight.dist = "euclidean",
                          obs.weight.dist.p = 2,
                          var.weight.dist = "euclidean",
                          var.weight.dist.p = 2,
                          k.obs = NULL,
                          k.var = NULL,
                          t = 1.01,
                          ncores = as.integer(1),
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

  if (obs.weight.dist %not.in% SUPPORTED_DISTANCES) {
    stop("Unsupported choice of ",
         sQuote("obs.weight.dist;"),
         " see the ", sQuote("method"),
         " argument of ",
         sQuote("stats::dist"),
         " for supported distances.")
  }

  if ((obs.weight.dist.p <= 0) || (length(obs.weight.dist.p) != 1L)) {
    stop(sQuote("obs.weight.dist.p"),
         " must be a positive scalar; see the ", sQuote("p"),
         " argument of ", sQuote("stats::dist"), " for details.")
  }

  if (var.weight.dist %not.in% SUPPORTED_DISTANCES) {
    stop("Unsupported choice of ",
         sQuote("var.weight.dist;"),
         " see the ", sQuote("method"),
         " argument of ",
         sQuote("stats::dist"),
         " for supported distances.")
  }

  if ((var.weight.dist.p <= 0) || (length(var.weight.dist.p) != 1L)) {
    stop(sQuote("var.weight.dist.p"),
         " must be a positive scalar; see the ", sQuote("p"),
         " argument of ", sQuote("stats::dist"), " for details.")
  }

  if (!is.null(k.obs)) {
    if (!is.integer(k.obs) || k.obs <= 0) {
      stop("If not NULL, ", sQuote("k.obs"), " must be a positive integer.")
    }
  }

  if (!is.null(k.var)) {
    if (!is.integer(k.var) || k.obs <= 0) {
      stop("If not NULL, ", sQuote("k.var"), " must be a positive integer.")
    }
  }

  if (!is.integer(ncores) || ncores <= 0L) {
    stop(sQuote("ncores"), " must be a positive integer.")
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
    phi = phi,
    weights.obs = weights.obs,
    weights.var = weights.var,
    obs.weight.dist = obs.weight.dist,
    obs.weight.dist.p = obs.weight.dist.p,
    var.weight.dist = var.weight.dist,
    var.weight.dist.p = var.weight.dist.p,
    k.obs = k.obs,
    k.var = k.var,
    t = t,
    ncores = ncores,
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
#' the variant of \code{CBASS} used, and the visualizations returned.
#' @param x an object of class \code{CBASS} as returned by \code{\link{CBASS}}
#' @param ... Additional unused arguments
#' @export
#' @examples
#' cbass_fit <- CBASS(X=presidential_speech[1:10,1:4])
#' print(cbass_fit)
print.CBASS <- function(x, ...) {
  alg_string <- switch(x$alg.type,
                       carp      = paste0("CBASS (t = ", round(x$t, 3), ")"),
                       carpl1    = paste0("CBASS(t = ", round(x$t, 3), ") [L1]"),
                       carlviz   = "CBASS-VIZ",
                       carpvizl1 = "CBASS-VIZ [L1]")

  cat("CBASS Fit Summary\n")
  cat("====================\n\n")
  cat("Algorithm: ", alg_string, "\n\n")

  cat("Available Visualizations:\n")
  cat(" - Static Dendrogram:   ", x$static, "\n")
  cat(" - Static Heatmap:      ", x$static, "\n")
  cat(" - Interactive Heatmap: ", x$interactive, "\n\n")

  cat("Number of Observations: ", x$n.obs, "\n")
  cat("Number of Variables:    ", x$p.var, "\n\n")

  cat("Pre-processing options:\n")
  cat(" - Global centering: ", x$X.center.global, "\n\n")

  cat("Observation RBF Kernel Weights:\n") # TODO: Add descriptions of what these parameters represent
  cat(" - phi = ", round(x$phi.obs, 3), "\n")
  cat(" - K   = ", x$k.obs, "\n\n")

  cat("Variable RBF Kernel Weights:\n") # TODO: Add descriptions of what these parameters represent
  cat(" - phi = ", round(x$phi.var, 3), "\n")
  cat(" - K   = ", x$k.var, "\n\n")

  cat("Raw Data:\n")
  print(x$X[1:min(5, x$n.obs), 1:min(5, x$p.var)])

  invisible(x)
}
