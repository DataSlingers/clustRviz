#' Compute CARP solution path
#'
#' \code{CARP} returns a fast approximation to the Convex Clustering
#' solution path along with visualizations such as dendrograms and
#' cluster paths. Visualizations may be static, interactive, or both.
#'
#' \code{CARP} solves the Convex Clustering problem via
#' Algorithmic Regularization Paths. A seqeunce of clustering
#' solutions is returned along with several visualizations.
#'
#'
#' @param X A n by p matrix with rows the observations and columns the variables
#' @param verbose either 0,1,or 2. Higher numbers indicate more verbosity while
#' computing.
#' @param interactive A logical. Should interactive paths and dendrograms be
#' returned?
#' @param static A logical. Should static paths and dendrograms be returned?
#' @param control A list. list of parameters for CARP fitting. see carp.control.
#' @param ... arguements to carp.control
#' @return X the original data matrix
#' @return carp.dend if static==TRUE, a dendrogam representation of the clustering solution path
#' @return carp.cluster.path.vis the CARP solution path
#' @return n.obs the number of observations
#' @return p.var the number of variables
#' @return phi A positive numner used for scaling in RBF kernel
#' @return k an integer >= 1. The number of neighbors used to create sparse weights
#' @return burn.in an integer. The number of initial iterations at a fixed
#' value of (small) lambda_k
#' @return alg.type Which CARP algorithm to perform. Choices are 'carpviz'
#' and 'carp;
#' @return X.center A logical. Should X be centered?
#' @return X.scale A logical. Should X be scaled?
#' @importFrom utils data
#' @importFrom dplyr %>%
#' @importFrom dplyr n
#' @importFrom dplyr tbl_df
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom stats var
#' @export
#' @examples
#' library(clustRviz)
#' data("presidential_speech")
#' Xdat <- presidential_speech[1:10,1:4]
#' carp.fit <- CARP(X=Xdat)
#' carp.fit
CARP <- function(X,
                 verbose = 1,
                 interactive = TRUE,
                 static = TRUE,
                 control = NULL,
                 ...) {
  n.obs <- NROW(X)
  p.var <- NCOL(X)
  Iter <- Cluster <- Lambda <- NULL
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
  extra.args <- list(...)
  if (length(extra.args)) {
    control.args <- names(formals(carp.control))
    indx <- match(names(extra.args), control.args, nomatch = 0L)
    if (any(indx == 0L)) {
      stop(gettextf("Argument %s not matched", names(extra.args)[indx == 0L]), domain = NA)
    }
  }
  internal.control <- carp.control(...)
  if (!is.null(control)) {
    internal.control[names(control)] <- control
  }
  obs.labels <- internal.control$obs.labels
  var.labels <- internal.control$var.labels
  X.center <- internal.control$X.center
  X.scale <- internal.control$X.scale
  k <- internal.control$k
  phi <- internal.control$phi
  rho <- internal.control$rho
  weights <- internal.control$weights
  weight.dist <- internal.control$weight.dist
  weight.dist.p <- internal.control$weight.dist.p
  ncores <- internal.control$ncores
  max.iter <- internal.control$max.iter
  burn.in <- internal.control$burn.in
  alg.type <- internal.control$alg.type
  t <- internal.control$t
  npcs <- internal.control$npcs
  dendrogram.scale <- internal.control$dendrogram.scale



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
      stop("obs.labels should hve length NROW(X)")
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
      stop("var.labels should be have length NCOL(X)")
    }
  }

  if (is.null(npcs)) {
    npcs <- min(4, p.var)
    npcs <- as.integer(npcs)
  } else {
    npcs <- as.integer(npcs)
    if (!is.integer(npcs) | npcs < 2) {
      stop("npcs should be an integer greater than or equal to 2.")
    }
    if (npcs > p.var) {
      stop("npcs should be less than or equal to NCOL(X)")
    }
  }
  if (!is.null(phi)) {
    if (phi <= 0) {
      stop("phi should be positive.")
    }
  }

  if (length(unique(p.labels) != length(p.labels))) {
    colnames(X) <- make.names(p.labels, unique = TRUE)
  } else {
    colnames(X) <- p.labels
  }
  if (length(unique(n.labels)) != length(n.labels)) {
    rownames(X) <- make.names(n.labels, unique = TRUE)
  } else {
    rownames(X) <- n.labels
  }

  # center and scale
  X.orig <- X
  if (X.center | X.scale) {
    X %>%
      scale(center = X.center, scale = X.scale) %>%
      t() -> X
  } else {
    X <- t(X)
  }

  # get weights
  if (is.null(weights)) {
    if (is.null(phi)) {
      phi.vec <- 10^(-10:10)
      sapply(phi.vec, function(phi) {
        stats::var(DenseWeights(X = t(X), phi = phi, method = weight.dist, p = weight.dist.p))
      }) %>%
        which.max() %>%
        phi.vec[.] -> phi
    }
    weights <- DenseWeights(t(X), phi = phi, method = weight.dist, p = weight.dist.p)
    if (is.null(k)) {
      k <- MinKNN(t(X), weights)
    }
    weights <- SparseWeights(X = t(X), dense.weights = weights, k = k)
  } else {
    if (length(weights) != choose(n.obs, 2)) {
      stop("Incorrect weight length")
    }
  }


  if (verbose.basic) message("Pre-computing weight-based edge sets")
  PreCompList <- suppressMessages(ConvexClusteringPreCompute(
    X = X,
    weights = weights,
    ncores = ncores,
    rho = rho
  ))
  cardE <- NROW(PreCompList$E)

  if (verbose.basic) message("Computing CARP Path")
  switch(
    alg.type,
    carpviz = {
      CARPL2_VIS_FRAC(
        x = X[TRUE],
        n = as.integer(n.obs),
        p = as.integer(p.var),
        lambda_init = 1e-8,
        weights = weights[weights != 0],
        uinit = as.matrix(PreCompList$uinit),
        vinit = as.matrix(PreCompList$vinit),
        premat = PreCompList$PreMat,
        IndMat = PreCompList$ind.mat,
        EOneIndMat = PreCompList$E1.ind.mat,
        ETwoIndMat = PreCompList$E2.ind.mat,
        rho = rho,
        max_iter = as.integer(max.iter),
        burn_in = as.integer(burn.in),
        verbose = verbose.deep,
        try_tol = 1e-5,
        ti = 10,
        t_switch = 1.01,
        keep = 1
      ) -> carp.sol.path
    },
    carp = {
      CARPL2_NF_FRAC(
        x = X[TRUE],
        n = as.integer(n.obs),
        p = as.integer(p.var),
        lambda_init = 1e-8,
        t = t,
        weights = weights[weights != 0],
        uinit = as.matrix(PreCompList$uinit),
        vinit = as.matrix(PreCompList$vinit),
        premat = PreCompList$PreMat,
        IndMat = PreCompList$ind.mat,
        EOneIndMat = PreCompList$E1.ind.mat,
        ETwoIndMat = PreCompList$E2.ind.mat,
        rho = rho,
        max_iter = as.integer(max.iter),
        burn_in = as.integer(burn.in),
        verbose = verbose.deep,
        keep = 1
      ) -> carp.sol.path
    },
    carpl1 = {
      CARPL1_NF_FRAC(
        x = X[TRUE],
        n = as.integer(n.obs),
        p = as.integer(p.var),
        lambda_init = 1e-8,
        t = t,
        weights = weights[weights != 0],
        uinit = as.matrix(PreCompList$uinit),
        vinit = as.matrix(PreCompList$vinit),
        premat = PreCompList$PreMat,
        IndMat = PreCompList$ind.mat,
        EOneIndMat = PreCompList$E1.ind.mat,
        ETwoIndMat = PreCompList$E2.ind.mat,
        rho = rho,
        max_iter = as.integer(max.iter),
        burn_in = as.integer(burn.in),
        verbose = verbose.deep,
        keep = 1
      ) -> carp.sol.path
    },
    carpvizl1 = {
      CARPL1_VIS_FRAC(
        x = X[TRUE],
        n = as.integer(n.obs),
        p = as.integer(p.var),
        lambda_init = 1e-8,
        weights = weights[weights != 0],
        uinit = as.matrix(PreCompList$uinit),
        vinit = as.matrix(PreCompList$vinit),
        premat = PreCompList$PreMat,
        IndMat = PreCompList$ind.mat,
        EOneIndMat = PreCompList$E1.ind.mat,
        ETwoIndMat = PreCompList$E2.ind.mat,
        rho = rho,
        max_iter = as.integer(max.iter),
        burn_in = as.integer(burn.in),
        verbose = verbose.deep,
        try_tol = 1e-5,
        ti = 10,
        t_switch = 1.01,
        keep = 1
      ) -> carp.sol.path
    }
  )

  if (verbose.basic) message("Post-processing")
  ISP(
    sp.path = carp.sol.path$v.zero.inds %>% t(),
    v.path = carp.sol.path$v.path,
    u.path = carp.sol.path$u.path,
    lambda.path = carp.sol.path$lambda.path,
    cardE = cardE
  ) -> carp.cluster.path
  carp.cluster.path$sp.path.inter %>% duplicated(fromLast = FALSE) -> sp.path.dups
  adj.path <- CreateAdjacencyPath(PreCompList$E, sp.path = carp.cluster.path$sp.path.inter, n.obs)
  clust.graph.path <- CreateClusterGraphPath(adj.path)
  clust.path <- GetClustersPath(clust.graph.path)
  clust.path.dups <- duplicated(clust.path, fromLast = FALSE)

  carp.cluster.path[["sp.path.dups"]] <- sp.path.dups
  carp.cluster.path[["adj.path"]] <- adj.path
  carp.cluster.path[["clust.graph.path"]] <- clust.graph.path
  carp.cluster.path[["clust.path"]] <- clust.path
  carp.cluster.path[["clust.path.dups"]] <- clust.path.dups

  if (static | interactive) {
    carp.dend <- CreateDendrogram(carp.cluster.path, n.labels, dendrogram.scale)
  } else {
    carp.dend <- NULL
  }
  if (interactive) {
    X.pca <- stats::prcomp(t(X), scale. = FALSE, center = FALSE)
    X.pca.rot <- X.pca$rotation[, 1:npcs]
    lapply(1:length(carp.cluster.path$clust.path), function(iter) {
      U <- t(matrix(carp.cluster.path$u.path.inter[, iter], ncol = n.obs)) %*% X.pca.rot
      names(U) <- paste("PC", 1:npcs)
      U %>%
        as.data.frame() %>%
        dplyr::tbl_df() %>%
        dplyr::mutate(
          Iter = iter,
          Obs = 1:n(),
          Cluster = carp.cluster.path$clust.path[[iter]]$membership,
          Lambda = carp.cluster.path$lambda.path.inter[iter],
          ObsLabel = n.labels
        )
    }) %>%
      do.call(rbind.data.frame, .) %>%
      dplyr::tbl_df() %>%
      dplyr::group_by(Iter) %>%
      dplyr::mutate(
        NCluster = length(unique(Cluster))
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        LambdaPercent = Lambda / max(Lambda)
      ) -> carp.cluster.path.vis
  } else {
    carp.cluster.path.vis <- NULL
  }
  carp.fit <- list(
    X = X.orig,
    carp.dend = carp.dend,
    carp.cluster.path.vis = carp.cluster.path.vis,
    carp.sol.path = carp.sol.path,
    cardE = cardE,
    n.obs = n.obs,
    p.var = p.var,
    phi = phi,
    k = k,
    burn.in = burn.in,
    alg.type = alg.type,
    t = t,
    X.center = X.center,
    X.scale = X.scale,
    static = static,
    interactive = interactive
  )
  class(carp.fit) <- "CARP"
  return(carp.fit)
}

#' Control for CARP fits
#'
#' Parameters for various CARP fitting options
#'
#'
#' @param obs.labels a vector of length n containing observations (row) labels
#' @param var.labels a vector of length p containing variable (column) labels
#' @param X.center A logical. Should X be centered?
#' @param X.scale A logical. Should X be scaled?
#' @param rho A positive number for augmented lagrangian. Not advisable to change.
#' @param weights A vector of positive number of length choose(n,2).
#' @param k an integer >= 1. The number of neighbors used to create sparse weights
#' @param weight.dist a string indicating the distance metric used to calculate weights.
#' @param weight.dist.p The power of the Minkowski distance, if used.
#' @param phi A positive numner used for scaling in RBF kernel
#' @param ncores an integer >= 1. The number of cores to use.
#' @param max.iter an integer. The maximum number of CARP iterations.
#' @param burn.in an integer. The number of initial iterations at a fixed
#' value of (small) lambda_k
#' @param alg.type Which CARP algorithm to perform. Choices are 'carpviz'
#' and 'carp;
#' @param t a number greater than 1. The size of the multiplicitive
#' regularization parameter update. Typical values are: 1.1, 1.05, 1.01, 1.005
#' Determines observation pairs' fusion amount.
#' @param npcs A integer >= 2. The number of principal components to compute
#' for path visualization.
#' @param dendrogram.scale a character string denoting how the scale of dendrogram
#' regularization proportions should be visualized. Choices are 'original'
#' or 'log'; if not provided, a choice is made based on the data provided.
#' @param ... unused generic arguements.
#' @return a list containing the CARP options.
#' @export
carp.control <- function(
                         obs.labels = NULL,
                         var.labels = NULL,
                         X.center = TRUE,
                         X.scale = FALSE,
                         phi = NULL,
                         rho = 1,
                         weights = NULL,
                         k = NULL,
                         weight.dist = "euclidean",
                         weight.dist.p = 2,
                         ncores = as.integer(1),
                         max.iter = as.integer(1e6),
                         burn.in = as.integer(50),
                         alg.type = "carpviz",
                         t = 1.05,
                         npcs = NULL,
                         dendrogram.scale = NULL,
                         ...) {
  if (!is.logical(X.center)) {
    stop("X.center should be either TRUE or FALSE")
  }
  if (!is.logical(X.scale)) {
    stop("X.scale should be either TRUE or FALSE")
  }
  if (rho < 0) {
    stop("rho should be non-negative")
  }
  if (!(weight.dist %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))) {
    stop("unrecognized weight.dist argument; see method arguement of stats::dist for options.")
  }
  if (weight.dist.p <= 0) {
    stop("weight.dist.p should be > 0; see p argument of stats::dist for details.")
  }
  if (!is.integer(ncores) | ncores <= 0) {
    stop("ncores should be a positive integer.")
  }
  if (!is.null(k)) {
    if (!is.integer(k) | k >= 0) {
      stop("k should be a positive integer.")
    }
  }
  if (!is.integer(max.iter) | max.iter <= 0) {
    stop("max.iter should be a positive integer.")
  }
  if (!is.integer(burn.in) | burn.in <= 0 | burn.in >= max.iter) {
    stop("burn.in should be a positive integer greater than max.iter.")
  }
  if (!(alg.type %in% c("carpviz", "carp", "carpl1", "carpvizl1"))) {
    stop("unrecognized alg.type. see help for details.")
  }
  if (t <= 1) {
    stop("t should be greater than 1.")
  }
  if (!is.null(dendrogram.scale)) {
    if (!(dendrogram.scale %in% c("original", "log"))) {
      stop("dendrogram.scale should be one of 'original' or 'log'")
    }
  }
  list(
    obs.labels = obs.labels,
    var.labels = var.labels,
    X.center = X.center,
    X.scale = X.scale,
    rho = rho,
    phi = phi,
    k = k,
    weights = weights,
    weight.dist = weight.dist,
    weight.dist.p = weight.dist.p,
    ncores = ncores,
    max.iter = max.iter,
    burn.in = burn.in,
    alg.type = alg.type,
    t = t,
    npcs = npcs,
    dendrogram.scale = dendrogram.scale
  )
}

#' Print \code{CARP} Results
#'
#' Prints a brief descripton of a fitted \code{CARP} object.
#'
#' Reports number of observations and variables of dataset, any preprocessing
#' done by the \code{\link{CARP}} function, regularization weight information,
#' the variant of \code{CARP} used, and the visualizations returned.
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
                       carlviz   = "CARP-VIZ",
                       carpvizl1 = "CARP-VIZ [L1]")

  cat("CARP Fit Summary\n")
  cat("====================\n\n")
  cat("Algorithm: ", alg.string, "\n\n")

  cat("Available Visualizations:\n")
  cat(" - Static Dendrogram:         ", x$static, "\n")
  cat(" - Static Cluster Path:       ", x$static, "\n")
  cat(" - Interactive Visualization: ", x$interactive, "\n\n")

  cat("Number of Observations: ", x$n.obs, "\n")
  cat("Number of Variables:    ", x$p.var, "\n\n")

  cat("Pre-processing options:\n")
  cat(" - Columnwise centering: ", x$X.center, "\n")
  cat(" - Columnwise scaling:   ", x$X.scale, "\n\n")

  cat("RBF Kernel Weights:\n") # TODO: Add descriptions of what these parameters represent
  cat(" - phi = ", round(x$phi, 3), "\n")
  cat(" - K   = ", x$k, "\n\n")

  cat("Raw Data:\n")
  print(x$X[1:min(5, x$n.obs), 1:min(5, x$p.var)])

  invisible(x)
}
