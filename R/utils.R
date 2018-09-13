

if (getRversion() >= "2.15.1") utils::globalVariables(c("."))

#' Log transformed word count of presidential speeches
#'
#' A dataset of the top 75 most variable log-transformed word counts for
#' each US president aggregated over several speeches
#' (Inaugural, State of the Union, etc.).
#' Stop words have been removed and words have been stemmed.
#'
#' @format A data.frame with 44 rows (one for each president) and 75 columns (log transformed word counts)
#' @details Grover Cleveland was elected president twice (1892 and 1884). For our purposes his speeches are combined.
#' @source \url{http://www.presidency.ucsb.edu}
"presidential_speech"

#' @importFrom Matrix spMatrix
CreateAdjacency <- function(E, sp.pattern, n) {
  adjmat <- Matrix::spMatrix(nrow = n, ncol = n)
  connected.pairs <- matrix(E[which(sp.pattern != 0), ], ncol = 2)
  adjmat@i <- as.integer(as.vector(connected.pairs[, 1]) - 1)
  adjmat@j <- as.integer(as.vector(connected.pairs[, 2]) - 1)
  adjmat@x <- as.double(rep(1, times = nrow(connected.pairs)))
  return(adjmat)
}

CreateClusterGraph <- function(AdjMatrix) {
  igraph::graph_from_adjacency_matrix(AdjMatrix, mode = "undirected")
}

#' @importFrom igraph components
GetClusters <- function(ClusterGraph) {
  igraph::components(ClusterGraph)
}

#' @importFrom Matrix which
#' @importFrom Matrix Matrix
#' @importFrom Matrix t
#' @importFrom Matrix Diagonal
ConvexClusteringPreCompute <- function(X, weights, rho, verbose = FALSE) {
  n <- ncol(X)
  p <- nrow(X)
  # Calcuate edge set
  weight.adj <- WeightAdjacency(weights, n)
  cardE <- sum(weight.adj)
  E <- Matrix::which(weight.adj != 0, arr.ind = TRUE)
  E <- E[order(E[, 1], E[, 2]), ]
  # Precompute indicies
  if (verbose) print("Ind List")
  ind.list <- lapply(1:cardE, function(l) {
    seq(from = (p) * (l - 1) + 1, length.out = p)
  })
  ind.mat <- Reduce(rbind, ind.list) - 1
  if (verbose) print("E1")
  E1.ind.mat <- Reduce(rbind, lapply(1:cardE, function(l) {
    inds <- E[l, ]
    seq(from = (p) * (inds[1] - 1) + 1, length.out = p)
  })) - 1
  if (verbose) print("E2")
  E2.ind.mat <- Reduce(rbind, lapply(1:cardE, function(l) {
    inds <- E[l, ]
    seq(from = (p) * (inds[2] - 1) + 1, length.out = p)
  })) - 1

  if (verbose) print("PreMat")
  # Precompute U-update matrix
  D <- Reduce(`+`, lapply(seq_len(cardE), function(l){
    pos.ind <- E[l, 1]
    neg.ind <- E[l, 2]

    d <- matrix(0, nrow = n, ncol = 1)
    d[pos.ind, 1] <- 1
    d[neg.ind, 1] <- -1
    tcrossprod(d)
  }))
  PreMat <- kronecker(D, Matrix::Diagonal(p)) + (1 / rho) * Matrix::Diagonal(n * p)

  uinit <- Matrix::Matrix(X[TRUE], nrow = p, ncol = n)
  if (verbose) print("Vinit")
  Vmat <- Matrix::Matrix(0, nrow = p, ncol = cardE)
  for (ind in 1:nrow(E)) {
    Vmat[, ind] <- uinit[, E[ind, 1]] - uinit[, E[ind, 2]]
  }
  vinit <- as.vector(Vmat[TRUE])
  ret <- list(
    ind.mat = ind.mat,
    E = E,
    E1.ind.mat = E1.ind.mat,
    E2.ind.mat = E2.ind.mat,
    PreMat = PreMat,
    uinit = uinit,
    vinit = vinit
  )
  return(ret)
}

#' @importFrom dplyr tbl_df
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr slice
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom dplyr left_join
#' @importFrom dplyr full_join
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr tibble
#' @importFrom dplyr bind_rows
#' @importFrom dplyr lead
#' @importFrom dplyr n
#' @importFrom tidyr gather
#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#' @importFrom purrr map
#' @importFrom purrr map2
#' @importFrom stringr str_replace
#' @importFrom stats na.omit
#' @importFrom zoo na.locf
ISP <- function(sp.path, v.path, u.path, lambda.path, cardE) {
  ColLab <- NULL
  SpValue <- NULL
  Iter <- NULL
  SpValueLag <- NULL
  HasChange <- NULL
  ColIndNum <- NULL
  NChanges <- NULL
  data <- NULL
  Rank <- NULL
  ColIndNum.x <- NULL
  ColIndNum.y <- NULL
  ColInd <- NULL
  Lambda <- NULL
  NewLambda <- NULL
  NewU <- NULL
  U <- NULL

  sp.path %>%
    dplyr::tbl_df() %>%
    dplyr::mutate(Iter = 1:n()) %>%
    tidyr::gather(ColLab, SpValue, -Iter) %>%
    dplyr::mutate(
      ColLab = factor(ColLab, levels = paste("V", 1:cardE, sep = ""), ordered = TRUE)
    ) %>%
    dplyr::arrange(Iter, ColLab) %>%
    dplyr::group_by(ColLab) %>%
    dplyr::mutate(
      SpValueLag = dplyr::lag(SpValue)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(Iter != 1) %>%
    # does the sparsity pattern change this iteration?
    dplyr::mutate(
      HasChange = SpValue - SpValueLag
    ) %>%
    # get iterations where sparsity has changed
    dplyr::filter(HasChange > 0) %>%
    dplyr::mutate(
      ColIndNum = as.numeric(stringr::str_replace(as.character(ColLab), "V", ""))
    ) %>%
    dplyr::select(Iter, ColIndNum) %>%
    dplyr::arrange(Iter, ColIndNum) %>%
    dplyr::group_by(Iter) %>%
    # How many changes in this iteration?
    dplyr::mutate(
      NChanges = n()
    ) -> change.frame
  change.frame %>%
    ungroup() %>%
    filter(
      Iter == length(lambda.path)
    ) %>%
    distinct(NChanges) %>%
    unlist() %>%
    unname() -> max.lam.changes
  if (length(max.lam.changes) > 0) {
    if ((max.lam.changes > 1)) {
      lambda.path <- rbind(lambda.path, 1.05 * lambda.path[length(lambda.path)])
      u.path <- cbind(u.path, u.path[, ncol(u.path)])
    }
  }
  change.frame %>%
    dplyr::filter(NChanges > 1) %>%
    dplyr::group_by(Iter) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      tst = purrr::map2(.x = Iter, .y = data, .f = function(x, y) {
        prev.mags <- apply(matrix(v.path[, x - 1], ncol = cardE)[, y$ColIndNum], 2, function(x) {
          sum(x^2)
        })
        data.frame(
          ColIndNum = y$ColIndNum,
          Rank = order(prev.mags)
        )
      })
    ) -> mc.frame

  if (nrow(mc.frame) == 0) {
    dplyr::tibble(
      Iter = 1:length(lambda.path)
    ) %>%
      dplyr::left_join(
        change.frame %>%
          dplyr::filter(NChanges == 1) %>%
          dplyr::select(-NChanges),
        by = c("Iter")
      ) %>%
      dplyr::arrange(Iter) %>%
      select(Iter, ColIndNum) -> IterRankCols
  } else {
    dplyr::tibble(
      Iter = 1:length(lambda.path)
    ) %>%
      dplyr::left_join(
        change.frame %>%
          dplyr::filter(NChanges > 1) %>%
          dplyr::group_by(Iter) %>%
          tidyr::nest() %>%
          dplyr::mutate(
            tst = purrr::map2(.x = Iter, .y = data, .f = function(x, y) {
              prev.mags <- apply(matrix(v.path[, x - 1], ncol = cardE)[, y$ColIndNum], 2, function(x) {
                sum(x^2)
              })
              data.frame(
                ColIndNum = y$ColIndNum,
                Rank = order(prev.mags)
              )
            })
          ) %>%
          dplyr::select(-data) %>%
          tidyr::unnest() %>%
          dplyr::ungroup() %>%
          dplyr::arrange(Iter, Rank) %>%
          dplyr::select(-Rank),
        by = c("Iter")
      ) %>%
      dplyr::left_join(
        change.frame %>%
          dplyr::filter(NChanges == 1) %>%
          dplyr::select(-NChanges),
        by = c("Iter")
      ) %>%
      dplyr::arrange(Iter) %>%
      dplyr::mutate(
        ColIndNum = ifelse(is.na(ColIndNum.x), ColIndNum.y, ColIndNum.x)
      ) %>%
      select(Iter, ColIndNum) -> IterRankCols
  }


  IterRankCols$ColInd <- zoo::na.locf(IterRankCols$ColIndNum, na.rm = FALSE)
  IterRankCols %>%
    dplyr::select(Iter, ColInd) -> IterRankCols





  lapply(1:nrow(IterRankCols), function(idx) {
    indvec <- rep(0, times = cardE)
    indvec[unique(stats::na.omit(IterRankCols$ColInd[1:idx]))] <- 1
    indvec
  }) %>%
    do.call(rbind, .) -> sp.path.inter2



  if (nrow(mc.frame) == 0) {
    dplyr::tibble(
      Iter = 1:length(lambda.path),
      Lambda = lambda.path[Iter]
    ) %>%
      dplyr::mutate(
        Iter = 1:n()
      ) %>%
      dplyr::select(Lambda) %>%
      unlist() %>%
      unname() -> lambda.path.inter2
  } else {
    dplyr::tibble(
      Iter = 1:length(lambda.path),
      Lambda = lambda.path[Iter]
    ) %>%
      dplyr::left_join(
        change.frame %>%
          dplyr::mutate(
            Lambda = lambda.path[Iter]
          ) %>%
          dplyr::filter(NChanges > 1) %>%
          dplyr::group_by(Iter) %>%
          tidyr::nest() %>%
          dplyr::mutate(
            NewLambda = purrr::map2(.x = Iter, .y = data, .f = function(x, y) {
              cur.lam <- unique(y$Lambda)
              next.lam <- lambda.path[x + 1]
              lam.seq <- seq(from = cur.lam, to = next.lam, length.out = nrow(y) + 1)
              lam.seq <- lam.seq[-length(lam.seq)]
              lam.seq
            })
          ) %>%
          dplyr::select(-data) %>%
          tidyr::unnest() %>%
          dplyr::arrange(Iter),
        by = c("Iter")
      ) %>%
      dplyr::mutate(
        Lambda = ifelse(is.na(NewLambda), Lambda, NewLambda),
        Iter = 1:n()
      ) %>%
      dplyr::select(Lambda) %>%
      unlist() %>%
      unname() -> lambda.path.inter2
  }

  seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
  if (nrow(mc.frame) == 0) {
    dplyr::tibble(
      Iter = 1:ncol(u.path)
    ) %>%
      dplyr::mutate(
        Up = purrr::map(.x = Iter, .f = function(x) {
          u.path[, x]
        })
      ) %>%
      dplyr::mutate(
        Iter = 1:n()
      ) %>%
      dplyr::group_by(Iter) %>%
      dplyr::ungroup() -> u.path.inter2
    u.path.inter2$Up %>%
      do.call(cbind, .) -> u.path.inter2
  } else {
    dplyr::tibble(
      Iter = 1:ncol(u.path)
    ) %>%
      dplyr::mutate(
        U = purrr::map(.x = Iter, .f = function(x) {
          u.path[, x]
        })
      ) %>%
      dplyr::left_join(
        change.frame %>%
          dplyr::filter(
            NChanges > 1
          ) %>%
          dplyr::select(-ColIndNum) %>%
          dplyr::ungroup() %>%
          dplyr::group_by(Iter) %>%
          tidyr::nest() %>%
          dplyr::mutate(
            NewU = purrr::map2(.x = Iter, .y = data, .f = function(x, y) {
              cur.u <- u.path[, x]
              next.u <- u.path[, x + 1]
              new.u <- matrix(seq2(from = cur.u, to = next.u, length.out = nrow(y) + 1), nrow = length(cur.u), byrow = TRUE)
              new.u <- new.u[, -ncol(new.u)]
              lapply(seq_len(ncol(new.u)), function(i) {
                new.u[, i]
              })
            })
          ) %>%
          dplyr::select(-data) %>%
          tidyr::unnest(),
        by = c("Iter")
      ) %>%
      dplyr::mutate(
        Iter = 1:n()
      ) %>%
      dplyr::group_by(Iter) %>%
      dplyr::mutate(
        Up = ifelse(is.null(NewU[[1]]), U, NewU)
      ) %>%
      dplyr::ungroup() -> u.path.inter2
    u.path.inter2$Up %>%
      do.call(cbind, .) -> u.path.inter2
  }

  list(sp.path.inter = sp.path.inter2, lambda.path.inter = lambda.path.inter2, u.path.inter = u.path.inter2)
}









#' Compute the minimum number of nearest neighbors required to fully
#' cluster all observations
#'
#' @importFrom Matrix which
#' @param X an n.obs x p.vars matrix
#' @param dense.weights a vector vector of dense weights such as returned by
#' \code{DenseWeights}
#' @return k an integer. The smallest number of nearest neighbors required to
#' fully cluster all observations.
#' @keywords internal
#' @examples
#' \dontrun{
#' library(clustRviz)
#' data("presidential_speech")
#' Xdat <- presidential_speech$X
#' Xdat.preprocessed <- scale(Xdat,center=TRUE,scale=FALSE)
#' dense.weights <- DenseWeights(X = Xdat.preprocessed,phi=1e-3)
#' MinKNN(
#' X = Xdat.preprocessed,
#' dense.weights = dense.weights
#' ) -> k.min
#' }
MinKNN <- function(X, dense.weights) {
  n <- nrow(X)
  weight.adj <- WeightAdjacency(dense.weights, n)
  cardE <- sum(weight.adj)
  E <- Matrix::which(weight.adj != 0, arr.ind = TRUE)
  E <- E[order(E[, 1], E[, 2]), ]
  k <- 0
  n.comp <- 2
  while (n.comp != 1) {
    k <- k + 1
    w.sp <- SparseWeights(X = X, dense.weights = dense.weights, k = k)
    CreateAdjacency(E, sp.pattern = as.numeric(w.sp != 0), n = n) -> adj.full
    CreateClusterGraph(adj.full) -> cluster.full
    n.comp <- GetClusters(cluster.full)$no
  }
  k
}



CreateDendrogram <- function(carp_cluster_path, n_labels, scale = NULL) {
  clust.path <- carp_cluster_path$clust.path[!carp_cluster_path$clust.path.dups]
  lapply(clust.path, function(x) {
    mem <- x$membership
    no <- x$no
    lapply(1:no, function(cl) {
      which(mem == cl)
    })
  }) -> cl.list
  cvx.hclust <- cvxhc(
    clust.path = cl.list,
    gamma.path = carp_cluster_path$lambda.path.inter[!carp_cluster_path$clust.path.dups],
    labels = n_labels
  )

  if (is.null(scale)) {
    min.prop <- min(cvx.hclust$height / max(cvx.hclust$height))
    if (min.prop < .01) {
      scale <- "log"
    } else {
      scale <- "original"
    }
  }
  if (scale == "original") {
    cvx.hclust$height <- cvx.hclust$height
  } else if (scale == "log") {
    cvx.hclust$height <- log(cvx.hclust$height + 1)
  } else {
    stop("CreateDendrogram scale argument must be either 'original' or 'log'")
  }
  max.dend.height <- max(cvx.hclust$height)
  cvx.hclust$height <- cvx.hclust$height * (1 / max.dend.height)
  cvx.hclust
}


#' Compute dense distanced-based gaussian kernel weights for use in CARP or CBASS
#'
#' @param X an n by p matrix with rows the observations and columns the variables.
#' @param phi a number. Scaling parameter used in the gaussian kernel.
#' @param method a string. Passed to the \code{dist} function. See \code{?dist}
#' @param p The power of Minkowski distance. Passed to the \code{dist} function. See \code{?dist}
#' @return dense.weights a numeric vector of weights
#' @importFrom stats dist
#' @keywords internal
#' @examples
#' \dontrun{
#' library(clustRviz)
#' data("presidential_speech")
#' Xdat <- presidential_speech$X
#' Xdat.preprocessed <- scale(Xdat,center=TRUE,scale=FALSE)
#' dense.weights <- DenseWeights(X = Xdat.preprocessed,phi=1e-3)
#' }
DenseWeights <- function(X,
                         phi = 1,
                         method = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
                         p = 2) {
  method <- match.arg(method)
  exp((-1) * phi * (stats::dist(X, method = method, p = p)[TRUE])^2)
}

#' (slowly) Compute sparse kNN weight vector from dense weights
#'
#' @param X an n by p matrix with rows the observations and columns the variables.
#' @param dense.weights a vector dense weights, such as those returned
#' by \code{DenseWeights}
#' @param k a positive integer less than nobs. The number of nearest neighbors.
#' @return sparse.weights a numeric vector of sparse weights
#' @keywords internal
#' @importFrom FNN get.knn
#' @importFrom utils combn
#' @importFrom dplyr %>%
SparseWeights_LEGACY <- function(X, dense.weights, k) {
  nobs <- nrow(X)
  tmp <- FNN::get.knn(X, k = k, algorithm = "brute")
  lapply(1:nobs, function(ind) {
    vec <- rep(0, times = nobs)
    vec[tmp$nn.index[ind, ]] <- 1
    vec[ind] <- 1
    vec
  }) %>%
    do.call(rbind, .) -> mat.pat

  apply(utils::combn(nobs, 2), 2, function(pair.ind) {
    vec1 <- mat.pat[pair.ind[1], ]
    vec2 <- mat.pat[pair.ind[2], ]
    andvec <- vec1 * vec2
    sum(as.logical(andvec[pair.ind[1]]) | as.logical(andvec[pair.ind[2]]))
  }) -> sp.mask
  dense.weights * sp.mask
}

#' Compute sparse kNN weight vector from dense weights
#'
#' @param X an n by p matrix with rows the observations and columns the variables.
#' @param dense.weights a vector dense weights, such as those returned
#' by \code{DenseWeights}
#' @param k a positive integer less than nobs. The number of nearest neighbors.
#' @return sparse.weights a numeric vector of sparse weights
#' @keywords internal
#' @importFrom FNN get.knn
#' @importFrom utils combn
#' @importFrom dplyr %>%
SparseWeights <- function(X, dense.weights, k) {
  nobs <- nrow(X)
  tmp <- FNN::get.knn(X, k = k, algorithm = "brute")
  lapply(1:nobs, function(ind) {
    vec <- rep(0, times = nobs)
    vec[tmp$nn.index[ind, ]] <- 1
    vec[ind] <- 1
    vec
  }) %>%
    do.call(rbind, .) -> mat.pat


  nnvec1 <- mat.pat
  nnvec1[lower.tri(nnvec1, diag = TRUE)] <- NA
  nnvec1 <- t(nnvec1)
  nnvec1.final <- nnvec1[TRUE][!is.na(nnvec1[TRUE])]
  nnvec1.final

  nnvec2 <- t(mat.pat)
  nnvec2[lower.tri(nnvec2, diag = TRUE)] <- NA
  nnvec2 <- t(nnvec2)
  nnvec2.final <- nnvec2[TRUE][!is.na(nnvec2[TRUE])]
  nnvec2.final

  sp.mask <- as.numeric(as.logical(nnvec1.final) | as.logical(nnvec2.final))
  dense.weights * sp.mask
}

#' Construct adjacency matrix induced by weights
#'
#' @param weights a vector of weights such as returned by \code{SparseWeights}
#' @param nobs the number of observations being clustered
#' @param weighted a logical. If \code{FALSE} created unweighted adjacency matrix
#' determined by support of weight vector. If TRUE create weighted adjacency;
#' default is FALSE.
#' @param upper a logical. If \code{TRUE} only return upper triangular of matrix.
#' If \code{FALSE} return symmetric adjacency
#' @return adj a sparse adjacency matrix
#' @keywords internal
#' @importFrom Matrix Matrix
#' @importFrom Matrix t
WeightAdjacency <- function(weights, nobs, weighted = FALSE, upper = TRUE) {
  adj <- Matrix::Matrix(data = 0, nrow = nobs, ncol = nobs, sparse = TRUE)
  if (!weighted) {
    adj[lower.tri(adj, diag = FALSE)] <- as.numeric(weights != 0)
  } else {
    adj[lower.tri(adj, diag = FALSE)] <- weights
  }
  adj <- Matrix::t(adj)
  if (!upper) {
    adj <- adj + Matrix::t(adj)
  }
  adj
}

#' Distances supported by \code{\link[stats]{dist}}
#' @noRd
SUPPORTED_DISTANCES <- c("euclidean",
                         "maximum",
                         "manhattan",
                         "canberra",
                         "binary",
                         "minkowski")

`%not.in%` <- Negate(`%in%`)
