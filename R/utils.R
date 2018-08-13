

if (getRversion() >= "2.15.1") utils::globalVariables(c("."))

#' Log transformed word count of presidential speeches
#'
#' A dataset of the top 75 most variable log-transformed word counts for
#' each US president aggregated over several speeches
#' (Inaugural, State of the Union, etc.).
#' Stop words have been removed and words have been stemmed.
#'
#' @format A list with elements
#' \describe{
#'   \item{labels}{a vector of row labels, one for each presidents}
#'   \item{X}{a 44 by 75 data matrix of log transformed word counts with
#'   the stemmed words the column names}
#'   ...
#' }
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
CreateAdjacencyPath <- function(E, sp.path, n) {
  lapply(1:nrow(sp.path), function(x) {
    CreateAdjacency(E, sp.path[x, ], n)
  })
}
CreateClusterGraph <- function(AdjMatrix) {
  igraph::graph_from_adjacency_matrix(AdjMatrix, mode = "undirected")
}
CreateClusterGraphPath <- function(AdjMatrixList) {
  lapply(AdjMatrixList, function(x) {
    CreateClusterGraph(x)
  })
}

#' @importFrom igraph components
GetClusters <- function(ClusterGraph) {
  igraph::components(ClusterGraph)
}
GetClustersPath <- function(ClusterGraphList) {
  lapply(ClusterGraphList, GetClusters)
}

#' @importFrom Matrix which
#' @importFrom Matrix Matrix
#' @importFrom Matrix t
#' @importFrom Matrix Diagonal
#' @importFrom parallel mclapply
ConvexClusteringPreCompute <- function(X, weights, rho, ncores = 2, verbose = FALSE) {
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
  # precompute matrix in u-update
  PreMat <- Reduce("+", parallel::mclapply(1:cardE, function(l) {
    pos.ind <- E[l, 1]
    neg.ind <- E[l, 2]

    d <- Matrix::Matrix(0, nrow = n, ncol = 1, sparse = TRUE)
    d[pos.ind, 1] <- 1
    d[neg.ind, 1] <- -1
    kronecker(d %*% Matrix::t(d), Matrix::Diagonal(p))
  }, mc.cores = ncores)) + (1 / rho) * Matrix::Diagonal(n * p)
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

iorder <- function(m) {
  N <- nrow(m) + 1
  iorder <- rep(0, N)
  iorder[1] <- m[N - 1, 1]
  iorder[2] <- m[N - 1, 2]
  loc <- 2
  for (i in seq(N - 2, 1))
  {
    for (j in seq(1, loc))
    {
      if (iorder[j] == i) {
        iorder[j] <- m[i, 1]
        if (j == loc) {
          loc <- loc + 1
          iorder[loc] <- m[i, 2]
        } else {
          loc <- loc + 1
          for (k in seq(loc, j + 2)) iorder[k] <- iorder[k - 1]
          iorder[j + 1] <- m[i, 2]
        }
      }
    }
  }
  -iorder
}

#' @importFrom dendextend is.dendrogram
#' @importFrom dendextend heights_per_k.dendrogram
#' @importFrom stats order.dendrogram
#' @importFrom stats cutree
#' @importFrom graphics par
#' @importFrom graphics grconvertX
#' @importFrom graphics grconvertY
#' @importFrom graphics rect
rect.dendrogram <- function(tree, k = NULL, which = NULL, x = NULL, h = NULL, border = 2,
                            cluster = NULL, horiz = FALSE, density = NULL, angle = 45,
                            text = NULL, text_cex = 1, text_col = 1, xpd = TRUE, lower_rect,
                            upper_rect = 0, prop_k_height = 0.5, stop_if_out = FALSE,
                            my.col.vec = NULL,
                            ...) {
  if (!dendextend::is.dendrogram(tree)) {
    stop("x is not a dendrogram object.")
  }
  if (length(h) > 1L | length(k) > 1L) {
    stop("'k' and 'h' must be a scalar(i.e.: of length 1)")
  }
  tree_heights <- dendextend::heights_per_k.dendrogram(tree)[-1]
  tree_order <- stats::order.dendrogram(tree)
  if (!is.null(h)) {
    if (!is.null(k)) {
      stop("specify exactly one of 'k' and 'h'")
    }
    ss_ks <- tree_heights < h
    k <- min(as.numeric(names(ss_ks))[ss_ks])
    k <- max(k, 2)
  }
  else if (is.null(k)) {
    stop("specify exactly one of 'k' and 'h'")
  }
  if (k < 2 | k > length(tree_heights)) {
    if (stop_if_out) {
      stop(gettextf("k must be between 2 and %d", length(tree_heights)),
        domain = NA
      )
    }
    else {
      warning(gettextf("k must be between 2 and %d", length(tree_heights)),
        domain = NA
      )
    }
  }
  if (is.null(cluster)) {
    cluster <- stats::cutree(tree, k = k)
  }
  clustab <- table(cluster)[unique(cluster[tree_order])]
  m <- c(0, cumsum(clustab))
  if (!is.null(x)) {
    if (!is.null(which)) {
      stop("specify exactly one of 'which' and 'x'")
    }
    which <- x
    for (n in seq_along(x)) which[n] <- max(which(m < x[n]))
  }
  else if (is.null(which)) {
    which <- 1L:k
  }
  if (any(which > k)) {
    stop(gettextf(
      "all elements of 'which' must be between 1 and %d",
      k
    ), domain = NA)
  }
  border <- rep_len(border, length(which))
  retval <- list()
  old_xpd <- graphics::par()["xpd"]
  graphics::par(xpd = xpd)
  my.col.vec <- rep(my.col.vec, length.out = length(which))
  for (n in seq_along(which)) {
    next_k_height <- tree_heights[names(tree_heights) ==
      k + 1]
    if (length(next_k_height) == 0) {
      next_k_height <- 0
      prop_k_height <- 1
    }
    if (!horiz) {
      xleft <- m[which[n]] + 0.66
      if (missing(lower_rect)) {
        lower_rect <- graphics::par("usr")[3L] - graphics::strheight("W") *
          (max(nchar(labels(tree))) + 1)
      }
      ybottom <- lower_rect
      xright <- m[which[n] + 1] + 0.33
      ytop <- tree_heights[names(tree_heights) == k] *
        prop_k_height + next_k_height * (1 - prop_k_height) +
        upper_rect
    }
    else {
      ybottom <- m[which[n]] + 0.66
      if (missing(lower_rect)) {
        lower_rect <- graphics::par("usr")[2L] + graphics::strwidth("X") *
          (max(nchar(labels(tree))) + 1)
      }
      xright <- lower_rect
      ytop <- m[which[n] + 1] + 0.33
      xleft <- tree_heights[names(tree_heights) == k] *
        prop_k_height + next_k_height * (1 - prop_k_height) +
        upper_rect
    }
    graphics::rect(xleft, ybottom, xright, ytop,
      border = border[n],

      density = density, angle = angle, my.col.vec[n], ...
    )
    if (!is.null(text)) {
      graphics::text((m[which[n]] + m[which[n] + 1] + 1) / 2, graphics::grconvertY(graphics::grconvertY(
        graphics::par("usr")[3L],
        "user", "ndc"
      ) + 0.02, "ndc", "user"), text[n],
      cex = text_cex, col = text_col
      )
    }
    retval[[n]] <- which(cluster == as.integer(names(clustab)[which[n]]))
  }
  graphics::par(xpd = old_xpd)
  invisible(retval)
}


#' @importFrom graphics rect
#' @importFrom graphics par
#' @importFrom stats cutree
my.rect.hclust <- function(tree, k = NULL, which = NULL, x = NULL, h = NULL, border = 2,
                           cluster = NULL, my.col.vec = NULL, lwd = NULL) {
  if (length(h) > 1L | length(k) > 1L) {
    stop("'k' and 'h' must be a scalar")
  }
  if (!is.null(h)) {
    if (!is.null(k)) {
      stop("specify exactly one of 'k' and 'h'")
    }
    k <- min(which(rev(tree$height) < h))
    k <- max(k, 2)
  }
  else if (is.null(k)) {
    stop("specify exactly one of 'k' and 'h'")
  }
  if (k < 2 | k > length(tree$height)) {
    graphics::rect(
      xleft = graphics::par("usr")[1L],
      xright = graphics::par("usr")[2L],
      ybottom = graphics::par("usr")[3L],
      ytop = graphics::par("usr")[4L],
      border = "red",
      col = my.col.vec[1],
      lwd = lwd
    )
  } else {
    if (is.null(cluster)) {
      cluster <- stats::cutree(tree, k = k)
    }
    clustab <- table(cluster)[unique(cluster[tree$order])]
    m <- c(0, cumsum(clustab))
    if (!is.null(x)) {
      if (!is.null(which)) {
        stop("specify exactly one of 'which' and 'x'")
      }
      which <- x
      for (n in seq_along(x)) which[n] <- max(which(m < x[n]))
    }
    else if (is.null(which)) {
      which <- 1L:k
    }
    if (any(which > k)) {
      stop(gettextf(
        "all elements of 'which' must be between 1 and %d",
        k
      ), domain = NA)
    }
    border <- rep_len(border, length(which))
    #### ADD
    my.col.vec <- rep(my.col.vec, length.out = length(which))
    ####
    retval <- list()
    for (n in seq_along(which)) {
      graphics::rect(m[which[n]] + 0.66, par("usr")[3L], m[which[n] +
        1] + 0.33,
      mean(rev(tree$height)[(k - 1):k]),
      border = border[n],
      col = my.col.vec[n],
      lwd = lwd
      )
      retval[[n]] <- which(cluster == as.integer(names(clustab)[which[n]]))
    }
  }

  # invisible(retval)
}

#' @importFrom stats dist
#' @importFrom stats hclust
#' @importFrom stats median
#' @importFrom stats reorder
#' @importFrom stats density
#' @importFrom stats sd
#' @importFrom stats as.dendrogram
#' @importFrom graphics abline
#' @importFrom graphics axis
#' @importFrom graphics hist
#' @importFrom graphics image
#' @importFrom graphics layout
#' @importFrom graphics lines
#' @importFrom graphics mtext
#' @importFrom graphics plot
#' @importFrom graphics plot.new
#' @importFrom graphics strheight
#' @importFrom graphics strwidth
#' @importFrom graphics text
#' @importFrom graphics title
#' @importFrom gtools invalid
my.heatmap.2 <- function(x, Rowv = TRUE,
                         Colv = if (symm) "Rowv" else TRUE,
                         distfun = stats::dist,
                         hclustfun = stats::hclust,
                         dendrogram = c("both", "row", "column", "none"),
                         symm = FALSE, scale = c("none", "row", "column"),
                         na.rm = TRUE, revC = identical(Colv, "Rowv"),
                         add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) || scale != "none",
                         col = "heat.colors", colsep, rowsep,
                         sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote,
                         notecex = 1, notecol = "cyan", na.color = par("bg"),
                         trace = c("column", "row", "both", "none"),
                         tracecol = "cyan", hline = stats::median(breaks),
                         vline = stats::median(breaks), linecol = tracecol,
                         margins = c(5, 5),
                         ColSideColors, RowSideColors, cexRow = 0.2 + 1 / log10(nr),
                         cexCol = 0.2 + 1 / log10(nc), labRow = NULL, labCol = NULL,
                         key = TRUE, keysize = 1.5,
                         density.info = c("histogram", "density", "none"),
                         denscol = tracecol, symkey = min(x < 0, na.rm = TRUE) || symbreaks,
                         densadj = 0.25, main = NULL, xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL,
                         k.row = 3, k.col = 4, border = "red",
                         Row.hclust = NULL,
                         Col.hclust = NULL, my.col.vec = NULL, lwd = 2, ...) {
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low) / (high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale)) {
    "none"
  } else {
    match.arg(scale)
  }
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col)) {
    col <- get(col, mode = "function")
  }
  if (!missing(breaks) && (scale != "none")) {
    warning(
      "Using scale=\"row\" or scale=\"column\" when breaks are",
      "specified can produce unpredictable results.",
      "Please consider using only one or the other."
    )
  }
  if (is.null(Rowv) || is.na(Rowv)) {
    Rowv <- FALSE
  }
  if (is.null(Colv) || is.na(Colv)) {
    Colv <- FALSE
  } else if (Colv == "Rowv" && !isTRUE(Rowv)) {
    Colv <- FALSE
  }
  if (length(di <- dim(x)) != 2 || !is.numeric(x)) {
    stop("`x' must be a numeric matrix")
  }
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1) {
    stop("`x' must have at least 2 rows and 2 columns")
  }
  if (!is.numeric(margins) || length(margins) != 2) {
    stop("`margins' must be a numeric vector of length 2")
  }
  if (missing(cellnote)) {
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  }
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
      c("both", "row"))) {
      if (is.logical(Colv) && (Colv)) {
        dendrogram <- "column"
      } else {
        dedrogram <- "none"
      }
      warning(
        "Discrepancy: Rowv is FALSE, while dendrogram is `",
        dendrogram, "'. Omitting row dendogram."
      )
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
      c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv)) {
        dendrogram <- "row"
      } else {
        dendrogram <- "none"
      }
      warning(
        "Discrepancy: Colv is FALSE, while dendrogram is `",
        dendrogram, "'. Omitting column dendogram."
      )
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- stats::as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd)) {
      stop("row dendrogram ordering gave index of wrong length")
    }
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- stats::as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd)) {
      stop("row dendrogram ordering gave index of wrong length")
    }
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc) {
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    }
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else {
      colInd <- rowInd
    }
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm) {
      x
    } else {
      t(x)
    } ))
    ddc <- stats::as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd)) {
      stop("column dendrogram ordering gave index of wrong length")
    }
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm) {
      x
    } else {
      t(x)
    } ))
    ddc <- stats::as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd)) {
      stop("column dendrogram ordering gave index of wrong length")
    }
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow)) {
    labRow <- if (is.null(rownames(x))) {
      (1:nr)[rowInd]
    } else {
      rownames(x)
    }
  } else {
    labRow <- labRow[rowInd]
  }
  if (is.null(labCol)) {
    labCol <- if (is.null(colnames(x))) {
      (1:nc)[colInd]
    } else {
      colnames(x)
    }
  } else {
    labCol <- labCol[colInd]
  }
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, stats::sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, stats::sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) <
    1) {
    if (missing(col) || is.function(col)) {
      breaks <- 16
    } else {
      breaks <- length(col) + 1
    }
  }
  if (length(breaks) == 1) {
    if (!symbreaks) {
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
        length = breaks
      )
    } else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function") {
    col <- col(ncol)
  }
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei)) {
    lhei <- c(keysize, 4)
  }
  if (missing(lwid) || is.null(lwid)) {
    lwid <- c(keysize, 4)
  }
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    if (!missing(ColSideColors)) {
      if (!is.character(ColSideColors) || length(ColSideColors) !=
        nc) {
        stop("'ColSideColors' must be a character vector of length ncol(x)")
      }
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] +
        1)
      lhei <- c(lhei[1], 0.2, lhei[2])
    }
    if (!missing(RowSideColors)) {
      if (!is.character(RowSideColors) || length(RowSideColors) !=
        nr) {
        stop("'RowSideColors' must be a character vector of length nrow(x)")
      }
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) -
        1), 1), lmat[, 2] + 1)
      lwid <- c(lwid[1], 0.2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  if (length(lhei) != nrow(lmat)) {
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  }
  if (length(lwid) != ncol(lmat)) {
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  }
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  graphics::layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  if (!missing(RowSideColors)) {
    par(mar = c(margins[1], 0, 0, 0.5))
    graphics::image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
  }
  if (!missing(ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2]))
    graphics::image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
  }
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr")) {
      ddr <- rev(ddr)
    }
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else {
    iy <- 1:nr
  }
  graphics::image(1:nc, 1:nr, x,
    xlim = 0.5 + c(0, nc), ylim = 0.5 +
      c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col,
    breaks = breaks, ...
  )
  retval$carpet <- x
  if (exists("ddr")) {
    retval$rowDendrogram <- ddr
  }
  if (exists("ddc")) {
    retval$colDendrogram <- ddc
  }
  retval$breaks <- breaks
  retval$col <- col
  if (!gtools::invalid(na.color) & any(is.na(x))) {
    mmat <- ifelse(is.na(x), 1, NA)
    graphics::image(1:nc, 1:nr, mmat,
      axes = FALSE, xlab = "", ylab = "",
      col = na.color, add = TRUE
    )
  }
  graphics::axis(1, 1:nc,
    labels = labCol, las = 2, line = -0.5, tick = 0,
    cex.axis = cexCol
  )
  if (!is.null(xlab)) {
    graphics::mtext(xlab, side = 1, line = margins[1] - 1.25)
  }
  graphics::axis(4, iy,
    labels = labRow, las = 2, line = -0.5, tick = 0,
    cex.axis = cexRow
  )
  if (!is.null(ylab)) {
    graphics::mtext(ylab, side = 4, line = margins[2] - 1.25)
  }
  if (!missing(add.expr)) {
    eval(substitute(add.expr))
  }
  if (!missing(colsep)) {
    for (csep in colsep) rect(
        xleft = csep + 0.5, ybottom = 0,
        xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) +
          1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor
      )
  }
  if (!missing(rowsep)) {
    for (rsep in rowsep) rect(
        xleft = 0, ybottom = (ncol(x) +
          1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) +
          1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1,
        col = sepcolor, border = sepcolor
      )
  }
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        graphics::abline(
          v = i - 0.5 + vline.vals, col = linecol,
          lty = 2
        )
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      graphics::lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        graphics::abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      graphics::lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote)) {
    graphics::text(
      x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
      col = notecol, cex = notecex
    )
  }
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    graphics::plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    tree <- Row.hclust
    cluster <- stats::cutree(tree, k = k.row)
    cluster
    clustab <- table(cluster)[unique(cluster[tree$order])]
    clustab
    m <- c(0, cumsum(clustab))
    m
    which <- 1:k.row
    which
    retval <- list()
    my.col.vec.row <- rep(my.col.vec, length.out = length(which))
    if (k.row == 1) {
      rect(
        xleft = par("usr")[1L],
        xright = par("usr")[2L],
        ybottom = par("usr")[3L],
        ytop = par("usr")[4L],
        border = "red",
        col = my.col.vec[1],
        lwd = lwd
      )
    } else {
      for (n in seq_along(which)) {
        rect(par("usr")[2L], m[which[n]] + 0.66,
          mean(rev(tree$height)[(k.row - 1):k.row]),
          m[which[n] + 1] + 0.33,
          border = "red",
          col = my.col.vec.row[n],
          lwd = lwd
        )
      }
    }
  }
  else {
    graphics::plot.new()
  }
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    graphics::plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    tree <- Col.hclust
    cluster <- stats::cutree(tree, k = k.col)
    cluster
    clustab <- table(cluster)[unique(cluster[tree$order])]
    clustab
    m <- c(0, cumsum(clustab))
    m
    which <- 1:k.col
    retval <- list()
    my.col.vec.col <- rep(my.col.vec, length.out = length(which))
    if (k.col == 1) {
      rect(
        xleft = par("usr")[1L],
        xright = par("usr")[2L],
        ybottom = par("usr")[3L],
        ytop = par("usr")[4L], border = "red", col = my.col.vec[1], lwd = lwd
      )
    } else {
      for (n in seq_along(which)) {
        rect(m[which[n]] + 0.66, par("usr")[3L], m[which[n] +
          1] + 0.33,
        mean(rev(tree$height)[(k.col - 1):k.col]),
        border = "red",
        col = my.col.vec.col[n], lwd = lwd
        )
      }
    }
  }
  else {
    graphics::plot.new()
  }
  if (!is.null(main)) {
    graphics::title(main, cex.main = 1.5 * op[["cex.main"]])
  }
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    z <- seq(min.raw, max.raw, length = length(col))
    graphics::image(
      z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
      xaxt = "n", yaxt = "n"
    )
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    graphics::axis(1, at = xv, labels = lv)
    if (scale == "row") {
      graphics::mtext(side = 1, "Row Z-Score", line = 2)
    } else if (scale == "column") {
      graphics::mtext(side = 1, "Column Z-Score", line = 2)
    } else {
      graphics::mtext(side = 1, "Value", line = 2)
    }
    if (density.info == "density") {
      dens <- stats::density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      graphics::lines(dens$x, dens$y / max(dens$y) * 0.95,
        col = denscol,
        lwd = 1
      )
      graphics::axis(2,
        at = pretty(dens$y) / max(dens$y) * 0.95,
        pretty(dens$y)
      )
      graphics::title("Color Key\nand Density Plot")
      par(cex = 0.5)
      graphics::mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- graphics::hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      graphics::lines(hx, hy / max(hy) * 0.95,
        lwd = 1, type = "s",
        col = denscol
      )
      graphics::axis(2, at = pretty(hy) / max(hy) * 0.95, pretty(hy))
      graphics::title("Color Key\nand Histogram")
      par(cex = 0.5)
      graphics::mtext(side = 2, "Count", line = 2)
    }
    else {
      graphics::title("Color Key")
    }
  }
  else {
    graphics::plot.new()
  }
  retval$colorTable <- data.frame(
    low = retval$breaks[-length(retval$breaks)],
    high = retval$breaks[-1], color = retval$col
  )
  invisible(retval)
}

cvxhc <- function(clust.path, gamma.path, labels) {
  h <- sort(gamma.path, decreasing = FALSE)
  # h = log(h + 1)
  h <- h[-1]
  N <- length(clust.path)
  n <- -(1:N) # Tracks group membership
  m <- matrix(0, nrow = N - 1, ncol = 2) # hclust merge output
  for (j in seq(1, N - 1))
  {

    # affected observation indicies
    old.clustering <- clust.path[[j]]
    new.clustering <- clust.path[[j + 1]]
    i <- new.clustering[!(new.clustering %in% old.clustering)][[1]]

    i
    p <- n[i]
    p
    n.neg <- sum(p < 0)
    n.pos <- sum(p > 0)
    n.neg
    n.pos
    # if n.neg > 0 and n.pos ==0, joining two singletons
    if ((n.neg > 0) & (n.pos == 0)) {
      i
      p <- n[i]
      p
    } else if ((n.neg > 0) & (n.pos > 0)) {
      i <- c(i[p < 0][1], i[p > 0][1])
      i
      p <- n[i]
    } else if ((n.neg == 0) & (n.pos > 0)) {
      tmp.clusts <- unique(p)
      i <- c(i[p == tmp.clusts[1]][1], i[p == tmp.clusts[2]][1])
      p <- n[i]
    }
    p

    # R's convention is to order each m[j,] pair as follows:
    p <- p[order(p)]

    m[j, ] <- p
    # Agglomerate this pair and all previous groups they belong to
    # into the current jth group:
    grp <- c(i, which(n %in% n[i[n[i] > 0]]))
    n[grp] <- j
  }
  structure(list(
    merge = m, height = h, order = iorder(m),
    labels = labels, method = "CVX",
    call = match.call(), dist.method = "euclidean"
  ),
  class = "hclust"
  )
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
