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

#' @importFrom dplyr tbl_df mutate select slice filter arrange left_join full_join
#' @importFrom dplyr group_by ungroup tibble bind_rows lead n
#' @importFrom tidyr gather nest unnest
#' @importFrom purrr map map2
#' @importFrom stringr str_replace
#' @importFrom stats na.omit
#' @importFrom zoo na.locf
#' @importFrom rlang .data
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
        ## We get the magnitude of the previous row differences by reconstructing V
        ## at the `x = Iter` iteration here. Note that, V = DX so the rows are the pairwise
        ## difference of interest, even though we refer to "Col" indices - this is a FIXME
        prev.mags <- rowSums(matrix(v.path[, x - 1], nrow = cardE)[y$ColIndNum,]^2)
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
              ## We get the magnitude of the previous row differences by reconstructing V
              ## at the `x = Iter` iteration here. Note that, V = DX so the rows are the pairwise
              ## difference of interest, even though we refer to "Col" indices - this is a FIXME
              prev.mags <- rowSums(matrix(v.path[, x - 1], nrow = cardE)[y$ColIndNum,]^2)
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
              new.u <- tcrossprod((next.u - cur.u) / NROW(y), seq(0, NROW(y) - 1)) + cur.u
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
    crv_error("CreateDendrogram scale argument must be either 'original' or 'log'")
  }
  max.dend.height <- max(cvx.hclust$height)
  cvx.hclust$height <- cvx.hclust$height * (1 / max.dend.height)
  cvx.hclust
}

#' @noRd
#' @importFrom rlang .data
#' @importFrom stats prcomp
#' @importFrom dplyr tibble %>% mutate group_by ungroup n_distinct
# Post-Process CARP and CBASS results
# This function takes a "correctly" oriented X
ConvexClusteringPostProcess <- function(X,
                                        edge_matrix,
                                        lambda_path,
                                        u_path,
                                        v_path,
                                        v_zero_indices,
                                        labels,
                                        dendrogram_scale,
                                        npcs,
                                        internal_transpose = FALSE){

  n         <- NROW(X)
  p         <- NCOL(X)
  num_edges <- NROW(edge_matrix)

  cluster_path <- ISP(sp.path     = t(v_zero_indices),
                      u.path      = u_path,
                      v.path      = v_path,
                      lambda.path = lambda_path,
                      cardE       = num_edges)

  cluster_path[["clust.path"]] <- get_cluster_assignments(edge_matrix, cluster_path$sp.path.inter, n)
  cluster_path[["clust.path.dups"]] <- duplicated(cluster_path[["clust.path"]], fromList = FALSE)

  U <- array(cluster_path$u.path.inter, dim = c(n, p, length(cluster_path[["clust.path.dups"]])))
  rownames(U) <- rownames(X)
  colnames(U) <- colnames(X)

  if (internal_transpose) {
    ## When looking at the column fusions from CBASS, we want U to be
    ## a 3 tensor of size p by n by K, each slice of which is really U^T
    ##
    ## It's not 100% clear to me why this combination of resizes and transposes works, but
    ## it does, so we're sticking with it...
    dim(U) <- c(p, n, dim(U)[3])
    U <- aperm(U, c(2, 1, 3))
    rownames(U) <- rownames(X)
    colnames(U) <- colnames(X)
  }

  cvx_dendrogram <- CreateDendrogram(cluster_path, labels, dendrogram_scale)

  X_pca <- stats::prcomp(X, scale. = FALSE, center = FALSE)
  rotation_matrix <- X_pca$rotation[, seq_len(npcs)]

  membership_info <- tibble(Iter = rep(seq_along(cluster_path$clust.path), each = n),
                            Obs  = rep(seq_len(n), times = length(cluster_path$clust.path)),
                            Cluster = as.vector(vapply(cluster_path$clust.path, function(x) x$membership, double(n))),
                            Lambda = rep(cluster_path$lambda.path.inter, each = n),
                            ObsLabel = rep(labels, times = length(cluster_path$clust.path))) %>%
                         group_by(.data$Iter) %>%
                         mutate(NCluster = n_distinct(.data$Cluster)) %>%
                         ungroup() %>%
                         mutate(LambdaPercent = .data$Lambda / max(.data$Lambda))

  list(U               = U,
       rotation_matrix = rotation_matrix,
       membership_info = membership_info,
       dendrogram      = cvx_dendrogram)
}

`%not.in%` <- Negate(`%in%`)

# From ?is.integer:
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

is_logical_scalar  <- function(x) {is.logical(x) && (length(x) == 1L) && (!is.na(x))}
is_numeric_scalar  <- function(x) {is.numeric(x) && (length(x) == 1L) && (!is.na(x))}
is_integer_scalar  <- function(x) is_numeric_scalar(x) && is.wholenumber(x)
is_percent_scalar  <- function(x) is_numeric_scalar(x) && (x >= 0) && (x <= 1)
is_positive_scalar <- function(x) is_numeric_scalar(x) && (x > 0)
is_positive_integer_scalar <- function(x) is_integer_scalar(x) && (x > 0)

is_character_scalar <- function(x) {is.character(x) && (length(x) == 1L) && (!is.na(x))}
is_nonempty_character_scalar <- function(x) {is_character_scalar(x) && nzchar(x)}

is_square <- function(x) {is.matrix(x) && (NROW(x) == NCOL(x))}

capitalize_string <- function(x){
  x <- gsub("_", " ", x)
  vapply(strsplit(x, " "),
         function(x) paste(paste0(toupper(substring(x, 1, 1)), substring(x, 2)), collapse = " "),
         character(1))
}

num_unique      <- function(x) length(unique(as.vector(x)))
num_unique_rows <- function(x) NROW(unique(x))
num_unique_cols <- function(x) num_unique_rows(t(x))

unscale_matrix <- function(X,
                           scale=attr(X, "scaled:scale", TRUE),
                           center=attr(X, "scaled:center", TRUE)){
  n <- NROW(X)
  p <- NCOL(X)

  X * matrix(scale, n, p, byrow=TRUE) + matrix(center, n, p, byrow=TRUE)
}

## A very thin wrapper around RColorBrewer::brewer.pal that doesn't warn with
## a few colors
#' @noRd
#' @importFrom RColorBrewer brewer.pal
my_palette <- function(n){
  if(n > 9){
    crv_warning("clustRviz default palette only has nine colors")
  }

  brewer.pal(9, "Set1")[seq_len(n)]
}
