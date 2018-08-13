#' Method for saving CARP and CBASS visualizations
#'
#' See \code{saveviz.CARP} and \code{saveviz.CBASS} for details
#'
#' @param x a CARP or CBASS object
#' @param ... additional arguements passed to \code{saveviz.CARP} or
#' \code{saveviz.CBASS}
#' @export
saveviz <- function(x, ...) {
  UseMethod("saveviz", x)
}


#' Save CARP visualizations
#'
#' This function will save dynamic CARP visualizations as gifs and static visualizations as image files
#'
#' @param x a CARP object returned by \code{CARP}
#' @param file.name name and location of image to be saved.
#' @param plot.type a string specifying the type of plot to produce. 'dendrogram'
#' produces CARP's cluster dendrogram and 'path' produces the
#' cluster path
#' @param image.type a string specifying the type of image output. 'dynamic'
#' produced a gif of the plot.type along the CARP path. 'static' produces
#' a single image of plot.type at a point along the CARP path.
#' @param static.image.type string specifying the graphics device with which to save
#' static images.
#' @param axis a character vector of length two with elements as 'PC1','PC2',..
#' etc. Specifics which principal component axis to display for the 'path'
#' plot.type
#' @param percent a number between 0 and 1. Specifies how far along the
#' CARP path the static visualizations should display in terms of
#' percent reguarlization.
#' @param k an interger between 1 and n.obs. Specifies how far along the
#' CARP path the static visualizations should display in term of
#' number of clusters.
#' @param percent.seq a vector of numbers between 0 and 1 specifying the
#' positions along the CARP path used to generate dynamic images.
#' @param dend.branch.width a positive number. Line width on dendrograms.
#' @param dend.labels.cex a positive number. Label size on dendrograms.
#' @param dynamic.width dynamic output width in pixels
#' @param dynamic.height dynamic output heigth in pixels
#' @param static.width static output width in inches
#' @param static.height static output height in inches
#' @param ... Unused additional generic arguements
#' @importFrom stats as.dendrogram
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_path
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 ggsave
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr distinct
#' @importFrom dplyr rename
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom dplyr select_
#' @importFrom dplyr %>%
#' @importFrom tools file_ext
#' @importFrom tools file_path_sans_ext
#' @importFrom grDevices adjustcolor
#' @importFrom grDevices png
#' @importFrom grDevices dev.off
#' @importFrom gganimate gganimate
#' @importFrom animation ani.options
#' @importFrom animation saveGIF
#' @importFrom RColorBrewer brewer.pal
#' @export
saveviz.CARP <- function(
                         x,
                         file.name,
                         plot.type = c("path", "dendrogram"),
                         image.type = c("dynamic", "static"),
                         static.image.type = c("png"),
                         axis = c("PC1", "PC2"),
                         dend.branch.width = 2,
                         dend.labels.cex = .6,
                         percent = NULL,
                         k = NULL,
                         percent.seq = seq(from = .05, to = 1, by = .05),
                         dynamic.width = 1200,
                         dynamic.height = 700,
                         static.width = 8,
                         static.height = 5,
                         ...) {
  Iter <- NULL
  Obs <- NULL
  V1 <- NULL
  V2 <- NULL
  ObsLabel <- NULL
  LambdaPercent <- NULL
  PlotIdx <- NULL
  FirstV1 <- NULL
  FirstV2 <- NULL
  FirstObsLabel <- NULL
  NCluster <- NULL

  plot.type <- match.arg(plot.type)
  image.type <- match.arg(image.type)
  static.image.type <- match.arg(static.image.type)
  if (image.type == "static") {
    n.not.null <- sum(
      c(
        !is.null(k),
        !is.null(percent)
      )
    )
    if (n.not.null != 1) {
      stop("Select exactly one of k or percent for static images")
    }
  }
  switch(
    plot.type,
    path = {
      plot.cols <- c(
        axis,
        "Iter",
        "Obs",
        "Cluster",
        "Lambda",
        "ObsLabel",
        "NCluster",
        "LambdaPercent"
      )
      plot.frame <- x$carp.cluster.path.vis[, plot.cols]
      names(plot.frame)[1:2] <- c("V1", "V2")
      plot.frame %>%
        filter(Iter == 1) %>%
        select(Obs, V1, V2, ObsLabel) %>%
        rename(
          FirstV1 = V1,
          FirstV2 = V2,
          FirstObsLabel = ObsLabel
        ) -> plot.frame.first.iter
      plot.frame %>%
        left_join(
          plot.frame.first.iter,
          by = c("Obs")
        ) -> plot.frame
      plot.frame.list <- list()
      switch(
        image.type,
        dynamic = {
          cur.file.ext <- tools::file_ext(file.name)
          if (cur.file.ext != "gif") {
            file.name <- paste(
              tools::file_path_sans_ext(file.name),
              "gif",
              sep = "."
            )
          }
          for (seq.idx in seq_along(percent.seq)) {
            percent <- percent.seq[seq.idx]
            plot.frame %>%
              dplyr::filter(LambdaPercent <= percent) %>%
              dplyr::filter(Iter > x$burn.in) %>%
              dplyr::mutate(
                PlotIdx = seq.idx
              ) -> plot.frame.list[[seq.idx]]
          }

          dplyr::bind_rows(plot.frame.list) -> plot.frame.ani
          plot.frame.ani %>%
            ggplot2::ggplot(ggplot2::aes(x = V1, y = V2, group = Obs, frame = PlotIdx)) +
            ggplot2::geom_path(
              ggplot2::aes(x = V1, y = V2),
              linejoin = "round",
              color = "red",
              size = 1
            ) +
            ggplot2::geom_point(
              ggplot2::aes(x = FirstV1, y = FirstV2),
              color = "black",
              size = I(4)
            ) +
            ggplot2::geom_text(
              ggplot2::aes(x = FirstV1, y = FirstV2, label = FirstObsLabel),
              size = I(6)
            ) +
            ggplot2::guides(color = FALSE, size = FALSE) +
            ggplot2::theme(axis.title = ggplot2::element_text(size = 25)) +
            ggplot2::theme(axis.text = ggplot2::element_text(size = 20)) +
            ggplot2::xlab(axis[1]) +
            ggplot2::ylab(axis[2]) -> p
          animation::ani.options(ani.width = dynamic.width, ani.height = dynamic.height)
          gganimate::gganimate(p, file.name)
        },
        static = {
          if (!is.null(percent)) {
            plot.frame %>%
              dplyr::filter(LambdaPercent <= percent) %>%
              dplyr::filter(Iter > x$burn.in) %>%
              ggplot2::ggplot(ggplot2::aes(x = V1, y = V2, group = Obs)) +
              ggplot2::geom_path(
                ggplot2::aes(x = V1, y = V2),
                linejoin = "round",
                color = "red",
                size = 1
              ) +
              ggplot2::geom_point(
                ggplot2::aes(x = FirstV1, y = FirstV2),
                color = "black",
                size = I(2)
              ) +
              ggplot2::geom_text(
                ggplot2::aes(x = FirstV1, y = FirstV2, label = FirstObsLabel),
                size = I(3)
              ) +
              ggplot2::guides(color = FALSE, size = FALSE) +
              ggplot2::theme(axis.title = ggplot2::element_text(size = 15)) +
              ggplot2::theme(axis.text = ggplot2::element_text(size = 10)) +
              ggplot2::xlab(axis[1]) +
              ggplot2::ylab(axis[2]) -> p
            ggsave(filename = file.name, plot = p, width = static.width, height = static.height, device = static.image.type)
          } else if (!is.null(k)) {
            plot.frame %>%
              dplyr::filter(NCluster >= k) %>%
              dplyr::filter(Iter > x$burn.in) %>%
              ggplot2::ggplot(ggplot2::aes(x = V1, y = V2, group = Obs)) +
              ggplot2::geom_path(
                ggplot2::aes(x = V1, y = V2),
                linejoin = "round",
                color = "red",
                size = 1
              ) +
              ggplot2::geom_point(
                ggplot2::aes(x = FirstV1, y = FirstV2),
                color = "black",
                size = I(2)
              ) +
              ggplot2::geom_text(
                ggplot2::aes(x = FirstV1, y = FirstV2, label = FirstObsLabel),
                size = I(3)
              ) +
              ggplot2::guides(color = FALSE, size = FALSE) +
              ggplot2::theme(axis.title = ggplot2::element_text(size = 15)) +
              ggplot2::theme(axis.text = ggplot2::element_text(size = 10)) +
              ggplot2::xlab(axis[1]) +
              ggplot2::ylab(axis[2]) -> p
            ggsave(filename = file.name, plot = p, width = static.width, height = static.height, device = static.image.type)
          }
        }
      )
    },
    dendrogram = {
      switch(
        image.type,
        dynamic = {
          cur.file.ext <- tools::file_ext(file.name)
          if (cur.file.ext != "gif") {
            file.name <- paste(
              tools::file_path_sans_ext(file.name),
              "gif",
              sep = "."
            )
          }
          animation::saveGIF({
            for (seq.idx in seq_along(percent.seq)) {
              percent <- percent.seq[seq.idx]
              x$carp.cluster.path.vis %>%
                dplyr::filter(LambdaPercent <= percent) %>%
                dplyr::select(NCluster) %>%
                unlist() %>%
                unname() %>%
                min() -> ncl
              x$carp.dend %>%
                stats::as.dendrogram() %>%
                dendextend::set("branches_lwd", dend.branch.width) %>%
                dendextend::set("labels_cex", dend.labels.cex) %>%
                plot(ylab = "Amount of Regularization", cex.lab = 1.5)
              my.cols <- grDevices::adjustcolor(c("grey", "black"), alpha.f = .2)
              par(mar = c(14, 7, 2, 1))
              my.rect.hclust(x$carp.dend, k = ncl, border = 2, my.col.vec = my.cols, lwd = 3)
            }
          }, movie.name = file.name, img.name = "dend", ani.width = dynamic.width, ani.height = dynamic.height, clean = TRUE)
        },
        static = {
          if (!is.null(percent)) {
            x$carp.cluster.path.vis %>%
              dplyr::filter(LambdaPercent <= percent) %>%
              dplyr::select(NCluster) %>%
              unlist() %>%
              unname() %>%
              min() -> ncl
            png(file.name, width = dynamic.width, height = dynamic.height)
            plot.new()
            par(mar = c(14, 7, 2, 1))
            x$carp.dend %>%
              stats::as.dendrogram() %>%
              dendextend::set("branches_lwd", 2) %>%
              dendextend::set("labels_cex", .6) %>%
              plot(ylab = "Amount of Regularization", cex.lab = 1.5)
            my.cols <- grDevices::adjustcolor(c("grey", "black"), alpha.f = .2)
            my.rect.hclust(x$carp.dend, k = ncl, border = 2, my.col.vec = my.cols, lwd = 3)
            grDevices::dev.off()
          } else if (!is.null(k)) {
            ncl <- k
            png(file.name, width = dynamic.width, height = dynamic.height)
            plot.new()
            par(mar = c(14, 7, 2, 1))
            x$carp.dend %>%
              stats::as.dendrogram() %>%
              dendextend::set("branches_lwd", 2) %>%
              dendextend::set("labels_cex", .6) %>%
              plot(ylab = "Amount of Regularization", cex.lab = 1.5)
            my.cols <- grDevices::adjustcolor(c("grey", "black"), alpha.f = .2)
            my.rect.hclust(x$carp.dend, k = ncl, border = 2, my.col.vec = my.cols, lwd = 3)
            grDevices::dev.off()
          }
        }
      )
    }
  )
}


#' Save CBASS visualizations
#'
#' This function will save dynamic CBASS visualizations as gifs and static visualizations as image files
#'
#' @param x a CBASS object returned by \code{CBASS}
#' @param file.name name and location of image to be saved.
#' @param plot.type a string specifying the type of plot to produce. 'obs.dendrogram'
#' and 'var.dendrogram' produces the CBASS cluster dendrogram for the observations and
#' variables respecitvely. 'heatmap' produces a heatmap of the CBASS solution.
#' @param image.type a string specifying the type of image output. 'dynamic'
#' produced a gif of the plot.type along the CARP path. 'static' produces
#' a single image of plot.type at a point along the CARP path.
#' @param static.image.type string specifying the graphics device with which to save
#' static images.
#' @param percent a number between 0 and 1. Specifies how far along the
#' CARP path the static visualizations should display via the amount of
#' regularization.
#' @param k.obs An interger between 1 and \code{n.obs}.
#' Specifies how far along the
#' CARP path the static visualizations should display via the number of unique
#' observation clusters.
#' @param k.var An interger between 1 and \code{p.var}.
#' Specifies how far along the
#' CARP path the static visualizations should display via the number of unique
#' variable clusters.
#' @param percent.seq a vector of numbers between 0 and 1 specifying the
#' positions along the CARP path used to generate dynamic images.
#' @param dend.branch.width a positive number. Line width on dendrograms.
#' @param dend.labels.cex a positive number. Label size on dendrograms.
#' @param heatrow.label.cex heatmap row label size
#' @param heatcol.label.cex heatmap column label size
#' @param dynamic.width dynamic output width in pixels
#' @param dynamic.height dynamic output heigth in pixels
#' @param static.width static output width in inches
#' @param static.height static output height in inches
#' @param ... Unused additional generic arguements
#' @importFrom stats as.dendrogram
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_path
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggrepel geom_text_repel
#' @importFrom purrr map_dfr
#' @importFrom dplyr tibble
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr distinct
#' @importFrom dplyr rename
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom dplyr select_
#' @importFrom dplyr %>%
#' @importFrom tools file_ext
#' @importFrom tools file_path_sans_ext
#' @importFrom grDevices adjustcolor
#' @importFrom grDevices png
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices dev.off
#' @importFrom gganimate gganimate
#' @importFrom animation ani.options
#' @importFrom animation saveGIF
#' @importFrom RColorBrewer brewer.pal
#' @export
saveviz.CBASS <- function(
                          x,
                          file.name,
                          plot.type = c("heatmap", "obs.dendrogram", "var.dendrogram"),
                          image.type = c("dynamic", "static"),
                          static.image.type = c("png"),
                          dend.branch.width = 2,
                          dend.labels.cex = .6,
                          heatrow.label.cex = 1.5,
                          heatcol.label.cex = 1.5,
                          percent = NULL,
                          k.obs = NULL,
                          k.var = NULL,
                          percent.seq = seq(from = .05, to = 1, by = .05),
                          dynamic.width = 1200,
                          dynamic.height = 700,
                          static.width = 8,
                          static.height = 5,
                          ...) {
  Iter <- NULL
  Obs <- NULL
  V1 <- NULL
  V2 <- NULL
  ObsLabel <- NULL
  LambdaPercent <- NULL
  PlotIdx <- NULL
  FirstV1 <- NULL
  FirstV2 <- NULL
  FirstObsLabel <- NULL
  NCluster <- NULL
  Lambda <- NULL
  NObsCl <- NULL
  NVarCl <- NULL
  Percent <- NULL
  cbass.fit <- NULL

  if (image.type == "static") {
    n.not.null <- sum(
      c(
        !is.null(k.obs),
        !is.null(k.var),
        !is.null(percent)
      )
    )
    if (n.not.null != 1) {
      stop("Select exactly one of k.obs, k.var, or percent")
    }
  }
  switch(
    plot.type,
    heatmap = {
      switch(
        image.type,
        dynamic = {
          ### Dynamic Heatmap
          cur.file.ext <- tools::file_ext(file.name)
          if (cur.file.ext != "gif") {
            file.name <- paste(
              tools::file_path_sans_ext(file.name),
              "gif",
              sep = "."
            )
          }
          if (x$X.center.global) {
            X.heat <- x$X
            X.heat <- X.heat - mean(X.heat)
            X.heat <- t(X.heat)
            X <- x$X
            X <- X - mean(X)
            X <- t(X)
          } else {
            X.heat <- t(x$X)
            X <- t(x$X)
          }
          colnames(X.heat) <- x$obs.labels
          rownames(X.heat) <- x$var.labels
          x$cbass.sol.path$lambda.path %>% as.vector() -> lam.seq
          lam.prop.seq <- lam.seq / max(lam.seq)
          nbreaks <- 50
          quant.probs <- seq(0, 1, length.out = nbreaks)
          breaks <- unique(stats::quantile(X[TRUE], probs = quant.probs))
          nbreaks <- length(breaks)
          heatcols <- grDevices::colorRampPalette(c("blue", "yellow"))(nbreaks - 1)
          my.cols <- grDevices::adjustcolor(c("black", "grey"), alpha.f = .3)


          animation::saveGIF({
            for (seq.idx in seq_along(percent.seq)) {
              percent.loop <- percent.seq[seq.idx]
              plt.iter <- which.min(abs(percent.loop - lam.prop.seq))
              # find lambda at iter
              cur.lam <- x$cbass.sol.path$lambda.path[plt.iter]
              cur.lam
              # find lambda closest in column path
              cur.col.lam.ind <- which.min(abs(x$cbass.cluster.path.obs$lambda.path.inter - cur.lam))
              # find clustering solution in column path
              cur.col.clust.assignment <- x$cbass.cluster.path.obs$clust.path[[cur.col.lam.ind]]$membership
              cur.col.clust.labels <- unique(cur.col.clust.assignment)
              cur.col.nclust <- length(cur.col.clust.labels)
              # find lambda closest in row path
              cur.row.lam.ind <- which.min(abs(x$cbass.cluster.path.var$lambda.path.inter - cur.lam))
              # find clustering solution in row path
              cur.row.clust.assignment <- x$cbass.cluster.path.var$clust.path[[cur.row.lam.ind]]$membership
              cur.row.clust.labels <- unique(cur.row.clust.assignment)
              cur.row.nclust <- length(cur.row.clust.labels)
              for (col.label.ind in seq_along(cur.col.clust.labels)) {
                cur.col.label <- cur.col.clust.labels[col.label.ind]
                col.inds <- which(cur.col.clust.assignment == cur.col.label)
                for (row.label.ind in seq_along(cur.row.clust.labels)) {
                  cur.row.label <- cur.row.clust.labels[row.label.ind]
                  row.inds <- which(cur.row.clust.assignment == cur.row.label)
                  mean.value <- mean(X[row.inds, col.inds])
                  X.heat[row.inds, col.inds] <- mean.value
                }
              }
              my.heatmap.2(
                x = X.heat,
                scale = "none",
                Colv = stats::as.dendrogram(x$cbass.dend.obs),
                Rowv = stats::as.dendrogram(x$cbass.dend.var),
                trace = "none",
                density.info = "none",
                key = FALSE,
                breaks = breaks,
                col = heatcols,
                symkey = F,
                Row.hclust = x$cbass.dend.var %>% stats::as.hclust(),
                Col.hclust = x$cbass.dend.obs %>% stats::as.hclust(),
                k.col = cur.col.nclust,
                k.row = cur.row.nclust,
                my.col.vec = my.cols,
                cexRow = heatrow.label.cex,
                cexCol = heatcol.label.cex,
                margins = c(10, 10)
              )
              par(mar = c(14, 7, 2, 1))
            }
          }, movie.name = file.name, img.name = "heatmap", ani.width = dynamic.width, ani.height = dynamic.height, clean = TRUE)

          ### END Dynamic Heatmap
        },
        static = {
          ### Static Heatmap
          if (x$X.center.global) {
            X.heat <- x$X
            X.heat <- X.heat - mean(X.heat)
            X.heat <- t(X.heat)
            X <- x$X
            X <- X - mean(X)
            X <- t(X)
          } else {
            X.heat <- t(x$X)
            X <- t(x$X)
          }
          colnames(X.heat) <- x$obs.labels
          rownames(X.heat) <- x$var.labels
          nbreaks <- 50
          quant.probs <- seq(0, 1, length.out = nbreaks)
          breaks <- unique(stats::quantile(X[TRUE], probs = quant.probs))
          nbreaks <- length(breaks)
          heatcols <- grDevices::colorRampPalette(c("blue", "yellow"))(nbreaks - 1)

          my.cols <- grDevices::adjustcolor(c("black", "grey"), alpha.f = .3)
          lam.vec <- x$cbass.sol.path$lambda.path %>% as.vector()
          max.lam <- max(lam.vec)
          lam.vec %>%
            purrr::map_dfr(.f = function(cur.lam) {
              # find lambda closest in column path
              cur.col.lam.ind <- which.min(abs(x$cbass.cluster.path.obs$lambda.path.inter - cur.lam))
              # find clustering solution in column path
              cur.col.clust.assignment <- x$cbass.cluster.path.obs$clust.path[[cur.col.lam.ind]]$membership
              cur.col.clust.labels <- unique(cur.col.clust.assignment)
              cur.col.nclust <- length(cur.col.clust.labels)
              # find lambda closest in rowumn path
              cur.row.lam.ind <- which.min(abs(x$cbass.cluster.path.var$lambda.path.inter - cur.lam))
              # find clustering solution in rowumn path
              cur.row.clust.assignment <- x$cbass.cluster.path.var$clust.path[[cur.row.lam.ind]]$membership
              cur.row.clust.labels <- unique(cur.row.clust.assignment)
              cur.row.nclust <- length(cur.row.clust.labels)
              dplyr::tibble(
                Lambda = cur.lam,
                NObsCl = cur.col.nclust,
                NVarCl = cur.row.nclust
              )
            }) %>%
            dplyr::mutate(
              Percent = Lambda / max.lam
            ) -> cut.table

          if (!is.null(k.obs)) {
            cut.table %>%
              dplyr::filter(NObsCl <= k.obs) %>%
              dplyr::slice(1) %>%
              dplyr::select(Lambda) %>%
              unlist() %>%
              unname() -> cur.lam
          } else if (!is.null(k.var)) {
            cut.table %>%
              dplyr::filter(NVarCl <= k.var) %>%
              dplyr::slice(1) %>%
              dplyr::select(Lambda) %>%
              unlist() %>%
              unname() -> cur.lam
          } else if (!is.null(percent)) {
            cut.table %>%
              dplyr::filter(Percent >= percent) %>%
              dplyr::slice(1) %>%
              dplyr::select(Lambda) %>%
              unlist() %>%
              unname() -> cur.lam
          } else {
            stop("Select exactly one of k.obs, k.var, or percent")
          }
          # find lambda closest in column path
          cur.col.lam.ind <- which.min(abs(x$cbass.cluster.path.obs$lambda.path.inter - cur.lam))
          # find clustering solution in column path
          cur.col.clust.assignment <- x$cbass.cluster.path.obs$clust.path[[cur.col.lam.ind]]$membership
          cur.col.clust.labels <- unique(cur.col.clust.assignment)
          cur.col.nclust <- length(cur.col.clust.labels)
          # find lambda closest in row path
          cur.row.lam.ind <- which.min(abs(x$cbass.cluster.path.var$lambda.path.inter - cur.lam))
          # find clustering solution in row path
          cur.row.clust.assignment <- x$cbass.cluster.path.var$clust.path[[cur.row.lam.ind]]$membership
          cur.row.clust.labels <- unique(cur.row.clust.assignment)
          cur.row.nclust <- length(cur.row.clust.labels)
          for (col.label.ind in seq_along(cur.col.clust.labels)) {
            cur.col.label <- cur.col.clust.labels[col.label.ind]
            col.inds <- which(cur.col.clust.assignment == cur.col.label)
            for (row.label.ind in seq_along(cur.row.clust.labels)) {
              cur.row.label <- cur.row.clust.labels[row.label.ind]
              row.inds <- which(cur.row.clust.assignment == cur.row.label)
              mean.value <- mean(X[row.inds, col.inds])
              X.heat[row.inds, col.inds] <- mean.value
            }
          }
          png(file.name, width = dynamic.width, height = dynamic.height)
          plot.new()
          par(mar = c(14, 7, 2, 1))
          my.heatmap.2(
            x = X.heat,
            scale = "none",
            Colv = stats::as.dendrogram(x$cbass.dend.obs),
            Rowv = stats::as.dendrogram(x$cbass.dend.var),
            trace = "none",
            density.info = "none",
            key = FALSE,
            breaks = breaks,
            col = heatcols,
            symkey = F,
            Row.hclust = x$cbass.dend.var %>% stats::as.hclust(),
            Col.hclust = x$cbass.dend.obs %>% stats::as.hclust(),
            k.col = cur.col.nclust,
            k.row = cur.row.nclust,
            my.col.vec = my.cols,
            cexRow = heatrow.label.cex,
            cexCol = heatcol.label.cex,
            margins = c(10, 10)
          )
          grDevices::dev.off()
          ### END Static Heatmap
        }
      )
    },
    obs.dendrogram = {
      switch(
        image.type,
        dynamic = {
          ### Dynamic Obs Dend
          cur.file.ext <- tools::file_ext(file.name)
          if (cur.file.ext != "gif") {
            file.name <- paste(
              tools::file_path_sans_ext(file.name),
              "gif",
              sep = "."
            )
          }
          lam.vec <- x$cbass.sol.path$lambda.path %>% as.vector()
          max.lam <- max(lam.vec)
          lam.vec %>%
            purrr::map_dfr(.f = function(cur.lam) {
              # find lambda closest in column path
              cur.col.lam.ind <- which.min(abs(x$cbass.cluster.path.obs$lambda.path.inter - cur.lam))
              # find clustering solution in column path
              cur.col.clust.assignment <- x$cbass.cluster.path.obs$clust.path[[cur.col.lam.ind]]$membership
              cur.col.clust.labels <- unique(cur.col.clust.assignment)
              cur.col.nclust <- length(cur.col.clust.labels)
              # find lambda closest in rowumn path
              cur.row.lam.ind <- which.min(abs(x$cbass.cluster.path.var$lambda.path.inter - cur.lam))
              # find clustering solution in rowumn path
              cur.row.clust.assignment <- x$cbass.cluster.path.var$clust.path[[cur.row.lam.ind]]$membership
              cur.row.clust.labels <- unique(cur.row.clust.assignment)
              cur.row.nclust <- length(cur.row.clust.labels)
              dplyr::tibble(
                Lambda = cur.lam,
                NObsCl = cur.col.nclust,
                NVarCl = cur.row.nclust
              )
            }) %>%
            dplyr::mutate(
              Percent = Lambda / max.lam
            ) -> cut.table

          animation::saveGIF({
            for (seq.idx in seq_along(percent.seq)) {
              percent <- percent.seq[seq.idx]
              cut.table %>%
                dplyr::filter(Percent <= percent) %>%
                dplyr::select(NObsCl) %>%
                unlist() %>%
                unname() %>%
                min() -> ncl
              cbass.fit$cbass.dend.obs %>%
                stats::as.dendrogram() %>%
                dendextend::set("branches_lwd", dend.branch.width) %>%
                dendextend::set("labels_cex", dend.labels.cex) %>%
                plot(ylab = "Amount of Regularization")
              my.cols <- grDevices::adjustcolor(c("black", "grey"), alpha.f = .3)
              par(mar = c(14, 7, 2, 1))
              my.rect.hclust(x$cbass.dend.obs, k = ncl, border = 2, my.col.vec = my.cols, lwd = 3)
            }
          }, movie.name = file.name, img.name = "obsdend", ani.width = dynamic.width, ani.height = dynamic.height, clean = TRUE)
          ### END Dynamic Obs Dend
        },
        static = {
          ### START Static Obs Dend
          if (!is.null(k.obs)) {
            cbass.fit.clustering <- clustering(cbass.fit, k.obs = k.obs)
          } else if (!is.null(k.var)) {
            cbass.fit.clustering <- clustering(cbass.fit, k.var = k.var)
          } else if (!is.null(percent)) {
            cbass.fit.clustering <- clustering(cbass.fit, percent = percent)
          } else {
            stop("Select exactly one of k.obs, k.var, or percent")
          }
          ncl <- length(unique(cbass.fit.clustering$clustering.assignment.obs))
          png(file.name, width = dynamic.width, height = dynamic.height)
          plot.new()
          par(mar = c(14, 7, 2, 1))
          cbass.fit$cbass.dend.obs %>%
            stats::as.dendrogram() %>%
            dendextend::set("branches_lwd", dend.branch.width) %>%
            dendextend::set("labels_cex", dend.labels.cex) %>%
            plot(ylab = "Amount of Regularization")
          my.cols <- grDevices::adjustcolor(c("black", "grey"), alpha.f = .3)
          my.rect.hclust(x$cbass.dend.obs, k = ncl, border = 2, my.col.vec = my.cols, lwd = 3)
          grDevices::dev.off()
          ### END Static Obs Dend
        }
      )
    },
    var.dendrogram = {
      switch(
        image.type,
        dynamic = {
          ### Dynamic Var Dend
          cur.file.ext <- tools::file_ext(file.name)
          if (cur.file.ext != "gif") {
            file.name <- paste(
              tools::file_path_sans_ext(file.name),
              "gif",
              sep = "."
            )
          }
          lam.vec <- x$cbass.sol.path$lambda.path %>% as.vector()
          max.lam <- max(lam.vec)
          lam.vec %>%
            purrr::map_dfr(.f = function(cur.lam) {
              # find lambda closest in column path
              cur.col.lam.ind <- which.min(abs(x$cbass.cluster.path.obs$lambda.path.inter - cur.lam))
              # find clustering solution in column path
              cur.col.clust.assignment <- x$cbass.cluster.path.obs$clust.path[[cur.col.lam.ind]]$membership
              cur.col.clust.labels <- unique(cur.col.clust.assignment)
              cur.col.nclust <- length(cur.col.clust.labels)
              # find lambda closest in rowumn path
              cur.row.lam.ind <- which.min(abs(x$cbass.cluster.path.var$lambda.path.inter - cur.lam))
              # find clustering solution in rowumn path
              cur.row.clust.assignment <- x$cbass.cluster.path.var$clust.path[[cur.row.lam.ind]]$membership
              cur.row.clust.labels <- unique(cur.row.clust.assignment)
              cur.row.nclust <- length(cur.row.clust.labels)
              dplyr::tibble(
                Lambda = cur.lam,
                NObsCl = cur.col.nclust,
                NVarCl = cur.row.nclust
              )
            }) %>%
            dplyr::mutate(
              Percent = Lambda / max.lam
            ) -> cut.table

          animation::saveGIF({
            for (seq.idx in seq_along(percent.seq)) {
              percent <- percent.seq[seq.idx]
              cut.table %>%
                dplyr::filter(Percent <= percent) %>%
                dplyr::select(NVarCl) %>%
                unlist() %>%
                unname() %>%
                min() -> ncl
              cbass.fit$cbass.dend.var %>%
                stats::as.dendrogram() %>%
                dendextend::set("branches_lwd", dend.branch.width) %>%
                dendextend::set("labels_cex", dend.labels.cex) %>%
                plot(ylab = "Amount of Regularization")
              my.cols <- grDevices::adjustcolor(c("black", "grey"), alpha.f = .3)
              par(mar = c(14, 7, 2, 1))
              my.rect.hclust(x$cbass.dend.var, k = ncl, border = 2, my.col.vec = my.cols, lwd = 3)
            }
          }, movie.name = file.name, img.name = "vardend", ani.width = dynamic.width, ani.height = dynamic.height, clean = TRUE)
          ### END Dynamic Var Dend
        },
        static = {
          ### Static Var Dend
          if (!is.null(k.obs)) {
            cbass.fit.clustering <- clustering(cbass.fit, k.obs = k.obs)
          } else if (!is.null(k.var)) {
            cbass.fit.clustering <- clustering(cbass.fit, k.var = k.var)
          } else if (!is.null(percent)) {
            cbass.fit.clustering <- clustering(cbass.fit, percent = percent)
          } else {
            stop("Select exactly one of k.obs, k.var, or percent")
          }
          ncl <- length(unique(cbass.fit.clustering$clustering.assignment.var))
          png(file.name, width = dynamic.width, height = dynamic.height)
          plot.new()
          par(mar = c(14, 7, 2, 1))
          cbass.fit$cbass.dend.var %>%
            stats::as.dendrogram() %>%
            dendextend::set("branches_lwd", dend.branch.width) %>%
            dendextend::set("labels_cex", dend.labels.cex) %>%
            plot(ylab = "Amount of Regularization")
          my.cols <- grDevices::adjustcolor(c("black", "grey"), alpha.f = .3)
          my.rect.hclust(x$cbass.dend.var, k = ncl, border = 2, my.col.vec = my.cols, lwd = 3)
          grDevices::dev.off()
          ### END Static Var Dend
        }
      )
    }
  )
}
