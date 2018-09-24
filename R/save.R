#' @rdname plot_carp
#' @export
saveviz <- function(x, ...) {
  UseMethod("saveviz", x)
}

#' @param file.name The name of the output file. The type of resulting image
#'                  is determined from the extension.
#' @param image.type The type of image output. If \code{image.type == "static"},
#'                   a fixed visualization of a single solution, with the level
#'                   of regularization being determined by \code{percent} and
#'                   \code{k}. If \code{image.type == "dynamic"}, an animated
#'                   visualization which varies along the solution path.
#' @param percent.seq A grid of values of \code{percent} along which to generate
#'                    dynamic visualizations (if \code{image.type == "dynamic"})
#' @param width The width of the output, given in \code{unit}s
#' @param height The height of the output, given in \code{unit}s
#' @param units The unit in which \code{width} and \code{height} are specified
#' @importFrom stats as.dendrogram
#' @importFrom ggplot2 ggplot aes geom_path geom_point geom_text guides theme element_text
#' @importFrom ggplot2 xlab ylab scale_color_manual ggsave
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr filter select distinct rename mutate left_join select %>%
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom grDevices adjustcolor png dev.off dev.cur dev.set
#' @importFrom gganimate gganimate
#' @importFrom animation ani.options saveGIF
#' @importFrom RColorBrewer brewer.pal
#' @rdname plot_carp
#' @export
saveviz.CARP <- function(x,
                         file.name,
                         type = c("dendrogram", "path"),
                         image.type = c("dynamic", "static"),
                         axis = c("PC1", "PC2"),
                         dend.branch.width = 2,
                         dend.labels.cex = .6,
                         percent,
                         k,
                         percent.seq = seq(from = .05, to = 1, by = .05),
                         width = 8,
                         height = 5,
                         units = c("in", "cm", "mm", "px"),
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

  type       <- match.arg(type)
  image.type <- match.arg(image.type)

  if(image.type == "static"){
    dev_cur <- dev.cur()
    on.exit(dev.set(dev_cur))
    switch(type,
           path = {
             p <- carp_path_plot(x, axis = axis, percent = percent, k = k)
             ggsave(filename = file.name,
                    plot = p,
                    width = width,
                    height = height)
           },
           dendrogram = {
             crv_new_dev_static(file.name, width = width, height = height, units = units)
             carp_dendro_plot(x,
                              percent = percent,
                              k = k,
                              dend.branch.width = dend.branch.width,
                              dend.labels.cex = dend.labels.cex,
                              ...)
             dev.off()
           })
    return(invisible(file.name))
  }

  switch(
    type,
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
      animation::ani.options(ani.width  = convert_units(width,  from = units, to = "px"),
                             ani.height = convert_units(height, from = units, to = "px"))
      gganimate::gganimate(p, file.name)
    },
    dendrogram = {
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
      }, movie.name = file.name,
         img.name = "dend",
         ani.width  = convert_units(width, from = units, to = "px"),
         ani.height = convert_units(height, from = units, to = "px"),
         clean = TRUE)
    }
  )
}

#' @param file.name The name of the output file. The type of resulting image
#'                  is determined from the extension.
#' @param image.type The type of image output. If \code{image.type == "static"},
#'                   a fixed visualization of a single solution, with the level
#'                   of regularization being determined by \code{percent} and
#'                   \code{k}. If \code{image.type == "dynamic"}, an animated
#'                   visualization which varies along the solution path.
#' @param percent.seq A grid of values of \code{percent} along which to generate
#'                    dynamic visualizations (if \code{image.type == "dynamic"})
#' @param width The width of the output, given in \code{unit}s
#' @param height The height of the output, given in \code{unit}s
#' @param units The unit in which \code{width} and \code{height} are specified
#' @importFrom stats as.dendrogram
#' @importFrom ggplot2 ggplot aes geom_path geom_point geom_text guides theme
#' @importFrom ggplot2 element_text xlab ylab scale_color_manual
#' @importFrom ggrepel geom_text_repel
#' @importFrom purrr map_dfr
#' @importFrom dplyr tibble filter select distinct rename mutate left_join select_ %>%
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom grDevices adjustcolor png colorRampPalette dev.off dev.cur dev.set
#' @importFrom gganimate gganimate
#' @importFrom animation ani.options saveGIF
#' @importFrom RColorBrewer brewer.pal
#' @rdname plot_cbass
#' @export
saveviz.CBASS <- function(x,
                          file.name,
                          type = c("heatmap", "obs.dendrogram", "var.dendrogram"),
                          image.type = c("dynamic", "static"),
                          dend.branch.width = 2,
                          dend.labels.cex = .6,
                          heatrow.label.cex = 1.5,
                          heatcol.label.cex = 1.5,
                          percent,
                          k.obs,
                          k.var,
                          percent.seq = seq(from = .05, to = 1, by = .05),
                          width = 8,
                          height = 5,
                          units = c("in", "cm", "mm", "px"),
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

  type       <- match.arg(type)
  image.type <- match.arg(image.type)

  if(image.type == "static"){
    dev_cur <- dev.cur()
    on.exit(dev.set(dev_cur))
    switch(type,
           heatmap = {
             crv_new_dev_static(file.name, width = width, height = height, units = units)
             cbass_heatmap_plot(x,
                                ...,
                                percent = percent,
                                k.obs = k.obs,
                                k.var = k.var,
                                heatrow.label.cex = heatrow.label.cex,
                                heatcol.label.cex = heatcol.label.cex)
             dev.off()
           },
           obs.dendrogram = {
             crv_new_dev_static(file.name, width = width, height = height, units = units)
             cbass_dendro_plot(x,
                               percent = percent,
                               k.obs = k.obs,
                               k.var = k.var,
                               dend.branch.width = dend.branch.width,
                               dend.labels.cex = dend.labels.cex,
                               type = "obs",
                               ...)
             dev.off()
           },
           var.dendrogram = {
             crv_new_dev_static(file.name, width = width, height = height, units = units)
             cbass_dendro_plot(x,
                               percent = percent,
                               k.obs = k.obs,
                               k.var = k.var,
                               dend.branch.width = dend.branch.width,
                               dend.labels.cex = dend.labels.cex,
                               type = "var",
                               ...)
             dev.off()
           })
    return(invisible(file.name))
  }

  switch(
    type,
    heatmap = {
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
      }, movie.name = file.name,
         img.name = "heatmap",
         ani.width  = convert_units(width, from = units, to = "px"),
         ani.height = convert_units(height, from = units, to = "px"),
         clean = TRUE)
    },
    obs.dendrogram = {
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
          x$cbass.dend.obs %>%
            stats::as.dendrogram() %>%
            dendextend::set("branches_lwd", dend.branch.width) %>%
            dendextend::set("labels_cex", dend.labels.cex) %>%
            plot(ylab = "Amount of Regularization")
          my.cols <- grDevices::adjustcolor(c("black", "grey"), alpha.f = .3)
          par(mar = c(14, 7, 2, 1))
          my.rect.hclust(x$cbass.dend.obs, k = ncl, border = 2, my.col.vec = my.cols, lwd = 3)
        }
      }, movie.name = file.name,
         img.name = "obsdend",
         ani.width  = convert_units(width, from = units, to = "px"),
         ani.height = convert_units(height, from = units, to = "px"),
         clean = TRUE)
      ### END Dynamic Obs Dend
    },
    var.dendrogram = {
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
          x$cbass.dend.var %>%
            stats::as.dendrogram() %>%
            dendextend::set("branches_lwd", dend.branch.width) %>%
            dendextend::set("labels_cex", dend.labels.cex) %>%
            plot(ylab = "Amount of Regularization")
          my.cols <- grDevices::adjustcolor(c("black", "grey"), alpha.f = .3)
          par(mar = c(14, 7, 2, 1))
          my.rect.hclust(x$cbass.dend.var, k = ncl, border = 2, my.col.vec = my.cols, lwd = 3)
        }
      }, movie.name = file.name,
         img.name = "vardend",
         ani.width  = convert_units(width, from = units, to = "px"),
         ani.height = convert_units(height, from = units, to = "px"),
         clean = TRUE)
      ### END Dynamic Var Dend
    }
  )
}
