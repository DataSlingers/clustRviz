#' @rdname plot_carp
#' @export
saveviz <- function(x, ...) {
  UseMethod("saveviz", x)
}

#' @param file.name The name of the output file. The type of resulting image
#'                  is determined from the extension.
#' @param dynamic Should the resulting animation be dynamic (animated) or not?
#'                If \code{TRUE}, a dynamic visualization which varies along the
#'                \code{CARP} solution path at a grid given by \code{percent.seq}
#'                is produced (as a \code{GIF}). If \code{FALSE}, a fixed visualization
#'                at a single solution (determined by either \code{percent} or \code{k}
#'                if supplied) is produced.
#' @param percent.seq A grid of values of \code{percent} along which to generate
#'                    dynamic visualizations (if \code{dynamic == TRUE})
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
#' @importFrom gganimate transition_manual anim_save
#' @importFrom animation ani.options saveGIF
#' @importFrom RColorBrewer brewer.pal
#' @rdname plot_carp
#' @export
saveviz.CARP <- function(x,
                         file.name,
                         type = c("dendrogram", "path"),
                         dynamic = TRUE,
                         axis = c("PC1", "PC2"),
                         dend.branch.width = 2,
                         dend.labels.cex = .6,
                         percent,
                         k,
                         percent.seq = seq(from = 0, to = 1, by = .05),
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

  if (!is_logical_scalar(dynamic)) {
    crv_error(sQuote("dynamic"), " must be either TRUE or FALSE.")
  }

  if (!all(vapply(percent.seq, is_percent_scalar, TRUE))) {
    crv_error("All elements of ", sQuote("percent.seq"), " must be between 0 and 1.")
  }

  if (!dynamic) {
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

  ## From here, we can safely assume we are making a dynamic plot (i.e., a GIF)
  file.name <- ensure_gif(file.name)

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
        ggplot2::ggplot(ggplot2::aes(x = V1, y = V2, group = Obs)) +
        gganimate::transition_manual(PlotIdx) +
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
      gganimate::anim_save(filename = file.name,animation = p)
    },
    dendrogram = {
      animation::saveGIF({
        for (pct in percent.seq) {
          carp_dendro_plot(x,
                           percent = pct,
                           dend.branch.width = dend.branch.width,
                           dend.labels.cex = dend.labels.cex,
                           ...)
        }
      }, movie.name = file.name,
         img.name = "carp_dendrogram",
         ani.width  = convert_units(width, from = units, to = "px"),
         ani.height = convert_units(height, from = units, to = "px"),
         clean = TRUE)
      invisible(file.name)
    }
  )
}

#' @param file.name The name of the output file. The type of resulting image
#'                  is determined from the extension.
#' @param dynamic Should the resulting animation be dynamic (animated) or not?
#'                If \code{TRUE}, a dynamic visualization which varies along the
#'                \code{CARP} solution path at a grid given by \code{percent.seq}
#'                is produced (as a \code{GIF}). If \code{FALSE}, a fixed visualization
#'                at a single solution (determined by \code{percent}, \code{k.obs} or
#'                \code{k.var} if supplied) is produced.
#' @param percent.seq A grid of values of \code{percent} along which to generate
#'                    dynamic visualizations (if \code{dynamic == TRUE})
#' @param width The width of the output, given in \code{unit}s
#' @param height The height of the output, given in \code{unit}s
#' @param units The unit in which \code{width} and \code{height} are specified
#' @importFrom stats quantile
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom grDevices adjustcolor colorRampPalette dev.off dev.cur dev.set
#' @importFrom animation ani.options saveGIF
#' @importFrom RColorBrewer brewer.pal
#' @rdname plot_cbass
#' @export
saveviz.CBASS <- function(x,
                          file.name,
                          type = c("heatmap", "obs.dendrogram", "var.dendrogram"),
                          dynamic = TRUE,
                          dend.branch.width = 2,
                          dend.labels.cex = .6,
                          heatrow.label.cex = 1.5,
                          heatcol.label.cex = 1.5,
                          percent,
                          k.obs,
                          k.var,
                          percent.seq = seq(from = 0, to = 1, by = .05),
                          width = 8,
                          height = 5,
                          units = c("in", "cm", "mm", "px"),
                          ...) {

  type       <- match.arg(type)

  if (!is_logical_scalar(dynamic)) {
    crv_error(sQuote("dynamic"), " must be either TRUE or FALSE.")
  }

  if (!all(vapply(percent.seq, is_percent_scalar, TRUE))) {
    crv_error("All elements of ", sQuote("percent.seq"), " must be between 0 and 1.")
  }

  if (!dynamic) {
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

  ## From here, we can safely assume we are making a dynamic plot (i.e., a GIF)
  file.name <- ensure_gif(file.name)

  switch(
    type,
    heatmap = {
      ## Calculate breaks and colors on the raw data so that they are consitent
      ## across frames. (If we use cbass_heatmap_plot's internal fitting, it will
      ## only look at the gradient for a single frame.)
      nbreaks     <- 50
      quant.probs <- seq(0, 1, length.out = nbreaks)
      breaks      <- unique(quantile(x$X[TRUE], probs = quant.probs))
      nbreaks     <- length(breaks)
      heatmap_col <- colorRampPalette(c("blue", "yellow"))(nbreaks - 1)

      animation::saveGIF({
        for (pct in percent.seq) {
          cbass_heatmap_plot(x,
                             percent = pct,
                             heatrow.label.cex = heatrow.label.cex,
                             heatcol.label.cex = heatcol.label.cex,
                             ...,
                             breaks = breaks,
                             heatmap_col = heatmap_col)
        }
      }, movie.name = file.name,
         img.name = "cbass_heatmap",
         ani.width  = convert_units(width, from = units, to = "px"),
         ani.height = convert_units(height, from = units, to = "px"),
         clean = TRUE)
      invisible(file.name)
    },
    obs.dendrogram = {
      animation::saveGIF({
        for (pct in percent.seq) {
          cbass_dendro_plot(x,
                            percent = pct,
                            dend.branch.width = dend.branch.width,
                            dend.labels.cex = dend.labels.cex,
                            type = "obs",
                            ...)
        }
      }, movie.name = file.name,
         img.name = "cbass_observation_dendrogram",
         ani.width  = convert_units(width, from = units, to = "px"),
         ani.height = convert_units(height, from = units, to = "px"),
         clean = TRUE)
      invisible(file.name)
    },
    var.dendrogram = {
      animation::saveGIF({
        for (pct in percent.seq) {
          cbass_dendro_plot(x,
                            percent = pct,
                            dend.branch.width = dend.branch.width,
                            dend.labels.cex = dend.labels.cex,
                            type = "obs",
                            ...)
        }
      }, movie.name = file.name,
         img.name = "cbass_variable_dendrogram",
         ani.width  = convert_units(width, from = units, to = "px"),
         ani.height = convert_units(height, from = units, to = "px"),
         clean = TRUE)
      invisible(file.name)
    }
  )
}
