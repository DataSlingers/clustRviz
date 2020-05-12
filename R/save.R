#' @rdname plot_carp
#' @export
saveviz <- function(x, ...) {
  UseMethod("saveviz", x)
}

#' @param file.name The name of the output file. The type of resulting image
#'                  is determined from the extension. If \code{dynamic = TRUE},
#'                  the extension is changed to \code{.gif} internally. If \code{type = "js"},
#'                  the extension is changed to \code{.html} internally.
#' @param dynamic Should the resulting animation be dynamic (animated) or not?
#'                If \code{TRUE}, a dynamic visualization which varies along the
#'                \code{CARP} solution path at a grid given by \code{percent.seq}
#'                is produced (as a \code{GIF}). If \code{FALSE}, a fixed visualization
#'                at a single solution (determined by either \code{percent} or \code{k}
#'                if supplied) is produced.
#'                Currently ignored when \code{type = "js"}.
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
#' @importFrom gganimate transition_manual anim_save animate
#' @importFrom animation ani.options saveGIF
#' @importFrom RColorBrewer brewer.pal
#' @rdname plot_carp
#' @export
saveviz.CARP <- function(x,
                         file.name,
                         type = c("dendrogram", "path", "js"),
                         dynamic = TRUE,
                         axis = c("PC1", "PC2"),
                         dend.branch.width = 2,
                         dend.labels.cex = .6,
                         percent,
                         k,
                         percent.seq = seq(from = 0, to = 1, by = .01),
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

  if ( (type == "js") && dynamic ){
    crv_warning("dynamic = TRUE is not supported with type = js.")
    dynamic <- FALSE
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
           },
           js = {
             file.name <- ensure_html(file.name)
             carp_heatmaply(x, percent = percent, k = k, file = file.name, ...)
           })
    return(invisible(file.name))
  }

  ## From here, we can safely assume we are making a dynamic plot (i.e., a GIF)
  file.name <- ensure_gif(file.name)

  switch(
    type,
    path = {
      p <- carp_dynamic_path_plot(x,
                                  axis = axis,
                                  percent.seq = percent.seq)
      animation::ani.options(ani.width  = convert_units(width,  from = units, to = "px"),
                             ani.height = convert_units(height, from = units, to = "px"))
      gganimate::anim_save(filename = file.name, animation = animate(p))
      invisible(file.name)
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
#'                  is determined from the extension.If \code{dynamic = TRUE},
#'                  the extension is changed to \code{.gif} internally. If \code{type = "js"},
#'                  the extension is changed to \code{.html} internally.
#' @param dynamic Should the resulting animation be dynamic (animated) or not?
#'                If \code{TRUE}, a dynamic visualization which varies along the
#'                \code{CARP} solution path at a grid given by \code{percent.seq}
#'                is produced (as a \code{GIF}). If \code{FALSE}, a fixed visualization
#'                at a single solution (determined by \code{percent}, \code{k.row} or
#'                \code{k.col} if supplied) is produced.
#'                Currently ignored when \code{type = "js"}.
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
                          type = c("heatmap", "row.dendrogram", "col.dendrogram", "js"),
                          dynamic = TRUE,
                          dend.branch.width = 2,
                          dend.labels.cex = .6,
                          heatrow.label.cex = 1.5,
                          heatcol.label.cex = 1.5,
                          percent,
                          k.row,
                          k.col,
                          percent.seq = seq(from = 0, to = 1, by = .01),
                          width = 8,
                          height = 5,
                          margins = c(5, 5),
                          units = c("in", "cm", "mm", "px"),
                          ...) {

  type       <- match.arg(type)

  if (!is_logical_scalar(dynamic)) {
    crv_error(sQuote("dynamic"), " must be either TRUE or FALSE.")
  }

  if (!all(vapply(percent.seq, is_percent_scalar, TRUE))) {
    crv_error("All elements of ", sQuote("percent.seq"), " must be between 0 and 1.")
  }

  if ( (type == "js") && dynamic ){
    crv_warning("dynamic = TRUE is not supported with type = js.")
    dynamic <- FALSE
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
                                k.row = k.row,
                                k.col = k.col,
                                heatrow.label.cex = heatrow.label.cex,
                                heatcol.label.cex = heatcol.label.cex,
                                margins = margins)
             dev.off()
           },
           row.dendrogram = {
             crv_new_dev_static(file.name, width = width, height = height, units = units)
             cbass_dendro_plot(x,
                               percent = percent,
                               k.row = k.row,
                               k.col = k.col,
                               dend.branch.width = dend.branch.width,
                               dend.labels.cex = dend.labels.cex,
                               type = "row",
                               ...)
             dev.off()
           },
           col.dendrogram = {
             crv_new_dev_static(file.name, width = width, height = height, units = units)
             cbass_dendro_plot(x,
                               percent = percent,
                               k.row = k.row,
                               k.col = k.col,
                               dend.branch.width = dend.branch.width,
                               dend.labels.cex = dend.labels.cex,
                               type = "col",
                               ...)
             dev.off()
           },
           js = {
             file.name <- ensure_html(file.name)
             cbass_heatmaply(x, percent = percent, k.row = k.row, k.col = k.col,
                             file = file.name, ...)
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
                             heatmap_col = heatmap_col,
                             margins = margins)
        }
      }, movie.name = file.name,
         img.name = "cbass_heatmap",
         ani.width  = convert_units(width, from = units, to = "px"),
         ani.height = convert_units(height, from = units, to = "px"),
         clean = TRUE)
      invisible(file.name)
    },
    row.dendrogram = {
      animation::saveGIF({
        for (pct in percent.seq) {
          cbass_dendro_plot(x,
                            percent = pct,
                            dend.branch.width = dend.branch.width,
                            dend.labels.cex = dend.labels.cex,
                            type = "row",
                            ...)
        }
      }, movie.name = file.name,
         img.name = "cbass_row_dendrogram",
         ani.width  = convert_units(width, from = units, to = "px"),
         ani.height = convert_units(height, from = units, to = "px"),
         clean = TRUE)
      invisible(file.name)
    },
    col.dendrogram = {
      animation::saveGIF({
        for (pct in percent.seq) {
          cbass_dendro_plot(x,
                            percent = pct,
                            dend.branch.width = dend.branch.width,
                            dend.labels.cex = dend.labels.cex,
                            type = "col",
                            ...)
        }
      }, movie.name = file.name,
         img.name = "cbass_column_dendrogram",
         ani.width  = convert_units(width, from = units, to = "px"),
         ani.height = convert_units(height, from = units, to = "px"),
         clean = TRUE)
      invisible(file.name)
    }
  )
}
