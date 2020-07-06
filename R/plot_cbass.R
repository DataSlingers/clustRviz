#' Visualize the results of Convex BiClustering (\code{CBASS})
#'
#' \code{plot.CBASS} provides a range of ways to visualize the results of convex
#' clustering, including: \itemize{
#' \item Dendrograms, illustrating the nested cluster hierarchy inferred from
#'       the convex clustering solution path (\code{type = "row.dendrogram"} and
#'       \code{type = "col.dendrogram"} for the row (observation) and column
#'       (feature / variable) clusterings, respectively);
#' \item Path plots, showing the coalescence of the estimated cluster centroids
#'       as the regularization parameter is increased (\code{type = "row.dendrogram"}
#'       and \code{type = "col.dendrogram"} for the row (observation) and column
#'       (feature / variable) clusterings, respectively);
#' \item A cluster heatmap, displaying row and column histograms, as well
#'       as the clustered data matrix in a single graphic (\code{type = "heatmap"});
#' \item An interactive Javascript cluster heatmap (\code{type = "heatmap"}); and
#' \item A \code{\link[shiny]{shiny}} app, which can display the clustering solutions
#'       as a "movie" or allow for interactive exploration (\code{type = "interactive"}).
#' }
#'
#' @param x An object of class \code{CBASS} as returned by \code{\link{CBASS}}
#' @param type A string indicating the type of visualization to show (see details above).
#' @param dynamic
#' @param interactive
#' @param axis A character vector of length two indicating which features or principal
#'             components to use as the axes in the path visualizations.
#'             Currently only features like \code{"PC1"} or \code{"PC2"} (indicating
#'             the first principal component projections) are supported.
#' @param percent A number between 0 and 1, giving the regularization level (as
#'                a fraction of the final regularization level used) at which to
#'                assign clusters in the static (\code{type = "dendrogram"} or \code{type = "path"})
#'                plots.
#' @param k.row An integer indicating the desired number of row clusters to be displayed
#'              in the static plots. (If plotting columns, the regularization level
#'              giving \code{k.row} rows clusters is used.)
#' @param k.col An integer indicating the desired number of column clusters to be displayed
#'              in the static plots. (If plotting rows, the regularization level
#'              giving \code{k.col} column clusters is used.)
#' @param ... Additional arguments, which are handled differently for different
#'            values of \code{type}.\itemize{
#'            \item When \code{type} is \code{"heatmap"}, \code{"row.path"},
#'                  \code{"col.path"}, or \code{"interactive"}, the presence of
#'                  unknown arguments triggers an error;
#'            \item when \code{type == "row.dendrogram"} or \code{type == "col.dendrogram"},
#'                  \code{...} is forwarded to \code{\link[stats]{plot.dendrogram}}; and
#'            \item when \code{type == "heatmap"}, \code{...} is forwarded to
#'                  \code{\link[heatmaply]{heatmaply}}.
#'            } See the documentation of the linked functions for details about
#'            additional supported arguments. \code{saveviz} passes arguments
#'            to the corresponding plot \code{type}.
#' @param dend.branch.width Branch width for dendrogram plots (ignored for
#'        other plot types) - must be positive.
#' @param dend.labels.cex Label size for dendrogram plots (ignored for other plot
#'        types) - must be positive.
#' @param dend.ylab.cex Y-axis size for dendrogram plots (ignored for other plot
#'        types) - must be positive.
#' @param heatrow.label.cex heatmap row label size
#' @param heatcol.label.cex heatmap column label size
#' @param margins A vector of length 2 specifying margin sizes. See the \code{margin}
#'                argument to \code{\link[gplots]{heatmap.2}}.
#' @param slider_y
#' @return The value of the return type depends on the \code{type} argument:\itemize{
#'  \item if \code{type \%in\% c("row.dendrogram", "col.dendrogram", "heatmap")},
#'         \code{x} is returned invisibly;
#'  \item if \code{type \%in\% c("row.path", "col.path")}, an object of class
#'        \code{\link[ggplot2]{ggplot}} which can be plotted directly (by invoking
#'        its print method) or modified further by the user is returned;
#'  \item if \code{type = "interactive"}, a \code{\link[shiny]{shiny}} app which can be activated
#'        by invoking its print method.
#' } \code{saveviz.CBASS} always returns \code{file.name} invisibly.
#' @details The \code{\link{saveviz.CBASS}} function provides a unified interface
#'          for exporting \code{CBASS} visualizations to files. For all plots,
#'          at most one of \code{percent}, \code{k.row}, and \code{k.col} must be supplied.
#' @importFrom shiny shinyApp fluidPage titlePanel tabsetPanel fluidRow animationOptions
#' @importFrom shiny tags column plotOutput renderPlot sliderInput
#' @importFrom stats as.dendrogram as.hclust quantile
#' @importFrom grDevices colorRampPalette adjustcolor
#' @export
#' @rdname plot_cbass
#' @examples
#' cbass_fit <- CBASS(presidential_speech)
#' plot(cbass_fit)
#' plot(cbass_fit, type = "heatmap")
#' \dontrun{
#' plot(cbass_fit, type='interactive')
#' }
plot.CBASS <- function(x,
                       ...,
                       type = c("heatmap",
                                "row.dendrogram",
                                "col.dendrogram",
                                "row.path",
                                "col.path"),
                       dynamic = FALSE,
                       interactive = FALSE,
                       percent,
                       k.row,
                       k.col,
                       percent.seq = seq(0, 1, 0.01),
                       dend.branch.width = 2,
                       dend.labels.cex = .6,
                       dend.ylab.cex = 1.2,
                       heatrow.label.cex = 1,
                       heatcol.label.cex = 1,
                       margins = c(5, 5),
                       axis = c("PC1", "PC2"),
                       slider_y = -0.3) {

  type <- match.arg(type)

  switch(
    type,
    row.dendrogram = {
      if (interactive == FALSE){
        if (dynamic == FALSE){
          cbass_dendro_plot(x,
                            percent = percent,
                            k.row = k.row,
                            k.col = k.col,
                            dend.branch.width = dend.branch.width,
                            dend.labels.cex = dend.labels.cex,
                            type = "row",
                            ...)
        } else {
          crv_error("Not implemented.")
        }
      } else {
        cbass_dendro_plotly(x,
                            dynamic = dynamic,
                            percent = percent,
                            percent.seq = percent.seq,
                            k.row = k.row,
                            k.col = k.col,
                            type = "row",
                            slider_y = slider_y,
                            ...)
      }
    },
    col.dendrogram = {
      if (interactive == FALSE){
        if (dynamic == FALSE){
          cbass_dendro_plot(x,
                            percent = percent,
                            k.row = k.row,
                            k.col = k.col,
                            dend.branch.width = dend.branch.width,
                            dend.labels.cex = dend.labels.cex,
                            dend.ylab.cex = dend.ylab.cex,
                            type = "col",
                            ...)
        } else {
          crv_error("Not implemented.")
        }
      } else {
        cbass_dendro_plotly(x,
                            dynamic = dynamic,
                            percent = percent,
                            percent.seq = percent.seq,
                            k.row = k.row,
                            k.col = k.col,
                            type = "col",
                            slider_y = slider_y,
                            ...)
      }
    },
    row.path = {
      if (interactive == FALSE){
        if (dynamic == FALSE){
          cbass_path_plot(x,
                          axis = axis,
                          percent = percent,
                          k.row = k.row,
                          k.col = k.col,
                          ...,
                          type = "row")
        } else {
          crv_error("Not implemented.")
        }
      } else {
        cbass_path_plotly(x,
                          axis = axis,
                          dynamic = dynamic,
                          percent = percent,
                          k.row = k.row,
                          k.col = k.col,
                          percent.seq = seq(0, 1, 0.01),
                          ...,
                          type = "row")
      }
    },
    col.path = {
      if (interactive == FALSE){
        if (dynamic == FALSE){
          cbass_path_plot(x,
                          axis = axis,
                          percent = percent,
                          k.row = k.row,
                          k.col = k.col,
                          ...,
                          type = "col")
        } else {
          crv_error("Not implemented.")
        }
      } else {
        cbass_path_plotly(x,
                          axis = axis,
                          dynamic = dynamic,
                          percent = percent,
                          k.row = k.row,
                          k.col = k.col,
                          percent.seq = seq(0, 1, 0.01),
                          ...,
                          type = "col")
      }
    },
    heatmap = {
      if (dynamic == FALSE){
        if (interactive == FALSE){
          cbass_heatmap_plot(x,
                             ...,
                             percent = percent,
                             k.row = k.row,
                             k.col = k.col,
                             heatrow.label.cex = heatrow.label.cex,
                             heatcol.label.cex = heatcol.label.cex,
                             margins = margins)
        } else {
          cbass_heatmaply(x,
                          ...,
                          percent = percent,
                          k.row = k.row,
                          k.col = k.col)
        }
      } else {
        crv_error("Not implemented.")
      }
    }
  )
}

#' @noRd
#' @importFrom rlang .data
#' @importFrom dplyr filter select left_join pull rename
#' @importFrom ggplot2 ggplot geom_path aes geom_point guides theme element_text xlab ylab
#' @importFrom ggrepel geom_text_repel
cbass_path_plot <- function(x,
                            ...,
                            axis,
                            percent,
                            k.row,
                            k.col,
                            type = c("row", "col")){

  type <- match.arg(type)

  dots <- list(...)
  if ( length(dots) != 0) {
    if (!is.null(names(dots))) {
      nm <- names(dots)
      nm <- nm[nzchar(nm)]
      crv_error("Unknown argument ", sQuote(nm[1]), " passed to ", sQuote("plot.CBASS."))
    } else {
      crv_error("Unknown argument passed to ", sQuote("plot.CBASS."))
    }
  }

  has_percent <- !missing(percent)
  has_k.row   <- !missing(k.row)
  has_k.col   <- !missing(k.col)

  n_args <- has_percent + has_k.row + has_k.col

  show_clusters <- (n_args == 1)

  if(n_args > 1){
    crv_error("At most one of ", sQuote("percent"), " ", sQuote("k.row"), " and ",
              sQuote("k.col"), " may be supplied.")
  }

  if (n_args == 0L) {
    percent <- 1 ## By default, show the whole path
  }

  plot_frame_full <- get_feature_paths(x, axis, type = type) %>% rename(V1 = axis[1], V2 = axis[2])

  if(has_k.row){

    if ( !is_integer_scalar(k.row) ){
      crv_error(sQuote("k"), " must be an integer scalar (vector of length 1).")
    }

    if( k.row <= 0 ) {
      crv_error(sQuote("k.row"), " must be positive.")
    }

    if( k.row > NROW(x$X) ){
      crv_error(sQuote("k.row"), " cannot be more than the rows in the original data set (", NROW(x$X), ").")
    }

    percent <- get_feature_paths(x, features = character(), type = "row") %>%
      select(.data$GammaPercent, .data$NCluster) %>%
      filter(.data$NCluster <= k.row) %>%
      select(.data$GammaPercent) %>%
      summarize(percent = min(.data$GammaPercent)) %>%
      pull
  }

  if(has_k.col){

    if ( !is_integer_scalar(k.col) ){
      crv_error(sQuote("k"), " must be an integer scalar (vector of length 1).")
    }

    if( k.col <= 0 ) {
      crv_error(sQuote("k.col"), " must be positive.")
    }

    if( k.col > NCOL(x$X) ){
      crv_error(sQuote("k.col"), " cannot be more than the columns in the original data set (", NCOL(x$X), ").")
    }

    percent <- get_feature_paths(x, features = character(), type = "col") %>%
      select(.data$GammaPercent, .data$NCluster) %>%
      filter(.data$NCluster <= k.col) %>%
      select(.data$GammaPercent) %>%
      summarize(percent = min(.data$GammaPercent)) %>%
      pull
  }

  if( !is_percent_scalar(percent) ){
    crv_error(sQuote("percent"), " must be a scalar between 0 and 1 (inclusive).")
  }

  ## If percent == min(GammaPercent), keep some (unshrunken) data to plot
  ## This comes up in static path plots when percent = 0 is given
  plot_frame_full <- plot_frame_full %>% filter(.data$GammaPercent <= max(percent, min(.data$GammaPercent)))

  plot_frame_init  <- plot_frame_full %>% filter(.data$Iter == min(.data$Iter))
  plot_frame_final <- plot_frame_full %>% filter(.data$Iter == max(.data$Iter)) %>%
                                          mutate(final_cluster = factor(.data$Cluster))

  plot_frame_full <- left_join(plot_frame_full,
                               plot_frame_final %>% select(.data$Obs, .data$final_cluster),
                               by = "Obs")

  ## FIXME -- It looks like we don't actually have full fusion in `plot_frame_final`
  ##          (even in points which should be in the same cluster...)

  g <- ggplot(mapping = aes(x = .data$V1, y = .data$V2, group = .data$Obs))

  if (show_clusters) {
    g <- g + geom_path(data = plot_frame_full, aes(color = .data$final_cluster), linejoin="round", size=1) +
             geom_point(data = plot_frame_final, aes(color = .data$final_cluster), size = 4)
  } else {
    g <- g + geom_path(data = plot_frame_full, color = "red", linejoin="round", size=1)
  }

  g + geom_point(data = plot_frame_init, color="black", size = 2) +
    geom_text_repel(data = plot_frame_init, mapping = aes(label = .data$ObsLabel), size = 3) +
    guides(color = FALSE, size = FALSE) +
    theme(axis.title = element_text(size = 15),
          axis.text  = element_text(size = 10)) +
    xlab(axis[1]) +
    ylab(axis[2])
}

#' @noRd
#' @importFrom rlang .data
#' @importFrom dplyr filter select summarize pull
#' @importFrom stats as.dendrogram
#' @importFrom dendextend set
#' @importFrom grDevices adjustcolor
cbass_dendro_plot <- function(x,
                             percent,
                             k.row,
                             k.col,
                             dend.branch.width = 2,
                             dend.labels.cex = .6,
                             dend.ylab.cex = 1.2,
                             type = c("row", "col"),
                             ...){

  type <- match.arg(type)

  if(dend.branch.width <= 0){
    crv_error(sQuote("dend.branch.width"), " must be positive.")
  }

  if(dend.labels.cex <= 0){
    crv_error(sQuote("dend.labels.cex"), " must be positive.")
  }

  has_percent <- !missing(percent)
  has_k.row   <- !missing(k.row)
  has_k.col   <- !missing(k.col)
  n_args      <- has_percent + has_k.row + has_k.col

  show_clusters <- (n_args == 1)

  if(n_args > 1){
    crv_error("At most one of ", sQuote("percent"), " ", sQuote("k.row"), " and ",
         sQuote("k.col"), " may be supplied.")
  }

  as.dendrogram(x, type = type) %>%
           set("branches_lwd", dend.branch.width) %>%
           set("labels_cex", dend.labels.cex) %>%
           plot(ylab = "Fraction of Regularization",
                cex.lab = dend.ylab.cex, yaxt = "n", ...)

  axis(2, at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0%", "25%", "50%", "75%", "100%"), las=2)

  if(show_clusters){
    labels <- get_cluster_labels(x, k.row = k.row, k.col = k.col, percent = percent, type = type)
    n_clusters <- nlevels(labels)

    my.cols <- adjustcolor(c("grey", "black"), alpha.f = .2)
    my.rect.hclust(as.hclust(x, type = type), k = n_clusters, border = 2, my.col.vec = my.cols, lwd = 3)
  }

  invisible(x)
}

#' @importFrom grDevices colorRampPalette adjustcolor
#' @importFrom stats quantile as.dendrogram as.hclust
cbass_heatmap_plot <- function(x,
                               ...,
                               percent,
                               k.row,
                               k.col,
                               heatrow.label.cex,
                               heatcol.label.cex,
                               margins,
                               breaks = NULL,
                               heatmap_col = NULL){

  dots <- list(...)

  if ( length(dots) != 0 ){
    crv_error("Unknown arguments passed to ", sQuote("plot.CBASS."))
  }

  has_percent <- !missing(percent)
  has_k.row   <- !missing(k.row)
  has_k.col   <- !missing(k.col)

  n_args <- has_percent + has_k.row + has_k.col

  if(n_args >= 2){
    crv_error("At most one of ", sQuote("percent,"), " ", sQuote("k.row"), " and ",
              sQuote("k.col"), " may be supplied.")
  }

  if(n_args == 0){
    U     <- get_clustered_data(x, percent = 0, refit = TRUE)
    n.row <- NROW(U)
    n.col <- NCOL(U)
  } else {
    U <- get_clustered_data(x, percent = percent, k.row = k.row, k.col = k.col, refit = TRUE)
    # Note that we can't assign k.row / k.col here since we might confuse it with user-supplied arguments
    n.row <- nlevels(get_cluster_labels(x, percent = percent, k.row = k.row, k.col = k.col, type = "row"))
    n.col <- nlevels(get_cluster_labels(x, percent = percent, k.row = k.row, k.col = k.col, type = "col"))
  }

  if (heatrow.label.cex < 0) {
    crv_error(sQuote("heatrow.label.cex"), " must be positive.")
  }

  if (heatcol.label.cex < 0) {
    crv_error(sQuote("heatcol.label.cex"), " must be positive.")
  }

  if (is.null(breaks)) {
    nbreaks <- 50
    quant.probs <- seq(0, 1, length.out = nbreaks)
    ## Rarely (when there are ties involved) the result of unique(quantile(...))
    ## isn't quite sorted, so we fix things up manually
    breaks <- sort(unique(quantile(as.vector(U), probs = quant.probs)))
  }

  if (is.null(heatmap_col)) {
    nbreaks <- length(breaks)
    heatmap_col <- colorRampPalette(c("blue", "yellow"))(nbreaks - 1)
  }

  my.cols  <- adjustcolor(c("black", "grey"), alpha.f = .3)

  my.heatmap.2(
    x = U,
    scale = "none",
    Rowv = as.dendrogram(x, type = "row"),
    Colv = as.dendrogram(x, type = "col"),
    trace = "none",
    density.info = "none",
    key = FALSE,
    breaks = breaks,
    col = heatmap_col,
    symkey = FALSE,
    Row.hclust = as.hclust(x, type = "row"),
    Col.hclust = as.hclust(x, type = "col"),
    k.col = n.col,
    k.row = n.row,
    my.col.vec = my.cols,
    cexRow = heatrow.label.cex,
    cexCol = heatcol.label.cex,
    margins = margins
  )

  invisible(x)
}

#' @noRd
#' Render CBASS results via the heatmaply package
#' @importFrom heatmaply heatmaply
cbass_heatmaply <- function(x,
                            ...,
                            percent,
                            k.row,
                            k.col){

  has_percent <- !missing(percent)
  has_k.row   <- !missing(k.row)
  has_k.col   <- !missing(k.col)

  n_args <- has_percent + has_k.row + has_k.col

  if(n_args >= 2){
    crv_error("At most one of ", sQuote("percent,"), " ", sQuote("k.row"), " and ",
              sQuote("k.col"), " may be supplied.")
  }

  if(n_args == 0){
    U     <- get_clustered_data(x, percent = 0, refit = TRUE)

    heatmaply(U,
              Rowv = as.hclust(x, type = "row"),
              Colv = as.hclust(x, type = "col"),
              ...)
  } else {
    U <- get_clustered_data(x, percent = percent, k.row = k.row, k.col = k.col, refit = TRUE)
    # Note that we can't assign k.row / k.col here since we might confuse it with user-supplied arguments
    n.row <- nlevels(get_cluster_labels(x, percent = percent, k.row = k.row, k.col = k.col, type = "row"))
    n.col <- nlevels(get_cluster_labels(x, percent = percent, k.row = k.row, k.col = k.col, type = "col"))

    heatmaply(U,
              Rowv = as.hclust(x, type = "row"),
              Colv = as.hclust(x, type = "col"),
              k_row = n.row,
              k_col = n.col,
              ...)
  }
}

#' @noRd
#' @importFrom rlang .data
#' @importFrom dplyr filter select left_join pull rename
#' @importFrom plotly add_markers add_paths add_text plot_ly hide_legend highlight style animation_slider
cbass_path_plotly <- function(x,
                              ...,
                              dynamic = FALSE,
                              axis,
                              percent,
                              percent.seq,
                              k.row,
                              k.col,
                              type = c("row", "col")){

  type <- match.arg(type)

  dots <- list(...)
  if ( length(dots) != 0) {
    if (!is.null(names(dots))) {
      nm <- names(dots)
      nm <- nm[nzchar(nm)]
      crv_error("Unknown argument ", sQuote(nm[1]), " passed to ", sQuote("plot.CBASS."))
    } else {
      crv_error("Unknown argument passed to ", sQuote("plot.CBASS."))
    }
  }

  has_percent <- !missing(percent)
  has_k.row   <- !missing(k.row)
  has_k.col   <- !missing(k.col)

  n_args <- has_percent + has_k.row + has_k.col

  show_clusters <- (n_args == 1)

  if(n_args > 1){
    crv_error("At most one of ", sQuote("percent"), " ", sQuote("k.row"), " and ",
              sQuote("k.col"), " may be supplied.")
  }

  if(n_args >= 1 & dynamic == TRUE){
    crv_error("Can't set ", sQuote("percent"), " ", sQuote("k.row"), " and ", sQuote("k.col"), " for dynamic plots")
  }

  if (n_args == 0L) {
    percent <- 1 ## By default, show the whole path
  }

  plot_frame_full <- get_feature_paths(x, axis, type = type) %>% rename(V1 = axis[1], V2 = axis[2])

  if(has_k.row){

    if ( !is_integer_scalar(k.row) ){
      crv_error(sQuote("k"), " must be an integer scalar (vector of length 1).")
    }

    if( k.row <= 0 ) {
      crv_error(sQuote("k.row"), " must be positive.")
    }

    if( k.row > NROW(x$X) ){
      crv_error(sQuote("k.row"), " cannot be more than the rows in the original data set (", NROW(x$X), ").")
    }

    percent <- get_feature_paths(x, features = character(), type = "row") %>%
      select(.data$GammaPercent, .data$NCluster) %>%
      filter(.data$NCluster <= k.row) %>%
      select(.data$GammaPercent) %>%
      summarize(percent = min(.data$GammaPercent)) %>%
      pull
  }

  if(has_k.col){

    if ( !is_integer_scalar(k.col) ){
      crv_error(sQuote("k"), " must be an integer scalar (vector of length 1).")
    }

    if( k.col <= 0 ) {
      crv_error(sQuote("k.col"), " must be positive.")
    }

    if( k.col > NCOL(x$X) ){
      crv_error(sQuote("k.col"), " cannot be more than the columns in the original data set (", NCOL(x$X), ").")
    }

    percent <- get_feature_paths(x, features = character(), type = "col") %>%
      select(.data$GammaPercent, .data$NCluster) %>%
      filter(.data$NCluster <= k.col) %>%
      select(.data$GammaPercent) %>%
      summarize(percent = min(.data$GammaPercent)) %>%
      pull
  }

  if( !is_percent_scalar(percent) ){
    crv_error(sQuote("percent"), " must be a scalar between 0 and 1 (inclusive).")
  }

  ## If percent == min(GammaPercent), keep some (unshrunken) data to plot
  ## This comes up in static path plots when percent = 0 is given
  plot_frame_full <- plot_frame_full %>% filter(.data$GammaPercent <= max(percent, min(.data$GammaPercent)))

  plot_frame_init  <- plot_frame_full %>% filter(.data$Iter == min(.data$Iter))
  plot_frame_final <- plot_frame_full %>% filter(.data$Iter == max(.data$Iter)) %>%
    mutate(final_cluster = factor(.data$Cluster))

  plot_frame_full <- left_join(plot_frame_full,
                               plot_frame_final %>% select(.data$Obs, .data$final_cluster),
                               by = "Obs")

  plot_frame_interactive <- bind_rows(lapply(percent.seq, function(pct){
    ## Make a list of things to plot at each "percent" and then combine
    plot_frame_full %>% filter(.data$GammaPercent <= pct & .data$GammaPercent > pct-0.01) %>%
      mutate(Regularization = pct)
  }))

  mytext <- paste(plot_frame_init$ObsLabel)

  if (dynamic == FALSE){
    if (show_clusters) {
      path_static <- plot_ly() %>%
        add_markers(
          data = plot_frame_init,
          ids = ~Obs,
          x = ~V1,
          y = ~V2,
          color = I("black"),
          size = 2) %>%
        add_markers(
          data = plot_frame_final,
          ids = ~Obs,
          x = ~V1,
          y = ~V2,
          color = I(~final_cluster),
          size = 3)

      for (i in seq_along(plot_frame_init$ObsLabel)){
        path_static <- path_static %>%
          add_paths(
            data = plot_frame_interactive[plot_frame_interactive$Obs==i,],
            color = I(~final_cluster),
            ids = ~Obs,
            x = ~V1,
            y = ~V2)
      }

      path_static %>%
        add_text(data = plot_frame_init,
                 x = ~V1,
                 y = ~V2,
                 text = ~ObsLabel,
                 #size = label_size,
                 inherit = FALSE) %>%
        hide_legend() %>%
        style(text=mytext, hoverinfo = "text", traces = 1) %>%
        style(hoverinfo = "none", traces = c(2:(length(mytext)+2)))
    } else{
      path_static <- plot_ly() %>%
        add_markers(
          data = plot_frame_init,
          ids = ~Obs,
          x = ~V1,
          y = ~V2,
          color = I("black"),
          size = 2) %>%
        add_markers(
          data = plot_frame_final,
          ids = ~Obs,
          x = ~V1,
          y = ~V2,
          color = I("red"),
          size = 3)

      for (i in seq_along(plot_frame_init$ObsLabel)){
        path_static <- path_static %>%
          add_paths(
            data = plot_frame_interactive[plot_frame_interactive$Obs==i,],
            color = I("red"),
            ids = ~Obs,
            x = ~V1,
            y = ~V2)
      }

      path_static %>%
        add_text(data = plot_frame_init,
                 x = ~V1,
                 y = ~V2,
                 text = ~ObsLabel,
                 #size = label_size,
                 inherit = FALSE) %>%
        hide_legend() %>%
        style(text=mytext, hoverinfo = "text", traces = 1) %>%
        style(hoverinfo = "none", traces = c(2:(length(mytext)+2)))
    }

  }
  else{
    plot_frame_animation <- bind_rows(lapply(percent.seq, function(pct){
      ## Make a list of things to plot at each "percent" and then combine
      plot_frame_full %>% filter(.data$GammaPercent <= pct) %>%
        mutate(Regularization = pct)
    }))

    path_dynamic<-plot_ly()

    for (i in seq_along(plot_frame_init$ObsLabel)){
      path_dynamic <- path_dynamic %>%
        add_paths(
          data = plot_frame_animation[plot_frame_animation$Obs==i,],
          color = I("red"),
          ids = ~Obs,
          #labels = ~ObsLabel,
          x = ~V1,
          y = ~V2,
          frame = ~Regularization*100)
    }

    path_dynamic %>%
      add_markers(
        data = plot_frame_interactive,
        ids = ~Obs,
        x = ~V1,
        y = ~V2,
        color = I("red"),
        size = 2.5,
        frame = ~Regularization*100) %>%
      add_markers(
        data = plot_frame_init,
        ids = ~Obs,
        x = ~V1,
        y = ~V2,
        color = I("black"),
        size = 0.5,
        inherit = FALSE) %>%
      add_text(data = plot_frame_init,
               x = ~V1,
               y = ~V2,
               text = ~ObsLabel,
               #size = label_size,
               inherit = FALSE) %>%
      hide_legend() %>%
      style(text=mytext, hoverinfo = "text", traces = c(length(mytext)+2)) %>%
      style(hoverinfo = "none", traces = c(1:(length(mytext)+1))) %>%
      animation_slider(currentvalue = list(prefix = "Regularization: ",suffix = "%")) #%>%
  }
}

#' @noRd
#' @importFrom rlang .data
#' @importFrom dplyr filter select summarize pull desc
#' @importFrom stats as.dendrogram is.leaf
#' @importFrom dendextend get_nodes_xy as.ggdend
#' @importFrom plotly add_segments add_markers add_text plot_ly hide_legend highlight animation_slider
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tibble as_tibble
cbass_dendro_plotly <- function(x,
                              dynamic = FALSE,
                              percent,
                              percent.seq,
                              k.row,
                              k.col,
                              type = c("row", "col"),
                              slider_y = slider_y,
                              ...){

  type <- match.arg(type)

  has_percent <- !missing(percent)
  has_k.row   <- !missing(k.row)
  has_k.col   <- !missing(k.col)
  n_args      <- has_percent + has_k.row + has_k.col

  show_clusters <- (n_args == 1)

  if(n_args > 1){
    crv_error("At most one of ", sQuote("percent"), " ", sQuote("k.row"), " and ",
              sQuote("k.col"), " may be supplied.")
  }

  d <- as.dendrogram(x, type=type)
  # get x/y locations of every node in the tree
  m <- get_nodes_xy(d)
  colnames(m) <- c("x", "y")
  allXY <- as_tibble(m)
  # get non-zero heights so we can split on them and find the relevant labels
  non0 <- allXY[["y"]][allXY[["y"]] > 0]
  # label is a list-column since non-zero heights have multiple labels
  # for now, we just have access to terminal node labels
  labs <- labels(d)
  allXY$label <- vector("list", nrow(allXY))
  allXY$label[[1]] <- labs
  allXY$label[allXY$y == 0] <- labs

  # collect all the *unique* non-trivial nodes
  nodes <- list()
  for (i in non0) {
    dsub <- cut(d, i)$lower
    for (j in seq_along(dsub)) {
      s <- dsub[[j]]
      if (is.leaf(s)) next
      if (any(vapply(nodes, function(x) identical(x, s), logical(1)))) next
      nodes[[length(nodes) + 1]] <- s
    }
  }

  heights <- vapply(nodes, function(x) attr(x, "height"), numeric(1))
  labs <- lapply(nodes, labels)

  # NOTE: this won't support nodes that have the same height
  # but that isn't possible, right?
  for (i in seq_along(heights)) {
    allXY$label[[which(allXY$y == heights[i])]] <- labs[[i]]
  }

  tidy_segments <- as.ggdend(d)$segments

  allTXT <- allXY[allXY$y == 0, ]

  allXY$members <- vapply(allXY$label, length, integer(1))
  allTXT$label <- as.character(allTXT$label)

  axis_x <- list( showgrid = F,
                  zeroline = F,
                  tickmode = "array",
                  tickvals = seq_along(allTXT$label),
                  ticktext = allTXT$label,
                  title = F)
  axis_y <- list( showgrid = T,
                  zeroline = F,
                  title = "Fraction of Regularization",
                  tickmode = "array",
                  tickvals = c(0,0.25,0.5,0.75,1),
                  ticktext = c("0%", "25%", "50%", "75%", "100%"),
                  range = c(-0.02, 1))

  # library(RColorBrewer)
  colors<-c(brewer.pal(8,"Set1"))

  if (dynamic == FALSE){
    if(show_clusters){
      labels <- get_cluster_labels(x, k.row = k.row, k.col = k.col, percent = percent, type = type)
      k <- nlevels(labels)
      if(has_percent){
        if (!is_percent_scalar(percent)) {
          crv_error(sQuote("percent"), " must be a scalar between 0 and 1 (inclusive).")
        }

        # if(percent == 0){ ## Don't bother showing boxes if percent is 0 and bail early
        #   return(invisible(x))
        # }
      } else {
        if (!is_integer_scalar(k)) {
          crv_error(sQuote("k"), " must be an integer scalar (vector of length 1).")
        }
        if ( k <= 0 ) {
          crv_error(sQuote("k"), " must be positive.")
        }
        if ( k > NROW(x$X) ) {
          crv_error(sQuote("k"), " cannot be more than the observations in the original data set (", NROW(x$X), ").")
        }

        percent <- x$row_fusions$cluster_membership %>%
          select(.data$GammaPercent, .data$NCluster) %>%
          filter(.data$NCluster <= k) %>%
          select(.data$GammaPercent) %>%
          summarize(percent = min(.data$GammaPercent)) %>%
          pull
      }

      cluster <- labels
      clustab <- table(cluster)[unique(cluster[labels(d)])]
      clustsum <- cumsum(clustab)
      m <- c(0, clustsum) + 0.5
      seg_x <- c(m,m[1],m[1])
      seg_y <- c(rep(percent,k+1),0,percent)
      seg_xend <- c(m,m[k+1],m[k+1])
      seg_yend <- c(rep(0,k+1),0,percent)
      x <- y <- xend <- yend <- NULL # for CMD check
      seg <- data.frame(percent=percent,x=seg_x,y=seg_y,xend=seg_xend,yend=seg_yend)

      dendro_static <- plot_ly(
        hoverinfo = "none") %>%
        add_segments(
          data = tidy_segments,
          x = ~x, y = ~y,
          xend = ~xend, yend = ~yend,
          color = I("black"),
          showlegend = FALSE)

      for (i in seq_along(clustsum)){
        dendro_static <- dendro_static %>%
          add_markers(
            data = allXY[allXY$y > 0 & allXY$y <= percent & allXY$x < m[i+1] & allXY$x > m[i], ],
            x = ~x, y = ~y, key = ~label, set = "A",
            name = "nodes",
            color = I(colors[(i-1)%%8+1]),
            text = ~paste0("members: ", members), hoverinfo = "text") %>%
          add_markers(
            data = allXY[allXY$y > percent, ],
            x = ~x, y = ~y, key = ~label, set = "A",
            name = "nodes",
            color = I("black"),
            text = ~paste0("members: ", members), hoverinfo = "text") %>%
          add_text(
            data = allTXT[allTXT$x < m[i+1] & allTXT$x > m[i], ], x = ~x, y = 0, text = "|", key = ~label, set = "A",
            textposition = "bottom center", textfont = list(size=25), color = I(colors[(i-1)%%8+1]),
          )
      }
      dendro_static %>%
        add_segments(data=seg, x=~x, y=~y, xend=~xend, yend=~yend, color = I("grey"))%>%
        plotly::layout(xaxis = axis_x, yaxis = axis_y) %>%
        hide_legend()
    } else{
      plot_ly(
        hoverinfo = "none") %>%
        add_segments(
          data = tidy_segments,
          x = ~x, y = ~y,
          xend = ~xend, yend = ~yend,
          color = I("black"),
          showlegend = FALSE) %>%
        add_markers(
          data = allXY[allXY$y > 0, ], x = ~x, y = ~y, key = ~label, set = "A",
          name = "nodes",
          color = I("black"),
          text = ~paste0("members: ", members), hoverinfo = "text") %>%
        add_text(
          data = allTXT, x = ~x, y = 0, text = "|", key = ~label, set = "A",
          textposition = "bottom center", textfont = list(size=25), color = I("black"),
        ) %>%
        plotly::layout(xaxis = axis_x, yaxis = axis_y) %>%
        hide_legend() %>%
        highlight(persistent = TRUE, dynamic = TRUE, off =NULL)
    }
    # invisible(x)
  } else {
    x <- y <- xend <- yend <- Regularization <- NULL # for CMD check
    dynamic_seg <- data.frame(percent = numeric(),
                              x = numeric(),
                              y = numeric(),
                              xend = numeric(),
                              yend = numeric())
    allXY_dynamic <- data.frame(x = numeric(),
                                y = numeric(),
                                label = integer(),
                                members = integer(),
                                color = character(),
                                Regularization = numeric())
    allTXT_dynamic <- data.frame(x = numeric(),
                                 y = numeric(),
                                 label = character(),
                                 Regularization = numeric())
    for (per in percent.seq){
      labels <- get_cluster_labels(x, percent = per, type = type)
      k <- nlevels(labels)
      cluster <- labels
      clustab <- table(cluster)[unique(cluster[labels(d)])]#[unique(cluster[tree$order])]
      clustsum <- cumsum(clustab)
      m <- c(0, clustsum) + 0.5
      seg_x <- c(m[1],m[1],m)
      seg_y <- c(0,per,rep(per,k+1))
      seg_xend <- c(m[k+1],m[k+1],m)
      seg_yend <- c(0,per,rep(0,k+1))
      x <- y <- xend <- yend <- NULL # for CMD check
      seg <- data.frame(Regularization=per,x=seg_x,y=seg_y,xend=seg_xend,yend=seg_yend)
      dynamic_seg <- rbind(dynamic_seg,seg)

      allXY_per <- allXY[allXY$y > 0,]
      allXY_per$label <- seq_along(allXY_per$label)

      allTXT_per <- allTXT

      for (i in seq_along(clustsum)){
        allXY_per[allXY_per$y <= per & allXY_per$x < m[i+1] & allXY_per$x > m[i], "color"] <- colors[(i-1)%%8+1]
        allTXT_per[allTXT_per$x < m[i+1] & allTXT_per$x > m[i], "color"] <- colors[(i-1)%%8+1]
      }
      allXY_per[allXY_per$y > per, "color"] <- "black"
      allXY_per$Regularization <- per
      allTXT_per$Regularization <- per

      allXY_dynamic <- rbind(allXY_dynamic, allXY_per)
      allTXT_dynamic <- rbind(allTXT_dynamic, allTXT_per)
    }

    count <- NULL # for CMD check
    c <- dynamic_seg %>% group_by(x,xend) %>% summarize(count = n())
    dynamic_seg<-arrange(merge(dynamic_seg,c),Regularization,desc(count),x,y)

    # avoid the situation when very short lines cannot be shown in the plotly
    ylength <- NULL # for CMD check
    tidy_segments_dynamic <- tidy_segments[,1:4] %>%
      mutate(ylength = abs(y-yend))
    tidy_segments_dynamic$ylength[tidy_segments_dynamic$ylength==0] <- 1

    dendro_dynamic <- plot_ly(
      hoverinfo = "none") %>%
      add_segments(
        data = tidy_segments_dynamic %>% arrange(desc(ylength)),
        x = ~x, y = ~y,
        xend = ~xend, yend = ~yend,
        color = I("black"),
        showlegend = FALSE)

    for (i in unique(allXY_dynamic$label)){
      dendro_dynamic <- dendro_dynamic %>%
        add_markers(
          data = allXY_dynamic[allXY_dynamic$label==i,], x = ~x, y = ~y,
          color = I(~color),
          frame = ~Regularization*100,
          text = ~paste0("members: ", members), hoverinfo = "text")
    }
    for (i in unique(allTXT_dynamic$x)){
      dendro_dynamic <- dendro_dynamic %>%
        add_text(
          data = allTXT_dynamic[allTXT_dynamic$x==i,], x = ~x, y = 0, text = "|", key = ~label, #set = "A",
          textposition = "bottom center", textfont = list(size=25), color = I(~color),
          frame = ~Regularization*100
        )
    }
    dendro_dynamic %>%
      add_segments(data=dynamic_seg, x=~x, y=~y, xend=~xend, yend=~yend, color = I("grey"),frame=~Regularization*100)%>%
      plotly::layout(#dragmode = "select",
        xaxis = axis_x, yaxis = axis_y) %>%
      hide_legend() %>%
      animation_slider(y=slider_y,currentvalue = list(prefix = "Regularization: ", suffix = "%"))
  }
}
