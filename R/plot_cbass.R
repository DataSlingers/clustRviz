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
#'            \item when \code{type == "js"}, \code{...} is forwarded to
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
                                "col.path",
                                "js",
                                "interactive"),
                       percent,
                       k.row,
                       k.col,
                       dend.branch.width = 2,
                       dend.labels.cex = .6,
                       dend.ylab.cex = 1.2,
                       heatrow.label.cex = 1,
                       heatcol.label.cex = 1,
                       margins = c(5, 5),
                       axis = c("PC1", "PC2")) {

  type <- match.arg(type)

  switch(
    type,
    row.dendrogram = {
      cbass_dendro_plot(x,
                        percent = percent,
                        k.row = k.row,
                        k.col = k.col,
                        dend.branch.width = dend.branch.width,
                        dend.labels.cex = dend.labels.cex,
                        type = "row",
                        ...)
    },
    col.dendrogram = {
      cbass_dendro_plot(x,
                        percent = percent,
                        k.row = k.row,
                        k.col = k.col,
                        dend.branch.width = dend.branch.width,
                        dend.labels.cex = dend.labels.cex,
                        dend.ylab.cex = dend.ylab.cex,
                        type = "col",
                        ...)
    },
    row.path = {
      cbass_path_plot(x,
                     axis = axis,
                     percent = percent,
                     k.row = k.row,
                     k.col = k.col,
                     ...,
                     type = "row")
    },
    col.path = {
      cbass_path_plot(x,
                      axis = axis,
                      percent = percent,
                      k.row = k.row,
                      k.col = k.col,
                      ...,
                      type = "col")
    },
    heatmap = {
      cbass_heatmap_plot(x,
                         ...,
                         percent = percent,
                         k.row = k.row,
                         k.col = k.col,
                         heatrow.label.cex = heatrow.label.cex,
                         heatcol.label.cex = heatcol.label.cex,
                         margins = margins)
    },
    js = {
      cbass_heatmaply(x,
                      ...,
                      percent = percent,
                      k.row = k.row,
                      k.col = k.col)
    },
    interactive = {
      dots <- list(...)
      if ( length(dots) != 0 ){
        crv_error("Unknown arguments passed to ", sQuote("plot.CARP."))
      }

      shinyApp(
        ui = fluidPage(tags$style(type = "text/css",
                                  ".recalculating { opacity: 1.0; }"),
                       titlePanel("CBASS Results [Convex BiClustering]"),
                       fluidRow(column(width = 2,
                                       sliderInput("regcent",
                                                   "Amount of Regularization",
                                                   min = 0,
                                                   max = 1,
                                                   value = 0.5,
                                                   step = 0.03,
                                                   animate = animationOptions(interval = 400,
                                                                              loop = TRUE))),
                                column(width = 10,
                                       plotOutput("heatmap",
                                                  height = "900px",
                                                  width = "1200px")))),
        server = function(input, output) {
          ## Calculate breaks and colors on the raw data so that they are consitent
          ## across frames. (If we use cbass_heatmap_plot's internal fitting, it will
          ## only look at the gradient for a single frame.)
          nbreaks     <- 50
          quant.probs <- seq(0, 1, length.out = nbreaks)
          breaks      <- unique(quantile(x$X[TRUE], probs = quant.probs))
          nbreaks     <- length(breaks)
          heatmap_col <- colorRampPalette(c("blue", "yellow"))(nbreaks - 1)

          output$heatmap <- renderPlot({
            cbass_heatmap_plot(x,
                               percent = input$regcent,
                               heatrow.label.cex = heatrow.label.cex,
                               heatcol.label.cex = heatcol.label.cex,
                               ...,
                               breaks = breaks,
                               heatmap_col = heatmap_col,
                               margins = margins)
          })
        }
      )
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
#' Render CBASS results via the heatmaply (interactive JS) package
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
