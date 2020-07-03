#' Visualize the results of Convex Clustering (\code{CARP})
#'
#' \code{plot.CARP} provides a range of ways to visualize the results of convex
#' clustering, including: \itemize{
#' \item A dendrogram, illustrating the nested cluster hierarchy inferred from
#'       the convex clustering solution path (\code{type = "dendrogram"});
#' \item A static path plot, showing the coalescence of the estimated cluster centroids
#'       at a fixed value of the regularization parameter is increased (\code{type = "path"});
#' \item A \code{\link[gganimate]{gganimate}} plot, showing the coalescence of the
#'       estimated cluster centroids as the regularization parameter is increased
#'       (\code{type = "dynamic_path"}); and
#' \item A \code{\link[shiny]{shiny}} app, which can display the clustering solutions
#'       as a "movie" or allow for interactive exploration (\code{type = "interactive"}).
#' }
#'
#' @param x An object of class \code{CARP} as returned by \code{\link{CARP}}
#' @param type A string indicating the type of visualization to show (see details above).
#' @param axis A character vector of length two indicating which features or principal
#'             components to use as the axes in the \code{type = "path"} visualization.
#'             Currently only features like \code{"PC1"} or \code{"PC2"} (indicating
#'             the first principal component projections) are supported.
#' @param percent A number between 0 and 1, giving the regularization level (as
#'                a fraction of the final regularization level used) at which to
#'                assign clusters in the static (\code{type = "dendrogram"} or \code{type = "path"})
#'                plots.
#' @param k An integer indicating the desired number of clusters to be displayed
#'          in the static plots. If no \code{CARP} iteration with exactly this
#'          many clusters is found, the first iterate with fewer than \code{k}
#'          clusters is used.
#' @param ... Additional arguments, which are handled differently for different
#'            values of \code{type}.\itemize{
#'            \item When \code{type} is \code{"path"}, \code{"dynamic_path"},
#'                  or \code{"interactive"}, the presence of
#'                  unknown arguments triggers an error;
#'            \item when \code{type == "dendrogram"}
#'                  \code{...} is forwarded to \code{\link[stats]{plot.dendrogram}}; and
#'            \item when \code{type == "js"}, \code{...} is forwarded to
#'                  \code{\link[heatmaply]{heatmaply}}.
#'            } See the documentation of the linked functions for details about
#'            additional supported arguments. \code{saveviz} passes arguments
#'            to the corresponding plot \code{type}.
#' @param max.nclust a positive integer. The maximum number of clusters
#' to display in the interactive plot.
#' @param min.nclust a positive value. The minimum number of clusters to
#' display in the interactive plot.
#' @param dend.branch.width Branch width for dendrogram plots (ignored for
#'        other plot types) - must be positive.
#' @param dend.labels.cex Label size for dendrogram plots (ignored for other plot
#'        types) - must be positive.
#' @param dend.ylab.cex Y-axis size for dendrogram plots (ignored for other plot
#'        types) - must be positive.
#' @return The value of the return type depends on the \code{type} argument:\itemize{
#'   \item if \code{type = "dendrogram"}, \code{x} is returned invisibly;
#'   \item if \code{type = "path"}, an object of class \code{\link[ggplot2]{ggplot}}
#'         which can be plotted directly (by invoking its print method) or modified
#'         further by the user is returned;
#'   \item if \code{type = "dynamic_path"}, an object of class \code{\link[gganimate:gganimate-package]{gganim}}
#'         is returned, and many be further manipulated by the user or plotted directly;
#'   \item if \code{type = "interactive"}, a \code{\link[shiny]{shiny}} app which can be activated
#'         by invoking its print method.
#' } \code{saveviz.CARP} always returns \code{file.name} invisibly.
#' @details The \code{\link{saveviz.CARP}} function provides a unified interface
#'          for exporting \code{CARP} visualizations to files. For all plots,
#'          at most one of \code{percent} and \code{k} may be supplied.
#' @importFrom shiny shinyApp fluidPage titlePanel tabsetPanel fluidRow
#' @importFrom shiny column plotOutput sliderInput uiOutput renderUI tags
#' @importFrom shiny selectInput animationOptions renderPlot tabPanel
#' @importFrom shiny br sidebarPanel mainPanel sidebarLayout
#' @importFrom stats as.dendrogram median
#' @importFrom ggplot2 ggplot aes geom_path geom_point geom_text guides theme
#' @importFrom ggplot2 element_text xlab ylab scale_color_manual
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr filter select distinct rename mutate left_join select_ %>%
#' @importFrom grDevices adjustcolor
#' @importFrom rlang %||%
#' @export
#' @rdname plot_carp
#' @examples
#' \dontrun{
#' carp_fit <- CARP(presidential_speech)
#' plot(carp_fit, type='interactive')
#' }
plot.CARP <- function(x,
                      ...,
                      type = c("dendrogram", "path", "js"),
                      dynamic = FALSE,
                      interactive = FALSE,
                      axis = c("PC1", "PC2"),
                      dend.branch.width = 2,
                      dend.labels.cex = .6,
                      dend.ylab.cex = 1.2,
                      percent,
                      k,
                      percent.seq = seq(0, 1, 0.01),
                      max.nclust = 9,
                      min.nclust = 1,
                      slider_y = -0.3) {

  type <- match.arg(type)
  switch(
    type,
    dendrogram = {
      if (interactive == FALSE) {
        if (dynamic == FALSE) {
          carp_dendro_plot(x,
                           percent = percent,
                           k = k,
                           dend.branch.width = dend.branch.width,
                           dend.labels.cex = dend.labels.cex,
                           dend.ylab.cex = dend.ylab.cex,
                           ...)
        }
      } else {
        carp_dendro_plotly(x,
                           dynamic = dynamic,
                           percent = percent,
                           k = k,
                           percent.seq = percent.seq,
                           slider_y = slider_y,
                           ...)
      }
    },
    path = {
      if (interactive == FALSE) {
        if (dynamic == FALSE) {
          carp_path_plot(x,
                         axis = axis,
                         percent = percent,
                         k = k,
                         ...)
        } else {
          carp_dynamic_path_plot(x,
                                 axis = axis,
                                 percent.seq = percent.seq)
        }
      } else {
        carp_path_plotly(x,
                         axis = axis,
                         dynamic = dynamic,
                         percent = percent,
                         k = k,
                         percent.seq = percent.seq,
                         ...)
      }
    },
    js = {
      carp_heatmaply(x,
                     ...,
                     percent = percent,
                     k = k)
    }
  )
}

#' @noRd
#' @importFrom rlang .data
#' @importFrom dplyr filter select left_join pull rename
#' @importFrom ggplot2 ggplot geom_path aes geom_point guides theme
#' @importFrom ggplot2 element_text xlab ylab scale_color_manual
#' @importFrom ggrepel geom_text_repel
carp_path_plot <- function(x,
                           ...,
                           axis,
                           percent,
                           k,
                           show_clusters = (n_args == 1L),
                           repel_labels = TRUE,
                           label_size = 3,
                           colors = NULL){

  dots <- list(...)
  if ( length(dots) != 0) {
    if (!is.null(names(dots))) {
      nm <- names(dots)
      nm <- nm[nzchar(nm)]
      crv_error("Unknown argument ", sQuote(nm[1]), " passed to ", sQuote("plot.CARP."))
    } else {
      crv_error("Unknown argument passed to ", sQuote("plot.CARP."))
    }
  }

  has_percent <- !missing(percent)
  has_k       <- !missing(k)

  n_args <- has_percent + has_k

  if(n_args > 1){
    crv_error("At most one of ", sQuote("percent"), " and ", sQuote("k"), " must be supplied.")
  }

  if (n_args == 0L){
    percent <- 1
    has_percent <- TRUE # We've now set the percent (whole path) for display
  }

  plot_frame_full <- get_feature_paths(x, axis) %>% rename(V1 = axis[1], V2 = axis[2])
  plot_frame_init <- plot_frame_full %>% filter(.data$Iter == min(.data$Iter))

  if (has_percent) {
    if (!is_percent_scalar(percent)) {
      crv_error(sQuote("percent"), " must be a scalar between 0 and 1 (inclusive).")
    }

    ## If percent == min(GammaPercent), keep some (unshrunken) data to plot
    ## This comes up in the default settings for the Shiny app or static when
    ## percent = 0
    plot_frame_full <- plot_frame_full %>% filter(.data$GammaPercent <= max(percent, min(.data$GammaPercent)))
  } else {
    # Get the first iteration at which we have k (or fewer) clusters
    # to avoid plotting "beyond" what we want

    if (!is_integer_scalar(k)) {
      crv_error(sQuote("k"), " must be an integer scalar (vector of length 1).")
    }
    if ( k <= 0 ) {
      crv_error(sQuote("k"), " must be positive.")
    }
    if ( k > NROW(x$X) ) {
      crv_error(sQuote("k"), " cannot be more than the observations in the original data set (", NROW(x$X), ").")
    }

    iter_first_k <- plot_frame_full %>% select(.data$Iter, .data$NCluster) %>%
                                        filter(.data$NCluster <= k) %>%
                                        summarize(iter_first_k = min(.data$Iter)) %>%
                                        pull

    plot_frame_full <- plot_frame_full %>% filter(.data$Iter <= iter_first_k)
  }


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

    if (!is.null(colors)) {
      g <- g + scale_color_manual(values = colors)
    }
  } else {
    ## We always show the centroids, even if we don't color things
    g <- g + geom_path(data = plot_frame_full, color = "red", linejoin = "round", size = 1) +
             geom_point(data = plot_frame_final, color = "red", size = 4)
  }

  g <- g + geom_point(data = plot_frame_init, color="black", size = 2) +
           guides(color = FALSE, size = FALSE) +
           theme(axis.title = element_text(size = 15),
                 axis.text  = element_text(size = 10)) +
           xlab(axis[1]) +
           ylab(axis[2])

  if (repel_labels) {
    g + geom_text_repel(data = plot_frame_init,
                        mapping = aes(label = .data$ObsLabel),
                        size = label_size)
  } else {
    g + geom_text(data = plot_frame_init,
                  mapping = aes(label = .data$ObsLabel),
                  size = label_size)
  }
}

#' @noRd
#' @importFrom rlang .data
#' @importFrom dplyr filter select summarize pull
#' @importFrom stats as.dendrogram
#' @importFrom dendextend set
#' @importFrom grDevices adjustcolor
carp_dendro_plot <- function(x,
                             percent,
                             k,
                             dend.branch.width = 2,
                             dend.labels.cex = .6,
                             dend.ylab.cex = 1.2,
                             show_clusters = (n_args == 1L),
                             base_colors = c("grey", "black"),
                             ...){

  if(dend.branch.width <= 0){
    crv_error(sQuote("dend.branch.width"), " must be positive.")
  }

  if(dend.labels.cex <= 0){
    crv_error(sQuote("dend.labels.cex"), " must be positive.")
  }

  has_percent <- !missing(percent)
  has_k       <- !missing(k)
  n_args      <- has_percent + has_k

  if(n_args > 1L){
    crv_error("At most one of ", sQuote("percent"), " and ", sQuote("k"), " may be supplied.")
  }

  as.dendrogram(x) %>%
    set("branches_lwd", dend.branch.width) %>%
    set("labels_cex", dend.labels.cex) %>%
    plot(ylab = "Fraction of Regularization",
         cex.lab = dend.ylab.cex, yaxt = "n",
         ...)

  axis(2, at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0%", "25%", "50%", "75%", "100%"), las=2)
  if(show_clusters){

    if(has_percent){
      if (!is_percent_scalar(percent)) {
        crv_error(sQuote("percent"), " must be a scalar between 0 and 1 (inclusive).")
      }

      if(percent == 0){ ## Don't bother showing boxes if percent is 0 and bail early
        return(invisible(x))
      }

      k <- get_feature_paths(x, features = character()) %>% filter(.data$GammaPercent <= percent) %>%
                                                            select(.data$NCluster) %>%
                                                            summarize(NCluster = min(.data$NCluster)) %>%
                                                            pull
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
    }

    my.cols <- adjustcolor(base_colors, alpha.f = .2)
    my.rect.hclust(as.hclust(x), k = k, border = 2, my.col.vec = my.cols, lwd = 3)
  }

  invisible(x)
}


#' @noRd
#' @importFrom rlang .data
#' @importFrom dplyr select filter rename left_join mutate bind_rows
#' @importFrom ggplot2 ggplot aes geom_path geom_point geom_text guides
#' @importFrom ggplot2 theme element_text xlab ylab
#' @importFrom gganimate transition_manual
carp_dynamic_path_plot <- function(x, axis, percent.seq){
  ## TODO - Combine this and carp_path_plot as much as possible
  plot_frame_full <- get_feature_paths(x, axis) %>% rename(V1 = axis[1], V2 = axis[2])

  plot_frame_first <- plot_frame_full %>% filter(.data$Iter == min(.data$Iter)) %>%
                                          select(.data$Obs, .data$V1, .data$V2, .data$ObsLabel) %>%
                                          rename(FirstV1 = .data$V1,
                                                 FirstV2 = .data$V2,
                                                 FirstObsLabel = .data$ObsLabel)

  plot_frame_animation <- bind_rows(lapply(percent.seq, function(pct){
    ## Make a list of things to plot at each "percent" and then combine
    plot_frame_full %>% filter(.data$GammaPercent <= pct) %>%
                        mutate(percent = pct)
  }))

  ggplot(plot_frame_animation,
         aes(x = .data$V1, y = .data$V2, group = .data$Obs)) +
    geom_path(linejoin = "round", color = "red", size = 1) +
    geom_point(data = plot_frame_first,
               aes(x = .data$FirstV1,
                   y = .data$FirstV2),
               color = "black",
               size = I(4)) +
    geom_text(data = plot_frame_first,
              aes(x = .data$FirstV1,
                  y = .data$FirstV2,
                  label = .data$FirstObsLabel),
              size = I(6)) +
    guides(color = FALSE, size = FALSE) +
    theme(axis.title = element_text(size = 25),
          axis.text = element_text(size = 20)) +
    xlab(axis[1]) + ylab(axis[2]) +
    transition_manual(.data$percent)
}

#' @noRd
#' Render CBASS results via the heatmaply (interactive JS) package
#' @importFrom heatmaply heatmaply
carp_heatmaply <- function(x,
                           ...,
                           percent,
                           k){

  has_percent <- !missing(percent)
  has_k       <- !missing(k)

  n_args <- has_percent + has_k

  if(n_args >= 2){
    crv_error("At most one of ", sQuote("percent,"), " and ",
              sQuote("k"), " may be supplied.")
  }

  if(n_args == 0){
    U     <- get_clustered_data(x, percent = 0, refit = TRUE)

    heatmaply(U,
              Rowv = as.hclust(x),
              dendrogram = "row",
              ...)
  } else {
    U <- get_clustered_data(x, percent = percent, k = k, refit = TRUE)
    k <- nlevels(get_cluster_labels(x, percent = percent, k = k))

    heatmaply(U,
              Rowv = as.hclust(x),
              k_row = k,
              dendrogram = "row",
              ...)
  }
}

#' @noRd
#' @importFrom rlang .data
#' @importFrom dplyr filter select left_join pull rename
#' @importFrom plotly add_markers add_paths add_text
carp_path_plotly <- function(x,
                           ...,
                           dynamic = FALSE,
                           axis,
                           percent,
                           percent.seq,
                           k,
                           show_clusters = (n_args == 1L)){

  dots <- list(...)
  if ( length(dots) != 0) {
    if (!is.null(names(dots))) {
      nm <- names(dots)
      nm <- nm[nzchar(nm)]
      crv_error("Unknown argument ", sQuote(nm[1]), " passed to ", sQuote("plot.CARP."))
    } else {
      crv_error("Unknown argument passed to ", sQuote("plot.CARP."))
    }
  }

  has_percent <- !missing(percent)
  has_k       <- !missing(k)

  n_args <- has_percent + has_k

  if(n_args > 1){
    crv_error("At most one of ", sQuote("percent"), " and ", sQuote("k"), " must be supplied.")
  }

  if(n_args >= 1 & dynamic == TRUE){
    crv_error("Can't set ", sQuote("percent"), " and ", sQuote("k"), " for dynamic plots")
  }

  if (n_args == 0L){
    percent <- 1
    has_percent <- TRUE # We've now set the percent (whole path) for display
  }

  plot_frame_full <- get_feature_paths(x, axis) %>% rename(V1 = axis[1], V2 = axis[2])
  plot_frame_init <- plot_frame_full %>% filter(.data$Iter == min(.data$Iter))

  if (has_percent) {
    if (!is_percent_scalar(percent)) {
      crv_error(sQuote("percent"), " must be a scalar between 0 and 1 (inclusive).")
    }

    ## If percent == min(GammaPercent), keep some (unshrunken) data to plot
    ## This comes up in the default settings for the static when percent = 0
    plot_frame_full <- plot_frame_full %>% filter(.data$GammaPercent <= max(percent, min(.data$GammaPercent)))
  } else {
    # Get the first iteration at which we have k (or fewer) clusters
    # to avoid plotting "beyond" what we want

    if (!is_integer_scalar(k)) {
      crv_error(sQuote("k"), " must be an integer scalar (vector of length 1).")
    }
    if ( k <= 0 ) {
      crv_error(sQuote("k"), " must be positive.")
    }
    if ( k > NROW(x$X) ) {
      crv_error(sQuote("k"), " cannot be more than the observations in the original data set (", NROW(x$X), ").")
    }

    iter_first_k <- plot_frame_full %>% select(.data$Iter, .data$NCluster) %>%
      filter(.data$NCluster <= k) %>%
      summarize(iter_first_k = min(.data$Iter)) %>%
      pull

    plot_frame_full <- plot_frame_full %>% filter(.data$Iter <= iter_first_k)
  }


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

      for (i in 1:length(plot_frame_init$ObsLabel)){
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

      for (i in 1:length(plot_frame_init$ObsLabel)){
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

    for (i in 1:length(plot_frame_init$ObsLabel)){
      path_dynamic <- path_dynamic %>%
        add_paths(
          data = plot_frame_animation[plot_frame_animation$Obs==i,],
          color = I("red"),
          ids = ~Obs,
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
        size = 2,
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
#' @importFrom dplyr filter select summarize pull
#' @importFrom stats as.dendrogram
#' @importFrom dendextend get_nodes_xy
#' @importFrom plotly add_segments add_maekers add_text
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tibble as_tibble
carp_dendro_plotly <- function(x,
                             dynamic = FALSE,
                             percent,
                             percent.seq,
                             k,
                             show_clusters = (n_args == 1L),
                             slider_y = -0.3,
                             ...
){

  has_percent <- !missing(percent)
  has_k       <- !missing(k)
  n_args      <- has_percent + has_k

  if(n_args > 1L){
    crv_error("At most one of ", sQuote("percent"), " and ", sQuote("k"), " may be supplied.")
  }

  d <- as.dendrogram(x)
  # get x/y locations of every node in the tree
  m <- dendextend::get_nodes_xy(d)
  colnames(m) <- c("x", "y")
  allXY <- tibble::as_tibble(m)
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

  heights <- sapply(nodes, function(x) attr(x, "height"))
  labs <- lapply(nodes, labels)

  # NOTE: this won't support nodes that have the same height
  # but that isn't possible, right?
  for (i in seq_along(heights)) {
    allXY$label[[which(allXY$y == heights[i])]] <- labs[[i]]
  }

  tidy_segments <- dendextend::as.ggdend(d)$segments

  allTXT <- allXY[allXY$y == 0, ]

  allXY$members <- sapply(allXY$label, length)
  allTXT$label <- as.character(allTXT$label)

  axis_x <- list( showgrid = F,
                  zeroline = F,
                  tickmode = "array",
                  tickvals = 1:length(allTXT$label),
                  ticktext = allTXT$label,
                  title = F)
  axis_y <- list( showgrid = T,
                  zeroline = F,
                  title = "Fraction of Regularization",
                  tickmode = "array",
                  tickvals = c(0,0.25,0.5,0.75,1),
                  ticktext = c("0%", "25%", "50%", "75%", "100%"),
                  range = c(-0.02, 1))

  colors<-c(brewer.pal(8,"Set1"))

  if (dynamic == FALSE){
    if(show_clusters){

      if(has_percent){
        if (!is_percent_scalar(percent)) {
          crv_error(sQuote("percent"), " must be a scalar between 0 and 1 (inclusive).")
        }

        # if(percent == 0){ ## Don't bother showing boxes if percent is 0 and bail early
        #   return(invisible(x))
        # }

        k <- get_feature_paths(x, features = character()) %>% filter(.data$GammaPercent <= percent) %>%
          select(.data$NCluster) %>%
          summarize(NCluster = min(.data$NCluster)) %>%
          pull
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

        percent <- get_feature_paths(x, features = character()) %>% filter(.data$NCluster == k) %>%
          select(.data$GammaPercent) %>%
          summarize(NCluster = mean(.data$GammaPercent)) %>%
          pull
      }

      tree <- as.hclust(x)
      cluster <- stats::cutree(tree, k = k)
      clustab <- table(cluster)[unique(cluster[tree$order])]
      clustsum <- cumsum(clustab)
      m <- c(0, clustsum) + 0.5
      seg_x <- c(m,m[1],m[1])
      seg_y <- c(rep(percent,k+1),0,percent)
      seg_xend <- c(m,m[k+1],m[k+1])
      seg_yend <- c(rep(0,k+1),0,percent)
      seg <- data.frame(percent=percent,x=seg_x,y=seg_y,xend=seg_xend,yend=seg_yend)

      dendro_static <- plot_ly(
        hoverinfo = "none") %>%
        add_segments(
          data = tidy_segments,
          x = ~x, y = ~y,
          xend = ~xend, yend = ~yend,
          color = I("black"),
          showlegend = FALSE)

      for (i in 1:length(clustsum)){
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
      k <- get_feature_paths(x, features = character()) %>% filter(.data$GammaPercent <= per) %>%
        select(.data$NCluster) %>%
        summarize(NCluster = min(.data$NCluster)) %>%
        pull
      tree <- as.hclust(x)
      cluster <- stats::cutree(tree, k = k)
      clustab <- table(cluster)[unique(cluster[tree$order])]
      clustsum <- cumsum(clustab)
      m <- c(0, clustsum) + 0.5
      seg_x <- c(m[1],m[1],m)
      seg_y <- c(0,per,rep(per,k+1))
      seg_xend <- c(m[k+1],m[k+1],m)
      seg_yend <- c(0,per,rep(0,k+1))
      seg <- data.frame(Regularization=per,x=seg_x,y=seg_y,xend=seg_xend,yend=seg_yend)
      dynamic_seg <- rbind(dynamic_seg,seg)

      allXY_per <- allXY[allXY$y > 0,]
      allXY_per$label <- 1:length(allXY_per$label)

      allTXT_per <- allTXT

      for (i in 1:length(clustsum)){
        allXY_per[allXY_per$y <= per & allXY_per$x < m[i+1] & allXY_per$x > m[i], "color"] <- colors[(i-1)%%8+1]
        allTXT_per[allTXT_per$x < m[i+1] & allTXT_per$x > m[i], "color"] <- colors[(i-1)%%8+1]
      }
      allXY_per[allXY_per$y > per, "color"] <- "black"
      allXY_per$Regularization <- per
      allTXT_per$Regularization <- per

      allXY_dynamic <- rbind(allXY_dynamic, allXY_per)
      allTXT_dynamic <- rbind(allTXT_dynamic, allTXT_per)
    }

    c <- dynamic_seg %>% group_by(x,xend) %>% summarise(count = n())
    dynamic_seg<-arrange(merge(dynamic_seg,c),Regularization,desc(count),x,y)

    # avoid the situation when very short lines cannot be shown in the plotly
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

