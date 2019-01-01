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
#' @importFrom ggplot2 ggplot aes_string geom_path geom_point geom_text guides theme
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
                      type = c("dendrogram", "path", "dynamic_path", "js", "interactive"),
                      axis = c("PC1", "PC2"),
                      dend.branch.width = 2,
                      dend.labels.cex = .6,
                      dend.ylab.cex = 1.2,
                      percent,
                      k,
                      percent.seq = seq(0, 1, length.out = 21),
                      max.nclust = 9,
                      min.nclust = 1) {

  type <- match.arg(type)
  switch(
    type,
    dendrogram = {
      carp_dendro_plot(x,
                       percent = percent,
                       k = k,
                       dend.branch.width = dend.branch.width,
                       dend.labels.cex = dend.labels.cex,
                       dend.ylab.cex = dend.ylab.cex,
                       ...)
    },
    path = {
      carp_path_plot(x,
                     axis = axis,
                     percent = percent,
                     k = k,
                     ...)
    },
    dynamic_path = {
      carp_dynamic_path_plot(x,
                             axis = axis,
                             percent.seq = percent.seq)
    },
    js = {
      carp_heatmaply(x,
                     ...,
                     percent = percent,
                     k = k)
    },
    interactive = {
      dots <- list(...)
      if ( length(dots) != 0 ){
        crv_error("Unknown arguments passed to ", sQuote("plot.CARP."))
      }

      shinyApp(
        ui = fluidPage(
          tags$style(
            type = "text/css",
            ".recalculating { opacity: 1.0; }"
          ),
          titlePanel("CARP Results [Convex Clustering]"),
          sidebarLayout(
            sidebarPanel(width = 2,
              uiOutput("choose_column1"),
              br(),
              uiOutput("choose_column2")
            ),
            mainPanel(
              tabsetPanel(
                tabPanel("Dynamic",
                  fluidRow(
                    column(width = 2,
                      sliderInput(
                        "regularization_movie",
                        "Amount of Regularization",
                        min = 0,
                        max = 1,
                        value = 0,
                        step = .03,
                        animate = animationOptions(interval = 1000,
                                                   loop = TRUE)
                      )
                    ),
                    column(width = 5, plotOutput("pcapathplot_movie", height = "700px")),
                    column(width = 5, plotOutput("dendplot_movie", height = "700px"))
                  )
                ),
                tabPanel("Static",
                  fluidRow(
                    column(width = 2,
                      sliderInput("regularization_static",
                                  "Number of Clusters",
                                  min = min.nclust,
                                  max = max.nclust,
                                  value = max.nclust,
                                  step = 1)),
                    column(width = 5, plotOutput("pcapathplot_static", height = "700px")),
                    column(width = 5, plotOutput("dendplot_static", height = "700px"))
                  )
                )
              )
            ))),
        server = function(input, output) {
          output$choose_column1 <- renderUI({
            selectInput("column1",
                        "Abscissa (X-Axis): ",
                        available_features(x),
                        selected = available_features(x)[1])
          })

          output$choose_column2 <- renderUI({
            selectInput("column2",
                        "Ordinate (Y-Axis): ",
                        available_features(x),
                        selected = available_features(x)[2])
          })

          output$dendplot_movie <- renderPlot({
            carp_dendro_plot(x,
                             percent = input$regularization_movie,
                             dend.branch.width = dend.branch.width,
                             dend.labels.cex = dend.labels.cex,
                             dend.ylab.cex = dend.ylab.cex,
                             ...)
          })

          output$pcapathplot_movie <- renderPlot({
            carp_path_plot(x,
                           axis = c(input$column1 %||% "PC1",
                                    input$column2 %||% "PC2"),
                           percent = input$regularization_movie,
                           show_clusters = FALSE,
                           repel_labels = FALSE,
                           label_size = 6) +
              theme(axis.title = element_text(size = 25),
                    axis.text  = element_text(size = 20))
          })

          output$dendplot_static <- renderPlot({
            carp_dendro_plot(x,
                             k = input$regularization_static,
                             base_colors = my_palette(n = input$regularization_static))
          })

          output$pcapathplot_static <- renderPlot({
            k <- input$regularization_static
            ## We need to manually construct a color map to keep the path plot
            ## and dendrogram colors lined-up
            my_colors <- my_palette(n = k)[order(unique(cutree(x$dendrogram, k = k)[x$dendrogram$order]))]

            carp_path_plot(x,
                           axis = c(input$column1 %||% "PC1",
                                    input$column2 %||% "PC2"),
                           k = k,
                           label_size = 6,
                           colors = my_colors) +
              theme(axis.title = element_text(size = 25),
                    axis.text  = element_text(size = 20))
          })
        }
      )
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

  g <- ggplot(mapping = aes_string(x = "V1", y = "V2", group = "Obs"))

  if (show_clusters) {
    g <- g + geom_path(data = plot_frame_full, aes_string(color = "final_cluster"), linejoin="round", size=1) +
             geom_point(data = plot_frame_final, aes_string(color = "final_cluster"), size = 4)

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
                        mapping = aes_string(label = "ObsLabel"),
                        size = label_size)
  } else {
    g + geom_text(data = plot_frame_init,
                  mapping = aes_string(label = "ObsLabel"),
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
#' @importFrom ggplot2 ggplot aes_string geom_path geom_point geom_text guides
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
         aes_string(x = "V1", y = "V2", group = "Obs")) +
    geom_path(linejoin = "round", color = "red", size = 1) +
    geom_point(data = plot_frame_first,
               aes_string(x = "FirstV1",
                          y = "FirstV2"),
               color = "black",
               size = I(4)) +
    geom_text(data = plot_frame_first,
              aes_string(x = "FirstV1",
                         y = "FirstV2",
                         label = "FirstObsLabel"),
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

