#' Visualize the results of Convex Clustering (\code{CARP})
#'
#' \code{plot.CARP} provides a range of ways to visualize the results of convex
#' clustering, including: \itemize{
#' \item A dendrogram, illustrating the nested cluster hierarchy inferred from
#'       the convex clustering solution path (\code{type = "dendrogram"});
#' \item A path plot, showing the coalescence of the estimated cluster centroids
#'       as the regularization parameter is increased (\code{type = "path"}); and
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
#'                assign clusters in the static (\code{type = "dendrogram" or \code{type = "path"})
#'                plots.
#' @param k An integer indicating the desired number of clusters to be displayed
#'          in the static plots. If no \code{CARP} iteration with exactly this
#'          many clusters is found, the first iterate with fewer than \code{k}
#'          clusters is used.
#' @param show_clusters A Boolean value indicating whether the cluster assignments
#'                      are indicated in the static plot types
#' @param ... Additional arguments. Currently an error when \code{type != "dendrogram"}
#'            and passed to \code{\link[stats]{plot.dendrogram}} when \code{type =
#'            "dendrogram"}.
#' @param max.nclust a positive integer. The maximum number of clusters
#' to display in the interactive plot.
#' @param min.nclust a positive value. The minimum number of clusters to
#' display in the interactive plot.
#' @param dend.branch.width a positive number. Line width on dendrograms.
#' @param dend.labels.cex a positive number. Label size on dendrograms.
#' @return The value of the return type depends on the \code{type} argument:\itemize{
#'   \item if \code{type = "dendrogram"}, \code{x} is returned invisibly;
#'   \item if \code{type = "path"}, an object of class \code{\link[ggplot2]{ggplot}}
#'         which can be plotted directly (by invoking its print method) or modified
#'         further by the user is returned;
#'   \item if \code{type = "interactive"}, a \code{shiny} app which can be activated
#'         by invoking its print method.
#' }
#' @details The \code{\link{saveviz.CARP}} function provides a unified interface
#'          for exporting \code{CARP} visualizations to files. For all plots,
#'          at most one of \code{percent} and \code{k} may be supplied.
#' @importFrom shiny shinyApp fluidPage titlePanel tabsetPanel fluidRow
#' @importFrom shiny column plotOutput sliderInput uiOutput renderUI tags
#' @importFrom shiny checkboxGroupInput animationOptions renderPlot
#' @importFrom stats as.dendrogram median
#' @importFrom ggplot2 ggplot aes geom_path geom_point geom_text guides theme
#' @importFrom ggplot2 element_text xlab ylab scale_color_manual
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr filter select distinct rename mutate left_join select_ %>%
#' @importFrom grDevices adjustcolor
#' @importFrom RColorBrewer brewer.pal
#' @export
#' @examples
#' \dontrun{
#' carp_fit <- CARP(presidential_speech)
#' plot(carp_fit, type='interactive')
#' }
plot.CARP <- function(x,
                      ...,
                      type = c("dendrogram", "path", "interactive"),
                      axis = c("PC1", "PC2"),
                      dend.branch.width = 2,
                      dend.labels.cex = .6,
                      percent,
                      k,
                      max.nclust = 9,
                      min.nclust = 1,
                      show_clusters) {

  LambdaPercent <- NULL
  Iter <- NULL
  V1 <- NULL
  V2 <- NULL
  Obs <- NULL
  ObsLabel <- NULL
  LambdaPercent <- NULL
  NCluster <- NULL
  Iter <- NULL
  Obs <- NULL
  Cluster <- NULL
  Var1 <- NULL
  Var2 <- NULL
  MaxVar1 <- NULL
  MaxVar2 <- NULL
  ObsLabel <- NULL
  PlotCluster <- NULL

  type <- match.arg(type)
  switch(
    type,
    dendrogram = {
      x$carp.dend %>%
        stats::as.dendrogram() %>%
        dendextend::set("branches_lwd", dend.branch.width) %>%
        dendextend::set("labels_cex", dend.labels.cex) %>%
        plot(ylab = "Amount of Regularization")
    },
    path = {
      carp_path_plot(x,
                     axis = axis,
                     percent = percent,
                     k = k,
                     show_clusters = show_clusters,
                     ...)
    },
    interactive = {
      dots <- list(...)
      if ( length(dots) != 0 ){
        stop("Unknown arguments passed to ", sQuote("plot.CARP."))
      }

      shiny::shinyApp(
        ui = shiny::fluidPage(
          shiny::tags$style(
            type = "text/css",
            ".recalculating { opacity: 1.0; }"
          ),
          shiny::titlePanel("Clustering Example"),
          shiny::tabsetPanel(
            shiny::tabPanel(
              "Movie",
              shiny::fluidRow(
                shiny::column(
                  2,
                  shiny::sliderInput(
                    "regcent_movie",
                    "Amount of Regularization",
                    min = 0,
                    max = 1,
                    value = 0,
                    # step=.07,animate = animationOptions(interval=700,loop=T)
                    step = .03, animate = shiny::animationOptions(interval = 1000, loop = T)
                  ),
                  shiny::uiOutput("choose_columns")
                ),
                shiny::column(
                  5,
                  shiny::plotOutput("pcapathplot_movie", height = "700px") # ,width = '900px')
                ),
                shiny::column(
                  5,
                  shiny::plotOutput("dendplot_movie", height = "700px")
                )
              )
            ),
            shiny::tabPanel(
              "Static",
              shiny::fluidRow(
                shiny::column(
                  2,
                  shiny::sliderInput("regcent_static",
                    "Number of Clusters",
                    min = min.nclust,
                    max = max.nclust,
                    value = max.nclust,
                    step = 1
                  ),
                  shiny::uiOutput("choose_columns_static")
                ),
                shiny::column(
                  5,
                  shiny::plotOutput("pcapathplot_static", height = "700px") # ,width = '900px')
                ),
                shiny::column(
                  5,
                  shiny::plotOutput("dendplot_static", height = "700px")
                )
              )
            )
          )
        ),
        server = function(input, output) {

          # Drop-down selection box for which data set
          # Check boxes
          output$choose_columns <- shiny::renderUI({

            # Get the data set with the appropriate name
            colnames <- paste("PC", 1:4, sep = "")

            # Create the checkboxes and select them all by default
            shiny::checkboxGroupInput("columns", "Choose columns",
              choices = colnames,
              selected = colnames[1:2]
            )
          })
          # Check boxes
          output$choose_columns_static <- shiny::renderUI({
            # Get the data set with the appropriate name
            colnames <- paste("PC", 1:4, sep = "")

            # Create the checkboxes and select them all by default
            shiny::checkboxGroupInput("columns_static", "Choose columns",
              choices = colnames,
              selected = colnames[1:2]
            )
          })


          output$dendplot_movie <- shiny::renderPlot({
            x$carp.cluster.path.vis %>%
              dplyr::filter(LambdaPercent <= input$regcent_movie) %>%
              dplyr::select(NCluster) %>%
              unlist() %>%
              unname() %>%
              min() -> ncl
            x$carp.dend %>%
              stats::as.dendrogram() %>%
              dendextend::set("branches_lwd", 2) %>%
              dendextend::set("labels_cex", .6) %>%
              plot(ylab = "Amount of Regularization", cex.lab = 1.5)
            my.cols <- grDevices::adjustcolor(c("grey", "black"), alpha.f = .2)
            my.rect.hclust(x$carp.dend, k = ncl, border = 2, my.col.vec = my.cols, lwd = 3)
          })

          output$pcapathplot_movie <- shiny::renderPlot({
            min.iter <- 5
            rename.list <- list(
              Obs = "Obs",
              Cluster = "Cluster",
              Iter = "Iter",
              ObsLabel = "ObsLabel",
              NCluster = "NCluster",
              Var1 = input$columns[1],
              Var2 = input$columns[2]
            )
            if (input$regcent_movie == 0) {
              cl.iter <- min.iter
            } else {
              x$carp.cluster.path.vis %>%
                dplyr::filter(LambdaPercent <= input$regcent_movie) %>%
                dplyr::select(Iter) %>%
                unlist() %>%
                unname() %>%
                max() -> cl.iter
            }
            rename.list <- list(
              Obs = "Obs",
              Cluster = "Cluster",
              Iter = "Iter",
              ObsLabel = "ObsLabel",
              NCluster = "NCluster",
              Var1 = input$columns[1],
              Var2 = input$columns[2]
            )
            x$carp.cluster.path.vis %>%
              dplyr::select_(.dots = rename.list) -> x$carp.cluster.path.vis.rename

            x$carp.cluster.path.vis.rename %>%
              dplyr::filter(Iter == cl.iter) %>%
              dplyr::select(Obs, Cluster, Var1, Var2) %>%
              dplyr::rename(
                MaxVar1 = Var1,
                MaxVar2 = Var2
              ) -> cl.assgn
            x$carp.cluster.path.vis.rename %>%
              dplyr::filter(Iter <= cl.iter) %>%
              dplyr::left_join(
                cl.assgn,
                by = c("Obs")
              ) -> tmp
            tmp %>%
              dplyr::filter(Iter > min.iter) %>%
              ggplot2::ggplot(aes(x = Var1, y = Var2, group = Obs)) +
              ggplot2::geom_path(
                ggplot2::aes(x = Var1, y = Var2),
                linejoin = "round",
                color = "red"
              ) +
              ggplot2::geom_point(
                ggplot2::aes(x = MaxVar1, y = MaxVar2),
                data = tmp %>% dplyr::filter(Iter == 1),
                color = "red",
                size = I(4)
              ) +
              ggplot2::geom_point(
                aes(x = Var1, y = Var2),
                data = tmp %>% dplyr::filter(Iter == 1),
                color = "black",
                size = I(4)
              ) +
              ggplot2::geom_text(
                ggplot2::aes(x = Var1, y = Var2, label = ObsLabel),
                size = I(6),
                data = tmp %>% dplyr::filter(Iter == 1)
              ) +
              ggplot2::guides(color = FALSE, size = FALSE) +
              ggplot2::theme(axis.title = ggplot2::element_text(size = 25)) +
              ggplot2::theme(axis.text = ggplot2::element_text(size = 20)) +
              ggplot2::xlab(input$columns[1]) +
              ggplot2::ylab(input$columns[2])
          })
          output$dendplot_static <- shiny::renderPlot({
            x$carp.dend %>%
              stats::as.dendrogram() %>%
              dendextend::set("branches_lwd", 2) %>%
              dendextend::set("labels_cex", .6) %>%
              plot(ylab = "Amount of Regularization", cex.lab = 1.5)
            my.cols <- grDevices::adjustcolor(RColorBrewer::brewer.pal(n = input$regcent_static, "Set1"), alpha.f = .2)
            my.rect.hclust(x$carp.dend, k = input$regcent_static, border = 2, my.col.vec = my.cols, lwd = 3)
          })
          output$pcapathplot_static <- shiny::renderPlot({
            ncl <- input$regcent_static
            my.cols <- grDevices::adjustcolor(RColorBrewer::brewer.pal(n = ncl, "Set1"))[order(unique(stats::cutree(x$carp.dend, k = ncl)[x$carp.dend$order]))]
            x$carp.cluster.path.vis %>%
              dplyr::distinct(Iter, NCluster) %>%
              dplyr::filter(NCluster == ncl) %>%
              dplyr::select(Iter) %>%
              unlist() %>%
              unname() %>%
              stats::median() -> cl.iter
            cl.iter <- floor(cl.iter)
            rename.list <- list(
              Obs = "Obs",
              Cluster = "Cluster",
              Iter = "Iter",
              ObsLabel = "ObsLabel",
              NCluster = "NCluster",
              Var1 = input$columns_static[1],
              Var2 = input$columns_static[2]
            )
            x$carp.cluster.path.vis %>%
              dplyr::select_(.dots = rename.list) -> carp.cluster.path.vis.rename

            carp.cluster.path.vis.rename %>%
              dplyr::filter(Iter == cl.iter) %>%
              dplyr::select(Obs, Cluster, Var1, Var2) %>%
              dplyr::rename(
                MaxVar1 = Var1,
                MaxVar2 = Var2,
                PlotCluster = Cluster
              ) %>%
              dplyr::mutate(
                PlotCluster = as.factor(PlotCluster)
              ) -> cl.assgn

            carp.cluster.path.vis.rename %>%
              dplyr::filter(Iter <= cl.iter) %>%
              dplyr::left_join(
                cl.assgn,
                by = c("Obs")
              ) -> tmp
            tmp %>%
              dplyr::filter(Iter > 50) %>%
              ggplot2::ggplot(aes(x = Var1, y = Var2, group = Obs)) +
              ggplot2::geom_path(
                ggplot2::aes(x = Var1, y = Var2, color = PlotCluster),
                linejoin = "round"
              ) +
              ggplot2::geom_point(
                ggplot2::aes(x = MaxVar1, y = MaxVar2, color = PlotCluster),
                data = tmp %>% dplyr::filter(Iter == 1),
                size = I(4)
              ) +
              ggplot2::geom_point(
                ggplot2::aes(x = Var1, y = Var2),
                data = tmp %>% dplyr::filter(Iter == 1),
                color = "black",
                size = I(4)
              ) +
              ggplot2::geom_text(
                aes(x = Var1, y = Var2, label = ObsLabel),
                size = I(6),
                data = tmp %>% dplyr::filter(Iter == 1)
              ) +
              ggplot2::scale_color_manual(values = my.cols) +
              ggplot2::guides(color = FALSE, size = FALSE) +
              ggplot2::theme(axis.title = ggplot2::element_text(size = 25)) +
              ggplot2::theme(axis.text = ggplot2::element_text(size = 20)) +
              ggplot2::xlab(input$columns_static[1]) +
              ggplot2::ylab(input$columns_static[2])
          })
        }
        # End Shiny App
      )
    }
  )
}

#' @noRd
#' @importFrom rlang .data
#' @importFrom dplyr filter select left_join pull
#' @importFrom ggplot2 ggplot geom_path aes geom_point guides theme element_text xlab ylab
#' @importFrom ggrepel geom_text_repel
carp_path_plot <- function(x, ..., axis, percent, k, show_clusters){

  dots <- list(...)
  if ( length(dots) != 0) {
    if (!is.null(names(dots))) {
      nm <- names(dots)
      nm <- nm[nzchar(nm)]
      stop("Unknown argument ", sQuote(nm[1]), " passed to ", sQuote("plot.CARP."))
    } else {
      stop("Unknown argument passed to ", sQuote("plot.CARP."))
    }
  }

  has_percent <- !missing(percent)
  has_k       <- !missing(k)

  n_args <- has_percent + has_k

  if(n_args > 1){
    stop("At most one of ", sQuote("percent"), " and ", sQuote("k"), " must be supplied.")
  }

  if (missing(show_clusters)) {
    show_clusters <- (n_args == 1L)
  }

  if (!is_logical_scalar(show_clusters)) {
    stop(sQuote("show_clusters"), " must be either TRUE or FALSE.")
  }

  if (show_clusters && (n_args != 1)) {
    stop("Exactly one of ", sQuote("percent"), " and ", sQuote("k"), " must be supplied if ", sQuote("show_clusters"), " is TRUE.")
  }

  if (n_args == 0L){
    percent <- 1
    has_percent <- TRUE # We've now set the percent (whole path) for display
  }

  plot_cols <- c(
    axis,
    "Iter",
    "Obs",
    "Cluster",
    "Lambda",
    "ObsLabel",
    "NCluster",
    "LambdaPercent"
  )

  if (any(plot_cols %not.in% colnames(x$carp.cluster.path.vis))) {
    missing_col <- plot_cols[which(plot_cols %not.in% colnames(x$carp.cluster.path.vis))][1]
    stop(sQuote(missing_col), " is not available for plotting.")
  }

  plot_frame_full <- x$carp.cluster.path.vis %>% select(plot_cols) %>%
                                                 filter(.data$Iter > x$burn.in)
  names(plot_frame_full)[1:2] <- c("V1", "V2")

  if (has_percent) {
    if (!is_percent_scalar(percent)) {
      stop(sQuote("percent"), " must be a scalar between 0 and 1 (inclusive).")
    }

    plot_frame_full <- plot_frame_full %>% filter(.data$LambdaPercent <= percent)
  } else {
    # Get the first iteration at which we have k (or fewer) clusters
    # to avoid plotting "beyond" what we want

    if (!is_integer_scalar(k)) {
      stop(sQuote("k"), " must be an integer scalar (vector of length 1).")
    }
    if ( k <= 0 ) {
      stop(sQuote("k"), " must be positive.")
    }
    if ( k > NROW(x$X) ) {
      stop(sQuote("k"), " cannot be more than the observations in the original data set (", NROW(x$X), ").")
    }

    iter_first_k <- plot_frame_full %>% select(.data$Iter, .data$NCluster) %>%
                                        filter(.data$NCluster <= k) %>%
                                        summarize(iter_first_k = min(.data$Iter)) %>%
                                        pull

    plot_frame_full <- plot_frame_full %>% filter(.data$Iter <= iter_first_k)
  }

  plot_frame_init  <- plot_frame_full %>% filter(.data$Iter == min(.data$Iter))
  plot_frame_final <- plot_frame_full %>% filter(.data$Iter == max(.data$Iter)) %>%
                                          mutate(final_cluster = factor(.data$Cluster))

  plot_frame_full <- left_join(plot_frame_full,
                               plot_frame_final %>% select(.data$Obs, .data$final_cluster),
                               by = "Obs")


  ## FIXME -- It looks like we don't actually have full fusion in `plot_frame_final`
  ##          (even in points which should be in the same cluster...)

  g <- ggplot(mapping = aes(x = V1, y = V2, group = Obs))

  if (show_clusters) {
    g <- g + geom_path(data = plot_frame_full, aes(color = final_cluster), linejoin="round", size=1) +
             geom_point(data = plot_frame_final, aes(color = final_cluster), size = 4)
  } else {
    g <- g + geom_path(data = plot_frame_full, color = "red", linejoin="round", size=1)
  }

  g + geom_point(data = plot_frame_init, color="black", size = 2) +
      geom_text_repel(data = plot_frame_init, mapping = aes(label = ObsLabel), size = 3) +
      guides(color = FALSE, size = FALSE) +
      theme(axis.title = element_text(size = 15),
            axis.text  = element_text(size = 10)) +
      xlab(axis[1]) +
      ylab(axis[2])
}
