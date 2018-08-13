#' Plot function for CARP objects
#'
#' \code{plot.CARP} displays both static and interactive visualizations of
#' the CARP clustering solution path, including dendrograms and clustering
#' paths
#'
#' Possible visualizations of the CARP solution path are: (i) a static
#' dendrogram representing the clustering solution path at various
#' levels of regularization; (ii) a static clustering path representing
#' the coalescence of observations projected via PCA; and (iii) a
#' interactive dendorgram and clustering path visualization, showing
#' how the clustering solutions change with increased reguarliazation.
#' The first interactive tab shows a movie of real-time cluster solutions
#' while the second tabs shows the cluster solution for various numbers
#' of clusters
#' @param x a CARP object returned by \code{CARP}
#' @param type a string specifying the type of plot to produce. 'dendrogram'
#' produces the static cluster dendrogram; 'path' produces the static
#' cluster path; and 'interactive' produces an interactive visualization of
#' both the cluster path and dendrogram
#' @param axis a character vector of length two with elements as 'PC1','PC2',..
#' etc. Specifics which principal component axis to display for the 'path'
#' visualization
#' @param percent a number between 0 and 1. Specifies how far along the
#' CARP path the 'path' visualization should display.
#' @param max.nclust a positive integer. The maximum number of clusters
#' to display in the interactive plot.
#' @param min.nclust a positive value. The minimum number of clusters to
#' display in the interactive plot.
#' @param ... Unused additional generic arguements
#' @param dend.branch.width a positive number. Line width on dendrograms.
#' @param dend.labels.cex a positive number. Label size on dendrograms.
#' @importFrom shiny shinyApp
#' @importFrom shiny fluidPage
#' @importFrom shiny titlePanel
#' @importFrom shiny tabsetPanel
#' @importFrom shiny fluidRow
#' @importFrom shiny animationOptions
#' @importFrom shiny column
#' @importFrom shiny plotOutput
#' @importFrom shiny sliderInput
#' @importFrom shiny uiOutput
#' @importFrom shiny renderUI
#' @importFrom shiny tags
#' @importFrom shiny checkboxGroupInput
#' @importFrom shiny renderPlot
#' @importFrom stats as.dendrogram
#' @importFrom stats median
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
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr distinct
#' @importFrom dplyr rename
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom dplyr select_
#' @importFrom dplyr %>%
#' @importFrom grDevices adjustcolor
#' @importFrom RColorBrewer brewer.pal
#' @export
#' @examples
#' \dontrun{
#' library(clustRviz)
#' data("presidential_speech")
#' Xdat <- presidential_speech[1:10,1:4]
#' carp.fit <- CARP(X=Xdat)
#' plot(carp.fit,type='interactive')
#' }
plot.CARP <- function(
                      x,
                      type = c("dendrogram", "path", "interactive"),
                      axis = c("PC1", "PC2"),
                      dend.branch.width = 2,
                      dend.labels.cex = .6,
                      percent = 1,
                      max.nclust = 9,
                      min.nclust = 1,
                      ...) {
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
          ggplot2::aes(x = V1, y = V2),
          data = plot.frame %>% dplyr::filter(Iter == 1),
          color = "black",
          size = I(2)
        ) +
        ggrepel::geom_text_repel(
          ggplot2::aes(x = V1, y = V2, label = ObsLabel),
          size = I(3),
          data = plot.frame %>% dplyr::filter(Iter == 1)
        ) +
        ggplot2::guides(color = FALSE, size = FALSE) +
        ggplot2::theme(axis.title = ggplot2::element_text(size = 15)) +
        ggplot2::theme(axis.text = ggplot2::element_text(size = 10)) +
        ggplot2::xlab(axis[1]) +
        ggplot2::ylab(axis[2])
    },
    interactive = {
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
