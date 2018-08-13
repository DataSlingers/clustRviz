#' Plot function for CBASS objects
#'
#' \code{plot.CBASS} displays both static and interactive visualizations of
#' the CBASS clustering solution path, including dendrograms and heatmaps.
#'
#' Possible visualizations of the CBASS solution path are: (i) a static
#' dendrogram representing the clustering solution path at various
#' levels of regularization (both observations and variables); and (ii) an
#' interactive cluster heatmap, showing
#' how the clustering solutions change with increased reguarliazation.
#' @param x a CBASS object returned by \code{CARP}
#' @param type a string specifying the type of plot to produce.
#' 'obs.dendrogram' produces the static cluster dendrogram for the observations;
#' 'var.dendrogram' produces the static cluster dendrogram for the variables;
#' 'heatmap' produces the static cluster heatmap with assiociated observation
#' and variable dendrograms;
#' and 'interactive' produces an interactive visualization of cluster
#' heatmap and its associated dendrograms.
#' @param ... Unused additional generic arguements
#' @param dend.branch.width a positive number. Line width on dendrograms.
#' @param dend.labels.cex a positive number. Label size on dendrograms.
#' @param heatrow.label.cex heatmap row label size
#' @param heatcol.label.cex heatmap column label size
#' @importFrom shiny shinyApp
#' @importFrom shiny fluidPage
#' @importFrom shiny titlePanel
#' @importFrom shiny tabsetPanel
#' @importFrom shiny fluidRow
#' @importFrom shiny animationOptions
#' @importFrom shiny tags
#' @importFrom shiny column
#' @importFrom shiny plotOutput
#' @importFrom shiny renderPlot
#' @importFrom shiny sliderInput
#' @importFrom stats as.dendrogram
#' @importFrom stats as.hclust
#' @importFrom stats quantile
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices adjustcolor
#' @export
#' @examples
#' \dontrun{
#' library(clustRviz)
#' data("presidential_speech")
#' Xdat <- presidential_speech[1:10,1:4]
#' cbass.fit <- CARP(X=Xdat)
#' plot(cbass.fit,type='interactive')
#' }
plot.CBASS <- function(
                       x,
                       type = c("obs.dendrogram", "var.dendrogram", "heatmap", "interactive"),
                       dend.branch.width = 2,
                       dend.labels.cex = .6,
                       heatrow.label.cex = 1.5,
                       heatcol.label.cex = 1.5,
                       ...) {
  type <- match.arg(type)
  switch(
    type,
    obs.dendrogram = {
      x$cbass.dend.obs %>%
        stats::as.dendrogram() %>%
        dendextend::set("branches_lwd", dend.branch.width) %>%
        dendextend::set("labels_cex", dend.labels.cex) %>%
        plot(ylab = "Amount of Regularization")
    },
    var.dendrogram = {
      x$cbass.dend.var %>%
        stats::as.dendrogram() %>%
        dendextend::set("branches_lwd", dend.branch.width) %>%
        dendextend::set("labels_cex", dend.labels.cex) %>%
        plot(ylab = "Amount of Regularization")
    },
    heatmap = {
      if (x$X.center.global) {
        X <- x$X
        X <- X - mean(X)
        X <- t(X)
      } else {
        X <- t(x$X.orig)
      }
      rownames(X) <- x$var.labels
      colnames(X) <- x$obs.labels
      nbreaks <- 50
      quant.probs <- seq(0, 1, length.out = nbreaks)
      breaks <- unique(stats::quantile(X[TRUE], probs = quant.probs))
      nbreaks <- length(breaks)
      heatcols <- grDevices::colorRampPalette(c("blue", "yellow"))(nbreaks - 1)

      my.cols <- grDevices::adjustcolor(c("black", "grey"), alpha.f = .3)
      my.heatmap.2(
        x = X,
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
        k.col = x$n.obs,
        k.row = x$p.var,
        my.col.vec = my.cols,
        cexRow = heatrow.label.cex,
        cexCol = heatcol.label.cex,
        margins = c(14, 8)
      )
    },
    interactive = {
      shiny::shinyApp(
        ui = shiny::fluidPage(
          shiny::tags$style(
            type = "text/css",
            ".recalculating { opacity: 1.0; }"
          ),


          shiny::titlePanel("BiClustering"),
          shiny::tabsetPanel(
            shiny::tabPanel(
              "Heatmap",
              shiny::fluidRow(
                shiny::column(
                  3,
                  shiny::sliderInput(
                    "regcent",
                    "Amount of Regularization",
                    min = 0,
                    max = 1,
                    value = .5,
                    step = .03, animate = shiny::animationOptions(interval = 300, loop = T)
                  )
                ),
                shiny::column(
                  9,
                  shiny::plotOutput("heatmap", height = "900px", width = "1200px")
                )
              )
            )
          )
        ),
        server = function(input, output) {
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
          output$heatmap <- shiny::renderPlot({
            plt.iter <- which.min(abs(input$regcent - lam.prop.seq))
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
          })
        }
      )
    }
  )
}
