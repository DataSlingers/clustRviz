#' Visualize the results of Convex BiClustering (\code{CBASS})
#'
#' \code{plot.CBASS} provides a range of ways to visualize the results of convex
#' clustering, including: \itemize{
#' \item Dendrograms, illustrating the nested cluster hierarchy inferred from
#'       the convex clustering solution path (\code{type = "obs.dendrogram"} and
#'       \code{type = "var.dendrogram"} for the observation and feature / variable
#'       clusterings, respectively);
#' \item Path plots, showing the coalescence of the estimated cluster centroids
#'       as the regularization parameter is increased (\code{type = "obs.dendrogram"}
#'       and \code{type = "var.dendrogram"} for the observation and feature / variable
#'       clusterings, respectively);
#' \item A cluster heatmap, displaying observation and feature histograms, as well
#'       as the clustered data matrix in a single graphic (\code{type = "heatmap"}); and
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
#' @param k.obs An integer indicating the desired number of observation clusters to be displayed
#'              in the static plots. (If plotting variables, the regularization level
#'              giving \code{k.obs} observation clusters is used.)
#' @param k.var An integer indicating the desired number of feature clusters to be displayed
#'              in the static plots. (If plotting variables, the regularization level
#'              giving \code{k.var} feature clusters is used.)
#' @param ... Additional arguments. Currently an error when \code{type} is not
#'            \code{"obs.dendrogram"} or \code{"var.dendrogram"}; passed to
#'            \code{\link[stats]{plot.dendrogram}} when \code{type == "obs.dendrogram"}
#'            or \code{type == "var.dendrogram"}.
#' @param dend.branch.width a positive number. Line width on dendrograms.
#' @param dend.labels.cex a positive number. Label size on dendrograms.
#' @param heatrow.label.cex heatmap row label size
#' @param heatcol.label.cex heatmap column label size
#' @return The value of the return type depends on the \code{type} argument:\itemize{
#'  \item if \code{type \%in\% c("obs.dendrogram", "var.dendrogram", "heatmap")},
#'         \code{x} is returned invisibly;
#'  \item if \code{type \%in\% c("obs.path", "var.path")}, an object of class
#'        \code{\link[ggplot2]{ggplot}} which can be plotted directly (by invoking
#'        its print method) or modified further by the user is returned;
#'  \item if \code{type = "interactive"}, a \code{\link[shiny]{shiny}} app which can be activated
#'        by invoking its print method.
#' } \code{saveviz.CBASS} always returns \code{file.name} invisibly.
#' @details The \code{\link{saveviz.CBASS}} function provides a unified interface
#'          for exporting \code{CBASS} visualizations to files. For all plots,
#'          at most one of \code{percent}, \code{k.obs}, and \code{k.var} must be supplied.
#' @importFrom shiny shinyApp fluidPage titlePanel tabsetPanel fluidRow animationOptions
#' @importFrom shiny tags column plotOutput renderPlot sliderInput
#' @importFrom stats as.dendrogram as.hclust quantile
#' @importFrom grDevices colorRampPalette adjustcolor
#' @export
#' @rdname plot_cbass
#' @examples
#' \dontrun{
#' cbass_fit <- CBASS(presidential_speech)
#' plot(cbass_fit, type='interactive')
#' }
plot.CBASS <- function(x,
                       ...,
                       type = c("heatmap",
                                "obs.dendrogram",
                                "var.dendrogram",
                                "obs.path",
                                "var.path",
                                "interactive"),
                       percent,
                       k.obs,
                       k.var,
                       dend.branch.width = 2,
                       dend.labels.cex = .6,
                       heatrow.label.cex = 1.5,
                       heatcol.label.cex = 1.5,
                       axis = c("PC1", "PC2")) {

  type <- match.arg(type)

  switch(
    type,
    obs.dendrogram = {
      cbass_dendro_plot(x,
                        percent = percent,
                        k.obs = k.obs,
                        k.var = k.var,
                        dend.branch.width = dend.branch.width,
                        dend.labels.cex = dend.labels.cex,
                        type = "obs",
                        ...)
    },
    var.dendrogram = {
      cbass_dendro_plot(x,
                        percent = percent,
                        k.obs = k.obs,
                        k.var = k.var,
                        dend.branch.width = dend.branch.width,
                        dend.labels.cex = dend.labels.cex,
                        type = "var",
                        ...)
    },
    obs.path = {
      cbass_path_plot(x,
                     axis = axis,
                     percent = percent,
                     k.obs = k.obs,
                     k.var = k.var,
                     ...,
                     type = "obs")
    },
    var.path = {
      cbass_path_plot(x,
                      axis = axis,
                      percent = percent,
                      k.obs = k.obs,
                      k.var = k.var,
                      ...,
                      type = "var")
    },
    heatmap = {
      cbass_heatmap_plot(x,
                         ...,
                         percent = percent,
                         k.obs = k.obs,
                         k.var = k.var,
                         heatrow.label.cex = heatrow.label.cex,
                         heatcol.label.cex = heatcol.label.cex)
    },
    interactive = {
      dots <- list(...)
      if ( length(dots) != 0 ){
        crv_error("Unknown arguments passed to ", sQuote("plot.CARP."))
      }

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

#' @noRd
#' @importFrom rlang .data
#' @importFrom dplyr filter select left_join pull
#' @importFrom ggplot2 ggplot geom_path aes geom_point guides theme element_text xlab ylab
#' @importFrom ggrepel geom_text_repel
cbass_path_plot <- function(x,
                            ...,
                            axis,
                            percent,
                            k.obs,
                            k.var,
                            type = c("obs", "var")){

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
  has_k.obs   <- !missing(k.obs)
  has_k.var   <- !missing(k.var)

  n_args <- has_percent + has_k.obs + has_k.var

  show_clusters <- (n_args == 1)

  if(n_args > 1){
    crv_error("At most one of ", sQuote("percent"), " ", sQuote("k.obs"), " and ",
              sQuote("k.var"), " may be supplied.")
  }

  if (n_args == 0L) {
    percent <- 1 ## By default, show the whole path
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

  path_obj <- if(type == "obs") x$cbass.cluster.path.vis.obs else x$cbass.cluster.path.vis.var

  if (any(plot_cols %not.in% colnames(path_obj))) {
    missing_col <- plot_cols[which(plot_cols %not.in% colnames(path_obj))][1]
    crv_error(sQuote(missing_col), " is not available for plotting.")
  }

  plot_frame_full <- path_obj %>% select(plot_cols) %>%
                                  filter(.data$Iter > x$burn.in)
  names(plot_frame_full)[1:2] <- c("V1", "V2")

  if(has_k.obs){

    if ( !is_integer_scalar(k.obs) ){
      crv_error(sQuote("k"), " must be an integer scalar (vector of length 1).")
    }

    if( k.obs <= 0 ) {
      crv_error(sQuote("k.obs"), " must be positive.")
    }

    if( k.obs > NROW(x$X) ){
      crv_error(sQuote("k.obs"), " cannot be more than the observations in the original data set (", NROW(x$X), ").")
    }

    percent <- x$cbass.cluster.path.vis.obs %>%
      select(.data$LambdaPercent, .data$NCluster) %>%
      filter(.data$NCluster <= k.obs) %>%
      select(.data$LambdaPercent) %>%
      summarize(percent = min(.data$LambdaPercent)) %>%
      pull
  }

  if(has_k.var){

    if ( !is_integer_scalar(k.var) ){
      crv_error(sQuote("k"), " must be an integer scalar (vector of length 1).")
    }

    if( k.var <= 0 ) {
      crv_error(sQuote("k.var"), " must be positive.")
    }

    if( k.var > NCOL(x$X) ){
      crv_error(sQuote("k.var"), " cannot be more than the features in the original data set (", NCOL(x$X), ").")
    }

    percent <- x$cbass.cluster.path.vis.var %>%
      select(.data$LambdaPercent, .data$NCluster) %>%
      filter(.data$NCluster <= k.var) %>%
      select(.data$LambdaPercent) %>%
      summarize(percent = min(.data$LambdaPercent)) %>%
      pull
  }

  if( !is_percent_scalar(percent) ){
    crv_error(sQuote("percent"), " must be a scalar between 0 and 1 (inclusive).")
  }

  plot_frame_full <- plot_frame_full %>% filter(.data$LambdaPercent <= percent)

  plot_frame_init  <- plot_frame_full %>% filter(.data$Iter == min(.data$Iter))
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
  } else {
    g <- g + geom_path(data = plot_frame_full, color = "red", linejoin="round", size=1)
  }

  g + geom_point(data = plot_frame_init, color="black", size = 2) +
    geom_text_repel(data = plot_frame_init, mapping = aes_string(label = "ObsLabel"), size = 3) +
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
                             k.obs,
                             k.var,
                             dend.branch.width = 2,
                             dend.labels.cex = .6,
                             type = c("obs", "var"),
                             ...){

  type <- match.arg(type)

  if(dend.branch.width <= 0){
    crv_error(sQuote("dend.branch.width"), " must be positive.")
  }

  if(dend.labels.cex <= 0){
    crv_error(sQuote("dend.labels.cex"), " must be positive.")
  }

  has_percent <- !missing(percent)
  has_k.obs   <- !missing(k.obs)
  has_k.var   <- !missing(k.var)
  n_args      <- has_percent + has_k.obs + has_k.var

  show_clusters <- (n_args == 1)

  if(n_args > 1){
    crv_error("At most one of ", sQuote("percent"), " ", sQuote("k.obs"), " and ",
         sQuote("k.var"), " may be supplied.")
  }

  dend <- if(type == "obs") x$cbass.dend.obs else x$cbass.dend.var

  ## Set better default margins
  par(mar = c(14, 7, 2, 1))

  dend %>% as.dendrogram %>%
           set("branches_lwd", dend.branch.width) %>%
           set("labels_cex", dend.labels.cex) %>%
           plot(ylab = "Amount of Regularization", cex.lab = 1.5, ...)

  if(show_clusters){
    labels <- get_cluster_labels(x, k.obs = k.obs, k.var = k.var, percent = percent, type = type)
    n_clusters <- nlevels(labels)

    my.cols <- adjustcolor(c("grey", "black"), alpha.f = .2)
    my.rect.hclust(dend, k = n_clusters, border = 2, my.col.vec = my.cols, lwd = 3)
  }

  invisible(x)
}

#' @importFrom grDevices colorRampPalette adjustcolor
#' @importFrom stats quantile as.dendrogram as.hclust
cbass_heatmap_plot <- function(x,
                               ...,
                               percent,
                               k.obs,
                               k.var,
                               heatrow.label.cex,
                               heatcol.label.cex,
                               breaks = NULL,
                               heatmap_col = NULL){

  dots <- list(...)

  if ( length(dots) != 0 ){
    crv_error("Unknown arguments passed to ", sQuote("plot.CBASS."))
  }

  has_percent <- !missing(percent)
  has_k.obs   <- !missing(k.obs)
  has_k.var   <- !missing(k.var)

  n_args <- has_percent + has_k.obs + has_k.var

  if(n_args >= 2){
    crv_error("At most one of ", sQuote("percent,"), " ", sQuote("k.obs"), " and ",
              sQuote("k.var"), " may be supplied.")
  }

  if(n_args == 0){
    U     <- get_clustered_data(x, percent = 0, refit = TRUE)
    k.col <- NCOL(U) ## FIXME - why is this backwards?
    k.row <- NROW(U)
  } else {
    U <- get_clustered_data(x, percent = percent, k.obs = k.obs, k.var = k.var, refit = TRUE)
    k.row <- nlevels(get_cluster_labels(x, percent = percent, k.obs = k.obs, k.var = k.var, type = "obs"))
    k.col <- nlevels(get_cluster_labels(x, percent = percent, k.obs = k.obs, k.var = k.var, type = "var"))
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
    breaks <- unique(quantile(as.vector(U), probs = quant.probs))
  }

  if (is.null(heatmap_col)) {
    nbreaks <- length(breaks)
    heatmap_col <- colorRampPalette(c("blue", "yellow"))(nbreaks - 1)
  }

  my.cols  <- adjustcolor(c("black", "grey"), alpha.f = .3)

  my.heatmap.2(
    x = U,
    scale = "none",
    Rowv = as.dendrogram(x$cbass.dend.obs),
    Colv = as.dendrogram(x$cbass.dend.var),
    trace = "none",
    density.info = "none",
    key = FALSE,
    breaks = breaks,
    col = heatmap_col,
    symkey = FALSE,
    Row.hclust = as.hclust(x$cbass.dend.obs),
    Col.hclust = as.hclust(x$cbass.dend.var),
    k.col = k.col,
    k.row = k.row,
    my.col.vec = my.cols,
    cexRow = heatrow.label.cex,
    cexCol = heatcol.label.cex,
    margins = c(14, 8)
  )

  invisible(x)
}
