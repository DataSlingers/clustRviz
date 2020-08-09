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
#' }
#'
#' @param x An object of class \code{CBASS} as returned by \code{\link{CBASS}}
#' @param type A string indicating the type of visualization to show (see details above).
#' @param dynamic A logical scalar.Should the resulting animation be dynamic (animated) or not?
#'                If \code{TRUE}, a dynamic visualization which varies along the CARP solution path at a
#'                grid given by \code{percent.seq} is produced. If \code{FALSE}, a fixed visualization at a single
#'                solution (determined by either \code{percent} or \code{k} if supplied) is produced.
#' @param interactive A logical scalar. Should the resulting animation be interactive or not?
#'                    If \code{TRUE}, an interactive visualization is produced by javascript(\code{\link{plotly}}).
#'                    If \code{FALSE}, a non-interactive visualization is produced by \code{\link[ggplot2]{ggplot}}.
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
#'            \item When \code{type} is \code{"heatmap"}, \code{"row.path"} or
#'                  \code{"col.path"}, the presence of
#'                  unknown arguments triggers an error;
#'            \item when \code{type == "row.dendrogram"} or \code{type == "col.dendrogram"},
#'                  \code{...} is forwarded to \code{\link[stats]{dendrogram}};
#'            \item when \code{type == "heatmap"}, \code{...} is forwarded to
#'                  \code{\link[heatmaply]{heatmaply}}.
#'            } See the documentation of the linked functions for details about
#'            additional supported arguments. \code{saveviz} passes arguments
#'            to the corresponding plot \code{type}.
#' @param slider_y A number to adjust the slider's vertical position for
#'                 interactive dendrogram plots (ignored for other plot types).
#' @param refit A logical scalar. Should "naive" centroids (TRUE) or the
#'              actual centroids estimated by convex clustering be used?
#'              When the default refit = FALSE, the estimated U from the convex
#'              clustering problem is used. The refit = TRUE returns actual
#'              centroids (mean) of all elements assigned to that cluster;
#'              Due to the global shrinkage imposed, these clusters are
#'              more "shrunk together" than the naive clusters. Only for the
#'              heatmap plots. (ignored for other plot types).
#' @return The value of the return type depends on the \code{type} argument:\itemize{
#'  \item if \code{type \%in\% c("row.dendrogram", "col.dendrogram", "heatmap")},
#'         \code{x} is returned invisibly;
#'  \item if \code{type \%in\% c("row.path", "col.path")}, an object of class
#'        \code{\link[ggplot2]{ggplot}} which can be plotted directly (by invoking
#'        its print method) or modified further by the user is returned;
#' } \code{saveviz.CBASS} always returns \code{file.name} invisibly.
#' @details The \code{\link{saveviz.CBASS}} function provides a unified interface
#'          for exporting \code{CBASS} visualizations to files. For all plots,
#'          at most one of \code{percent}, \code{k.row}, and \code{k.col} must be supplied.
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
                       slider_y = -0.3,
                       refit = FALSE) {

  type <- match.arg(type)

  if (!is_logical_scalar(dynamic)) {
    crv_error(sQuote("dynamic"), " must be a logical scalar.")
  }
  if (!is_logical_scalar(interactive)) {
    crv_error(sQuote("interactive"), " must be a logical scalar.")
  }

  switch(
    type,
    row.dendrogram = {
      if (!interactive){
        if (!dynamic){
          cbass_dendro_plot(x,
                            percent = percent,
                            k.row = k.row,
                            k.col = k.col,
                            type = "row",
                            ...)
        } else {
          cbass_dynamic_dendro_plot(x,
                                    percent.seq = percent.seq,
                                    type = "row",
                                    ...)
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
      if (!interactive){
        if (!dynamic){
          cbass_dendro_plot(x,
                            percent = percent,
                            k.row = k.row,
                            k.col = k.col,
                            type = "col",
                            ...)
        } else {
          cbass_dynamic_dendro_plot(x,
                                    percent.seq = percent.seq,
                                    type = "col",
                                    ...)
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
      if (!interactive){
        if (!dynamic){
          cbass_path_plot(x,
                          axis = axis,
                          percent = percent,
                          k.row = k.row,
                          k.col = k.col,
                          ...,
                          type = "row")
        } else {
          cbass_dynamic_path_plot(x,
                                  axis = axis,
                                  percent.seq = percent.seq,
                                  ...,
                                  type = "row")
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
      if (!interactive){
        if (!dynamic){
          cbass_path_plot(x,
                          axis = axis,
                          percent = percent,
                          k.row = k.row,
                          k.col = k.col,
                          ...,
                          type = "col")
        } else {
          cbass_dynamic_path_plot(x,
                                  axis = axis,
                                  percent.seq = percent.seq,
                                  ...,
                                  type = "col")
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
      if (!dynamic){
        if (!interactive){
          cbass_heatmap_plot(x,
                             ...,
                             percent = percent,
                             k.row = k.row,
                             k.col = k.col,
                             refit = refit)
        } else {
          cbass_heatmaply(x,
                          ...,
                          percent = percent,
                          k.row = k.row,
                          k.col = k.col)
        }
      } else {
        if (interactive){
          cbass_heatmaply_dynamic(x,
                                  ...,
                                  percent.seq = percent.seq,
                                  slider_y = slider_y)
        } else {
          cbass_dynamic_heatmap_plot(x,
                                     percent.seq = percent.seq,
                                     refit = refit,
                                     ...)
        }
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
      crv_error(sQuote("k.row"), " must be an integer scalar (vector of length 1).")
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
      crv_error(sQuote("k.col"), " must be an integer scalar (vector of length 1).")
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
#' @importFrom dplyr select filter rename left_join mutate bind_rows
#' @importFrom ggplot2 ggplot aes geom_path geom_point geom_text guides
#' @importFrom ggplot2 theme element_text xlab ylab
#' @importFrom gganimate transition_manual
cbass_dynamic_path_plot <- function(x, ..., axis, percent.seq, type = c("row", "col")){
  ## TODO - Combine this and cbass_path_plot as much as possible
  plot_frame_full <- get_feature_paths(x, axis, type = type) %>% rename(V1 = axis[1], V2 = axis[2])

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
#' @importFrom dplyr filter select summarize pull
#' @importFrom stats as.dendrogram
#' @importFrom dendextend set color_branches as.ggdend
#' @importFrom ggplot2 ggplot geom_segment scale_x_continuous scale_y_continuous theme
cbass_dendro_plot <- function(x,percent,
                              k.row,
                              k.col,
                              type = c("row", "col"),
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

  d <- as.dendrogram(x, type = type)

  if(show_clusters){
    if(has_percent){
      if (!is_percent_scalar(percent)) {
        crv_error(sQuote("percent"), " must be a scalar between 0 and 1 (inclusive).")
      }
    } else {
      if (!is_integer_scalar(k.row)|!is_integer_scalar(k.col)) {
        crv_error(sQuote("k.row"), "and", sQuote("k.col"), " must be an integer scalar (vector of length 1).")
      }
      if ( k.row <= 0 | k.col <= 0 ) {
        crv_error(sQuote("k.row"), "and", sQuote("k.col"), " must be positive.")
      }
      if ( k.row > NROW(x$X) | k.col > NCOL(x$X) ) {
        crv_error(sQuote("k.row"), "and", sQuote("k.col"), " cannot be more than the observations in the original data set (", NROW(x$X), ").")
      }

      percent <- get_feature_paths(x, features = character(), type = type) %>% filter(.data$NCluster == k) %>%
        select(.data$GammaPercent) %>%
        summarize(NCluster = mean(.data$GammaPercent)) %>%
        pull
    }

    labels <- get_cluster_labels(x, k.row = k.row, k.col = k.col, percent = percent, type = type)
    k <- nlevels(labels)

    c <- color_branches(d, k = k)
    dend <- as.ggdend(c)
    segs <- dend$segments

    adjust <- get_feature_paths(x, features = character(), type = type) %>% filter(.data$NCluster == 1) %>%
      select(.data$GammaPercent) %>%
      summarize(GammaPercent = min(.data$GammaPercent)) %>%
      pull
    segs$y <- segs$y*adjust
    segs$yend <- segs$yend*adjust

    tree <- as.hclust(x, type = type)
    cluster <- stats::cutree(tree, k = k)
    clustab <- table(cluster)[unique(cluster[tree$order])]
    clustsum <- cumsum(clustab)
    m <- c(0, clustsum) + 0.5
    line_x <- c(m[1],m[1],m)
    line_y <- c(0,percent,rep(percent,k+1))
    line_xend <- c(m[k+1],m[k+1],m)
    line_yend <- c(0,percent,rep(0,k+1))
    lines <- data.frame(x=line_x,y=line_y,xend=line_xend,yend=line_yend)

    p <- ggplot() +
      geom_segment(data = segs, aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend, color = .data$col), show.legend = FALSE) +
      geom_segment(data = lines, aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend, color = NA), show.legend = FALSE) +
      labs(title = paste0('Fraction of Regularization: ', percent * 100, '%\nNumber of Clusters: ', k))
  } else {
    dend <- as.ggdend(d)
    segs <- dend$segments
    segs$col <- "black"

    p <- ggplot() +
      geom_segment(data = segs, aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend), show.legend = FALSE)
  }

  label <- dend$label$label

  p +
    labs(y='Fraction of Regularization') +
    scale_x_continuous(breaks=1:length(label),labels= levels(label)) +
    scale_y_continuous(limits = c(0,1), breaks=c(0,0.25,0.5,0.75,1),labels = c("0%", "25%", "50%", "75%", "100%"))+
    theme(axis.text.x=element_text(hjust=1,vjust=0.5,angle=90),
          axis.title.x=element_blank(),
          panel.background = element_blank(),
          panel.grid=element_blank(),
          panel.border=element_blank())
}

#' @noRd
#' @importFrom rlang .data
#' @importFrom dplyr filter select summarize pull
#' @importFrom stats as.dendrogram
#' @importFrom dendextend set color_branches as.ggdend
#' @importFrom ggplot2 ggplot geom_segment scale_x_continuous scale_y_continuous theme labs
#' @importFrom gganimate transition_time
cbass_dynamic_dendro_plot <- function(x,
                                     percent.seq,
                                     type = c("row", "col"),
                                     ...){
  adjust <- get_feature_paths(x, features = character(), type = type) %>% filter(.data$NCluster == 1) %>%
    select(.data$GammaPercent) %>%
    summarize(GammaPercent = min(.data$GammaPercent)) %>%
    pull

  segs_dynamic <- data.frame()
  lines_dynamic <- data.frame()
  for (per in percent.seq){
    labels <- get_cluster_labels(x, percent = per, type = type)
    k <- nlevels(labels)

    # segments
    d <- as.dendrogram(x, type = type)
    c <- color_branches(d, k = k)
    dend <- as.ggdend(c)
    segs <- dend$segments
    label <- dend$label$label

    segs$y <- segs$y*adjust
    segs$yend <- segs$yend*adjust

    segs_dynamic <- rbind(segs_dynamic,cbind(segs, reg = per*100))

    # lines
    tree <- as.hclust(x, type = type)
    cluster <- stats::cutree(tree, k = k)
    clustab <- table(cluster)[unique(cluster[tree$order])]
    clustsum <- cumsum(clustab)
    m <- c(0, clustsum) + 0.5
    line_x <- c(m[1],m[1],m)
    line_y <- c(0,per,rep(per,k+1))
    line_xend <- c(m[k+1],m[k+1],m)
    line_yend <- c(0,per,rep(0,k+1))
    lines <- data.frame(reg=per*100,x=line_x,y=line_y,xend=line_xend,yend=line_yend)
    lines_dynamic <- rbind(lines_dynamic,lines)
  }

  ggplot() +
    geom_segment(data = segs_dynamic, aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend, color = .data$col), show.legend = FALSE) +
    geom_segment(data = lines_dynamic, aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend, color = NA), show.legend = FALSE) +
    labs(y='Fraction of Regularization') +
    scale_x_continuous(breaks=1:length(label),labels= levels(label)) +
    scale_y_continuous(limits = c(0,1), breaks=c(0,0.25,0.5,0.75,1),labels = c("0%", "25%", "50%", "75%", "100%")) +
    theme(axis.text.x=element_text(hjust=1,vjust=0.5,angle=90),
          axis.title.x=element_blank(),
          panel.background = element_blank(),
          panel.grid=element_blank(),
          panel.border=element_blank()) +
    transition_manual(.data$reg) +
    labs(title = paste0('Fraction of Regularization: ', '{current_frame}', '%'))
}

#' @noRd
#' @importFrom rlang .data
#' @importFrom dplyr filter select summarize pull
#' @importFrom stats as.dendrogram
#' @importFrom dendextend set color_branches as.ggdend
#' @importFrom ggplot2 ggplot geom_segment scale_x_continuous scale_y_continuous theme geom_tile scale_fill_gradient2 labs
cbass_heatmap_plot <- function(x,
                              percent,
                              k.row,
                              k.col,
                              show_clusters = (n_args == 1L),
                              refit = FALSE,
                              ...){
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

  d_row <- as.dendrogram(x, type = "row")
  d_col <- as.dendrogram(x, type = "col")

  adjust_row <- get_feature_paths(x, features = character(), "row") %>% filter(.data$NCluster == 1) %>%
    select(.data$GammaPercent) %>%
    summarize(GammaPercent = min(.data$GammaPercent)) %>%
    pull

  adjust_col <- get_feature_paths(x, features = character(), "col") %>% filter(.data$NCluster == 1) %>%
    select(.data$GammaPercent) %>%
    summarize(GammaPercent = min(.data$GammaPercent)) %>%
    pull

  if(show_clusters){
    k.row <- nlevels(get_cluster_labels(x, k.row = k.row, k.col = k.col, percent = percent, type = "row"))
    if(!has_percent){
      percent <- get_feature_paths(x, features = character(), type = "row") %>% filter(.data$NCluster == k.row) %>%
        select(.data$GammaPercent) %>%
        summarize(NCluster = mean(.data$GammaPercent)) %>%
        pull
    }
    k.col <- nlevels(get_cluster_labels(x, percent = percent, type = "col"))

    # heatmap
    U <- get_clustered_data(x, percent = percent, refit = refit)

    # dendrogram
    segs_row <- adjusted_dendrogram(d_row, rev = TRUE, k = k.row, cluster = TRUE, adjust = adjust_row)
    segs_col <- adjusted_dendrogram(d_col, rev = FALSE, k = k.col, cluster = TRUE, adjust = adjust_col)

    lines_row <- dendrogram_box(x, rev = TRUE, k = k.row, type = "row", percent = percent)
    lines_col <- dendrogram_box(x, rev = FALSE, k = k.col, type = "col", percent = percent)

    p <- ggplot() +
      geom_segment(data = segs_row, aes(x = (.data$y/3+1)*length(cn)+0.5, xend = (.data$yend/3+1)*length(cn)+0.5, y = .data$x, yend = .data$xend, color = .data$col), show.legend = FALSE) +
      geom_segment(data = segs_col, aes(y = (.data$y/3+1)*length(rn)+0.5, yend = (.data$yend/3+1)*length(rn)+0.5, x = .data$x, xend = .data$xend, color = .data$col), show.legend = FALSE) +
      geom_segment(data = lines_row, aes(x = (.data$y/3+1)*length(cn)+0.5, xend = (.data$yend/3+1)*length(cn)+0.5, y = .data$x, yend = .data$xend, color = NA), show.legend = FALSE) +
      geom_segment(data = lines_col, aes(y = (.data$y/3+1)*length(rn)+0.5, yend = (.data$yend/3+1)*length(rn)+0.5, x = .data$x, xend = .data$xend, color = NA), show.legend = FALSE) +
      labs(title = paste0('Fraction of Regularization: ', percent * 100, '%'))
  } else {
    # heatmap
    U     <- get_clustered_data(x, percent = 0, refit = refit)
    k.row <- NROW(U)
    k.col <- NCOL(U)

    # dendrogram
    segs_row <- adjusted_dendrogram(d_row, rev = TRUE, k = k.row, cluster = FALSE, adjust = adjust_row)
    segs_col <- adjusted_dendrogram(d_col, rev = FALSE, k = k.col, cluster = FALSE, adjust = adjust_col)

    p <- ggplot() +
      geom_segment(data = segs_row, aes(x = (.data$y/3+1)*length(cn)+0.5, xend = (.data$yend/3+1)*length(cn)+0.5, y = .data$x, yend = .data$xend), show.legend = FALSE) +
      geom_segment(data = segs_col, aes(y = (.data$y/3+1)*length(rn)+0.5, yend = (.data$yend/3+1)*length(rn)+0.5, x = .data$x, xend = .data$xend), show.legend = FALSE)
  }

  # heatmap
  h <- heatmapr(
    U,
    Rowv = as.hclust(x),
    dendrogram = "row",
    k_row = k.row,
    k_col = k.col
  )
  data_mat <- h$matrix$data
  dat <- data.frame(value = as.vector(data_mat),
                    x = as.vector(col(data_mat)),
                    y = as.vector(row(data_mat)))

  rn <- rownames(h$matrix$data)
  cn <- colnames(h$matrix$data)
  r <- range(dat$value)

  p +
    geom_tile(data =  dat, aes(x = .data$x, y = .data$y,fill = .data$value)) +
    scale_fill_gradient2(low = "#313695", mid = "#FFFFBF", high = "#A50026", midpoint = (floor(r[1])+ceiling(r[2]))/2, limits = c(floor(r[1]),ceiling(r[2]))) +
    scale_x_continuous(breaks=1:length(cn),labels= cn) +
    scale_y_continuous(breaks=1:length(rn),labels= rn) +
    theme(axis.text.x=element_text(hjust=1,vjust=0.5,angle=90),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background = element_blank(),
          panel.grid=element_blank(),
          panel.border=element_blank())
}

### utils for the heatmap

adjusted_dendrogram <- function(d, rev = FALSE, k, cluster = FALSE, adjust){
  if (cluster == TRUE){
    c <- color_branches(d, k = k)
  } else {
    c <- d
  }
  if (rev == TRUE) {
    c <- rev(c)
  }
  dend <- as.ggdend(c)
  segs <- dend$segments

  if (cluster == FALSE){
    segs$col <- "black"
  }

  segs$y <- segs$y*adjust
  segs$yend <- segs$yend*adjust

  return(segs)
}

dendrogram_box <- function(x, rev = FALSE, k, type, percent){
  tree <- as.hclust(x, type = type)
  if (rev == TRUE) {
    tree <- rev(tree)
  }
  cluster <- stats::cutree(tree, k = k)
  clustab <- table(cluster)[unique(cluster[tree$order])]
  clustsum <- cumsum(clustab)
  m <- c(0, clustsum) + 0.5
  line_x <- c(m[1],m[1],m)
  line_y <- c(0,percent,rep(percent,k+1))
  line_xend <- c(m[k+1],m[k+1],m)
  line_yend <- c(0,percent,rep(0,k+1))
  lines <- data.frame(x=line_x,y=line_y,xend=line_xend,yend=line_yend)

  return(lines)
}

#' @noRd
#' @importFrom rlang .data
#' @importFrom dplyr filter select summarize pull
#' @importFrom stats as.dendrogram
#' @importFrom dendextend set color_branches as.ggdend
#' @importFrom ggplot2 ggplot geom_segment scale_x_continuous scale_y_continuous theme geom_tile scale_fill_gradient2 labs
#' @importFrom gganimate transition_manual
cbass_dynamic_heatmap_plot <- function(x,
                                      percent.seq,
                                      refit = FALSE,
                                      ...){
  adjust_row <- get_feature_paths(x, features = character(), "row") %>% filter(.data$NCluster == 1) %>%
    select(.data$GammaPercent) %>%
    summarize(GammaPercent = min(.data$GammaPercent)) %>%
    pull

  adjust_col <- get_feature_paths(x, features = character(), "col") %>% filter(.data$NCluster == 1) %>%
    select(.data$GammaPercent) %>%
    summarize(GammaPercent = min(.data$GammaPercent)) %>%
    pull

  data_mat_dynamic <- data.frame()
  segs_row_dynamic <- data.frame()
  segs_col_dynamic <- data.frame()
  lines_row_dynamic <- data.frame()
  lines_col_dynamic <- data.frame()
  for (per in percent.seq){
    U <- get_clustered_data(x, percent = per, refit = refit)
    k.row <- nlevels(get_cluster_labels(x, percent = per, type = "row"))
    k.col <- nlevels(get_cluster_labels(x, percent = per, type = "col"))
    h <- heatmapr(
      U,
      Rowv = as.hclust(x),
      dendrogram = "row",
      k_row = k.row,
      k_col = k.col
    )
    # data for heatmap
    data_mat <- h$matrix$data
    data_mat_dynamic <- rbind(data_mat_dynamic,cbind(value = as.vector(data_mat),
                                                     x = as.vector(col(data_mat)),
                                                     y = as.vector(row(data_mat)),
                                                     reg = as.vector(per)*100))
    rn <- rownames(h$matrix$data)
    cn <- colnames(h$matrix$data)

    # segments
    d_row <- as.dendrogram(x, type = "row")
    d_col <- as.dendrogram(x, type = "col")
    segs_row <- adjusted_dendrogram(d_row, rev = TRUE, k = k.row, cluster = TRUE, adjust = adjust_row)
    segs_col <- adjusted_dendrogram(d_col, rev = FALSE, k = k.col, cluster = TRUE, adjust = adjust_col)
    lines_row <- dendrogram_box(x, rev = TRUE, k = k.row, type = "row", percent = per)
    lines_col <- dendrogram_box(x, rev = FALSE, k = k.col, type = "col", percent = per)

    segs_row_dynamic <- rbind(segs_row_dynamic,cbind(segs_row, reg = per*100))
    segs_col_dynamic <- rbind(segs_col_dynamic,cbind(segs_col, reg = per*100))
    lines_row_dynamic <- rbind(lines_row_dynamic,cbind(lines_row, reg = per*100))
    lines_col_dynamic <- rbind(lines_col_dynamic,cbind(lines_col, reg = per*100))
  }

  r <- range(data_mat_dynamic$value)

  ggplot(data =  data_mat_dynamic, aes(x = .data$x, y = .data$y)) +
    geom_tile(aes(fill = .data$value)) +
    scale_fill_gradient2(low = "#313695", mid = "#FFFFBF", high = "#A50026", midpoint = (floor(r[1])+ceiling(r[2]))/2, limits = c(floor(r[1]),ceiling(r[2]))) +
    scale_x_continuous(breaks=1:length(cn),labels= cn) +
    scale_y_continuous(breaks=1:length(rn),labels= rn) +
    geom_segment(data = segs_row_dynamic, aes(x = (.data$y/3+1)*length(cn)+0.5, xend = (.data$yend/3+1)*length(cn)+0.5, y = .data$x, yend = .data$xend, color = .data$col), show.legend = FALSE) +
    geom_segment(data = segs_col_dynamic, aes(y = (.data$y/3+1)*length(rn)+0.5, yend = (.data$yend/3+1)*length(rn)+0.5, x = .data$x, xend = .data$xend, color = .data$col), show.legend = FALSE) +
    geom_segment(data = lines_row_dynamic, aes(x = (.data$y/3+1)*length(cn)+0.5, xend = (.data$yend/3+1)*length(cn)+0.5, y = .data$x, yend = .data$xend, color = NA), show.legend = FALSE) +
    geom_segment(data = lines_col_dynamic, aes(y = (.data$y/3+1)*length(rn)+0.5, yend = (.data$yend/3+1)*length(rn)+0.5, x = .data$x, xend = .data$xend, color = NA), show.legend = FALSE) +
    theme(axis.text.x=element_text(hjust=1,vjust=0.5,angle=90),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background = element_blank(),
          panel.grid=element_blank(),
          panel.border=element_blank()) +
    transition_manual(.data$reg) +
    labs(title = paste0('Fraction of Regularization: ', '{current_frame}', '%'))
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

  if(n_args >= 1 & dynamic){
    crv_error("Can't set ", sQuote("percent"), " ", sQuote("k.row"), " and ", sQuote("k.col"), " for dynamic plots")
  }

  if (n_args == 0L) {
    percent <- 1 ## By default, show the whole path
  }

  plot_frame_full <- get_feature_paths(x, axis, type = type) %>% rename(V1 = axis[1], V2 = axis[2])

  if(has_k.row){

    if ( !is_integer_scalar(k.row) ){
      crv_error(sQuote("k.row"), " must be an integer scalar (vector of length 1).")
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
      crv_error(sQuote("k.col"), " must be an integer scalar (vector of length 1).")
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

  if (!dynamic){
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

  if (!dynamic){
    if(show_clusters){
      if(has_percent){
        if (!is_percent_scalar(percent)) {
          crv_error(sQuote("percent"), " must be a scalar between 0 and 1 (inclusive).")
        }

        # if(percent == 0){ ## Don't bother showing boxes if percent is 0 and bail early
        #   return(invisible(x))
        # }
      } else {
        if (!is_integer_scalar(k.row)|!is_integer_scalar(k.col)) {
          crv_error(sQuote("k.row"), "and", sQuote("k.col"), " must be an integer scalar (vector of length 1).")
        }
        if ( k.row <= 0 | k.col <= 0 ) {
          crv_error(sQuote("k.row"), "and", sQuote("k.col"), " must be positive.")
        }
        if ( k.row > NROW(x$X) | k.col > NCOL(x$X) ) {
          crv_error(sQuote("k.row"), "and", sQuote("k.col"), " cannot be more than the observations in the original data set (", NROW(x$X), ").")
        }

        percent <- get_feature_paths(x, features = character(), type = type) %>% filter(.data$NCluster == k) %>%
          select(.data$GammaPercent) %>%
          summarize(NCluster = mean(.data$GammaPercent)) %>%
          pull
      }

      labels <- get_cluster_labels(x, k.row = k.row, k.col = k.col, percent = percent, type = type)
      k <- nlevels(labels)

      cluster <- labels
      clustab <- table(cluster)[unique(cluster[labels(d)])]
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

    c <- dynamic_seg %>% group_by(.data$x, .data$xend) %>% summarize(count = n())
    dynamic_sec <- arrange(merge(dynamic_seg, c), .data$Regularization, desc(.data$count), .data$x, .data$y)

    # avoid the situation when very short lines cannot be shown in the plotly
    tidy_segments_dynamic <- tidy_segments[,1:4] %>%
      mutate(ylength = abs(.data$y-.data$yend))
    tidy_segments_dynamic$ylength[tidy_segments_dynamic$ylength==0] <- 1

    dendro_dynamic <- plot_ly(
      hoverinfo = "none") %>%
      add_segments(
        data = tidy_segments_dynamic %>% arrange(desc(.data$ylength)),
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

#' @noRd
#' @importFrom heatmaply heatmapr
#' @importFrom dendextend as.ggdend
#' @importFrom plotly plot_ly add_heatmap add_segments animation_slider
#' @importFrom tibble tibble
cbass_heatmaply_dynamic <- function(x,
                                   ...,
                                   percent.seq = percent.seq,
                                   slider_y = slider_y){
  data_mat_dynamic <- data.frame()
  seg_dynamic_row <- data.frame()
  seg_dynamic_col <- data.frame()
  for (per in percent.seq){
    U <- get_clustered_data(x, percent = per, refit = TRUE)
    k.row <- nlevels(get_cluster_labels(x, percent = per, type = "row"))
    k.col <- nlevels(get_cluster_labels(x, percent = per, type = "col"))
    h <- heatmapr(
      U,
      Rowv = as.hclust(x, type = "row"),
      Colv = as.hclust(x, type = "col"),
      dendrogram = "both",
      k_row = k.row,
      k_col = k.col
    )
    # data for heatmap
    data_mat <- h$matrix$data
    data_mat_dynamic <- rbind(data_mat_dynamic,cbind(value = as.vector(data_mat),
                                                     x = as.vector(col(data_mat)),
                                                     y = as.vector(row(data_mat)),
                                                     percent = as.vector(per)))

    # data for dendrogram
    dend_row <- h$rows
    dend_col <- h$cols
    segs_row <- as.ggdend(dend_row)$segments
    segs_col <- as.ggdend(dend_col)$segments
    segs_row$col[is.na(segs_row$col)] <- "black" # default value for NA is "black"
    segs_col$col[is.na(segs_col$col)] <- "black" # default value for NA is "black"
    seg_dynamic_row <- rbind(seg_dynamic_row,cbind(seg = seq_along(segs_row$x), segs_row, percent = per))
    seg_dynamic_col <- rbind(seg_dynamic_col,cbind(seg = seq_along(segs_col$x), segs_col, percent = per))

    rn <- rownames(h$matrix$data)
    cn <- colnames(h$matrix$data)
  }

  colors_row <- sort(unique(seg_dynamic_row$col))
  colors_col <- sort(unique(seg_dynamic_col$col))

  p <- plot_ly(colorbar = list(title = "")) %>%
    add_heatmap(data = data_mat_dynamic,
                z = ~value, x = ~x, y = ~y,
                text = ~paste0("column: ", cn[x], "<br>", "row: ", rn[y], "<br>", "value: ", value),
                showlegend = FALSE,
                hoverinfo = "text",
                frame = ~percent * 100) %>%
    plotly::layout(
      xaxis = list(
        title = "",
        tickvals = seq_along(cn), ticktext = cn,
        range = c(0.5, 4/3 * length(cn) + 0.5),
        showticklabels = TRUE,
        showgrid = FALSE),
      yaxis = list(
        title = "",
        tickvals = seq_along(rn), ticktext = rn,
        range = c(0.5, 4/3 * length(rn) + 0.5),
        showticklabels = TRUE,
        showgrid = FALSE))

  for (i in seq_along(levels(factor(seg_dynamic_row$seg)))){
    p <- p %>%
      add_segments(data = seg_dynamic_row[seg_dynamic_row$seg==i,],
                   x = ~(y/3+1)*length(cn)+0.5, xend = ~(yend/3+1)*length(cn)+0.5,
                   y = ~x, yend = ~xend,
                   color = ~col,
                   showlegend = FALSE,
                   colors = colors_row,
                   hoverinfo = "none",
                   frame = ~percent*100)
  }
  for (i in seq_along(levels(factor(seg_dynamic_col$seg)))){
    p <- p %>%
      add_segments(data = seg_dynamic_col[seg_dynamic_col$seg==i,],
                   x = ~x, xend = ~xend,
                   y = ~(y/3+1)*length(rn)+0.5, yend = ~(yend/3+1)*length(rn)+0.5,
                   color = ~col,
                   showlegend = FALSE,
                   colors = colors_col,
                   hoverinfo = "none",
                   frame = ~percent*100)
  }

  p %>%
    animation_slider(y = slider_y,
                     currentvalue = list(prefix = "Regularization: ", suffix = "%"))
}
