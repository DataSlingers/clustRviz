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
#'       (\code{dynamic = TRUE})
#' }
#'
#' @param x An object of class \code{CARP} as returned by \code{\link{CARP}}
#' @param type A string indicating the type of visualization to show (see details above).
#' @param dynamic A logical scalar.Should the resulting animation be dynamic (animated) or not?
#'                If \code{TRUE}, a dynamic visualization which varies along the CARP solution path at a
#'                grid given by \code{percent.seq} is produced. If \code{FALSE}, a fixed visualization at a single
#'                solution (determined by either \code{percent} or \code{k} if supplied) is produced.
#' @param interactive A logical scalar. Should the resulting animation be interactive or not?
#'                    If \code{TRUE}, an interactive visualization is produced by javascript(\code{\link{plotly}}).
#'                    If \code{FALSE}, a non-interactive visualization is produced by \code{\link[ggplot2]{ggplot}}.
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
#' @param percent.seq A grid of values of percent along which to generate dynamic
#'                    visualizations (if dynamic == TRUE)
#' @param ... Additional arguments, which are handled differently for different
#'            values of \code{type}.\itemize{
#'            \item When \code{type} is \code{"path"}, the presence of
#'                  unknown arguments triggers an error;
#'            \item when \code{type == "dendrogram"}
#'                  \code{...} is forwarded to \code{\link[stats]{dendrogram}}; and
#'            \item when \code{type == "heatmap"} and \code{interactive == TRUE},
#'                  \code{...} is forwarded to \code{\link[heatmaply]{heatmaply}}.
#'            } See the documentation of \code{\link[ggplot2]{ggplot}} (for \code{interactive == FALSE})
#'            or \code{\link[plotly]{plotly}} (for \code{interactive == TRUE}) for details about additional
#'            supported arguments to the corresponding plot \code{type}.
#' @param slider_y A number to adjust the slider's vertical position for
#'                 interactive dendrogram and interactive heatmap plots
#'                 (ignored for other plot types).
#' @param refit A logical scalar. Should "naive" centroids (TRUE) or the
#'              actual centroids estimated by convex clustering be used?
#'              When the default refit = FALSE, the estimated U from the convex
#'              clustering problem is used. The refit = TRUE returns actual
#'              centroids (mean) of all elements assigned to that cluster;
#'              Due to the global shrinkage imposed, these clusters are
#'              more "shrunk together" than the naive clusters. Only for the
#'              heatmap plots. (ignored for other plot types).
#' @return The value of the return type depends on the \code{interactive} and \code{dynamic} arguments:\itemize{
#'   \item if \code{interactive = FALSE} and \code{dynamic = FALSE}, an object of class \code{\link[ggplot2]{ggplot}}
#'         is returned;
#'   \item if \code{interactive = FALSE} and \code{dynamic = TRUE}, an object of class \code{\link[gganimate:gganimate-package]{gganim}}
#'         is returned;
#'   \item if \code{interactive = TRUE}, an object of class \code{\link[plotly]{plotly}}
#'         is returned.
#' } All the plots can be plotted directly (by invoking its print method) or further
#'   manipulated by the user.
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
#' carp_fit <- CARP(presidential_speech)
#' plot(carp_fit, type = "path")
#' plot(carp_fit, type = "dendrogram")
#' plot(carp_fit, type = "heatmap")
#' plot(carp_fit, type = "heatmap", dynamic = TRUE)
#' \dontrun{
#' plot(carp_fit, type='heatmap', interactive=TRUE)
#' }
plot.CARP <- function(x,
                      ...,
                      type = c("dendrogram", "path", "heatmap"),
                      dynamic = FALSE,
                      interactive = FALSE,
                      axis = c("PC1", "PC2"),
                      percent,
                      k,
                      percent.seq = seq(0, 1, 0.01),
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
    dendrogram = {
      if (!interactive) {
        if (!dynamic) {
          carp_dendro_plot(x,
                           percent = percent,
                           k = k,
                           ...)
        } else {
          carp_dynamic_dendro_plot(x,
                                   percent.seq = percent.seq,
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
      if (!interactive) {
        if (!dynamic) {
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
    heatmap = {
      if (interactive){
        if (!dynamic){
          carp_heatmaply(x,
                         ...,
                         percent = percent,
                         k = k)
        } else {
          carp_heatmaply_dynamic(x,
                                 ...,
                                 percent.seq = percent.seq,
                                 slider_y = slider_y,
                                 refit = refit)
        }
      } else {
        if (!dynamic){
          carp_heatmap_plot (x,
                             percent = percent,
                             k = k,
                             refit = refit,
                             ...)
        } else {
          carp_dynamic_heatmap_plot(x,
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
    ## This comes up in the default settings for the static when
    ## percent = 0
    plot_frame_full <- plot_frame_full %>% filter(.data$GammaPercent <= max(percent, min(.data$GammaPercent)))
    k <- plot_frame_full %>% summarize(k = min(.data$NCluster)) %>% pull
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
    percent <- get_feature_paths(x, features = character()) %>% filter(.data$NCluster == k) %>%
      select(.data$GammaPercent) %>%
      summarize(percent = mean(.data$GammaPercent)) %>%
      pull
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
             geom_point(data = plot_frame_final, aes(color = .data$final_cluster), size = 4) +
            labs(title = paste0('Fraction of Regularization: ', round(percent * 100), '%\nNumber of Clusters: ', k))

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
#' @importFrom dplyr filter select summarize pull
#' @importFrom stats as.dendrogram cutree
#' @importFrom dendextend set color_branches as.ggdend
#' @importFrom ggplot2 ggplot geom_segment scale_x_continuous scale_y_continuous theme
#' @importFrom ggplot2 element_text xlab ylab aes element_blank labs
carp_dendro_plot <- function(x,
                             percent,
                             k,
                             show_clusters = (n_args == 1L),
                             ...){

  has_percent <- !missing(percent)
  has_k       <- !missing(k)
  n_args      <- has_percent + has_k

  if(n_args > 1L){
    crv_error("At most one of ", sQuote("percent"), " and ", sQuote("k"), " may be supplied.")
  }

  d <- as.dendrogram(x)

  adjust <- get_feature_paths(x, features = character()) %>% filter(.data$NCluster > 1) %>%
    select(.data$GammaPercent) %>%
    summarize(GammaPercent = max(.data$GammaPercent)) %>%
    pull

  if(show_clusters){

    if(has_percent){
      if (!is_percent_scalar(percent)) {
        crv_error(sQuote("percent"), " must be a scalar between 0 and 1 (inclusive).")
      }

      k <- get_feature_paths(x, features = character()) %>% filter(.data$GammaPercent >= percent) %>%
        select(.data$NCluster) %>%
        summarize(NCluster = max(.data$NCluster)) %>%
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
        summarize(percent = mean(.data$GammaPercent)) %>%
        pull
    }

    segs <- adjusted_dendrogram(d, rev = FALSE, k = k, cluster = TRUE, adjust = adjust)
    lines <- dendrogram_box(x, rev = FALSE, k = k, percent = percent)

    p <- ggplot() +
      geom_segment(data = segs, aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend, color = .data$col), show.legend = FALSE) +
      geom_segment(data = lines, aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend, color = NA), show.legend = FALSE) +
      labs(title = paste0('Fraction of Regularization: ', round(percent * 100), '%\nNumber of Clusters: ', k))
    } else {

    segs <- adjusted_dendrogram(d, rev = FALSE, k = k, cluster = FALSE, adjust = adjust)

    p <- ggplot() +
      geom_segment(data = segs, aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend), show.legend = FALSE)
    }

  dend <- as.ggdend(d)
  label <- dend$label$label

  p +
    labs(y='Fraction of Regularization') +
    scale_x_continuous(breaks=seq_along(label),labels= levels(label)) +
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
#' @importFrom stats as.dendrogram cutree
#' @importFrom dendextend set color_branches as.ggdend
#' @importFrom ggplot2 ggplot geom_segment scale_x_continuous scale_y_continuous theme geom_tile scale_fill_gradient2 labs
carp_heatmap_plot <- function(x,
                             percent,
                             k,
                             show_clusters = (n_args == 1L),
                             refit = FALSE,
                             ...){
  has_percent <- !missing(percent)
  has_k       <- !missing(k)
  n_args      <- has_percent + has_k

  if(n_args > 1L){
    crv_error("At most one of ", sQuote("percent"), " and ", sQuote("k"), " may be supplied.")
  }

  d <- as.dendrogram(x)

  adjust <- get_feature_paths(x, features = character()) %>% filter(.data$NCluster > 1) %>%
    select(.data$GammaPercent) %>%
    summarize(GammaPercent = max(.data$GammaPercent)) %>%
    pull

  if(show_clusters){

    if(has_percent){
      if (!is_percent_scalar(percent)) {
        crv_error(sQuote("percent"), " must be a scalar between 0 and 1 (inclusive).")
      }

      k <- get_feature_paths(x, features = character()) %>% filter(.data$GammaPercent >= percent) %>%
        select(.data$NCluster) %>%
        summarize(NCluster = max(.data$NCluster)) %>%
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
        summarize(percent = mean(.data$GammaPercent)) %>%
        pull
    }
    # heatmap
    U <- get_clustered_data(x, percent = percent, refit = refit)

    # dendrogram
    segs <- adjusted_dendrogram(d, rev = TRUE, k = k, cluster = TRUE, adjust = adjust)

    lines <- dendrogram_box(x, rev = TRUE, k = k, percent = percent)

    p <- ggplot() +
      geom_segment(data = segs, aes(x = (.data$y/3+1)*length(cn)+0.5, xend = (.data$yend/3+1)*length(cn)+0.5, y = .data$x, yend = .data$xend, color = .data$col), show.legend = FALSE) +
      geom_segment(data = lines, aes(x = (.data$y/3+1)*length(cn)+0.5, xend = (.data$yend/3+1)*length(cn)+0.5, y = .data$x, yend = .data$xend, color = NA), show.legend = FALSE) +
      labs(title = paste0('Fraction of Regularization: ', round(percent * 100), '%\nNumber of Clusters: ', k))
  } else {
    # heatmap
    U <- get_clustered_data(x, percent = 0, refit = refit)
    k <- NROW(U)

    # dendrogram
    segs <- adjusted_dendrogram(d, rev = TRUE, k = k, cluster = FALSE, adjust = adjust)

    p <- ggplot() +
      geom_segment(data = segs, aes(x = (.data$y/3+1)*length(cn)+0.5, xend = (.data$yend/3+1)*length(cn)+0.5, y =.data$ x, yend = .data$xend), show.legend = FALSE)
  }

  # heatmap
  h <- heatmapr(
    U,
    Rowv = as.hclust(x),
    dendrogram = "row",
    k_row = k
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
    scale_x_continuous(breaks=seq_along(cn),labels= cn) +
    scale_y_continuous(breaks=seq_along(rn),labels= rn) +
    theme(axis.text.x=element_text(hjust=1,vjust=0.5,angle=90),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background = element_blank(),
          panel.grid=element_blank(),
          panel.border=element_blank())
  }

#' @noRd
#' @importFrom rlang .data
#' @importFrom dplyr select filter rename left_join mutate bind_rows
#' @importFrom ggplot2 ggplot aes geom_path geom_point geom_text guides
#' @importFrom ggplot2 theme element_text xlab ylab
#' @importFrom ggrepel geom_text_repel
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
                        mutate(percent = pct*100)
  }))

  ggplot(plot_frame_animation,
         aes(x = .data$V1, y = .data$V2, group = .data$Obs)) +
    geom_path(linejoin = "round", color = "red", size = 1) +
    geom_point(data = plot_frame_first,
               aes(x = .data$FirstV1,
                   y = .data$FirstV2),
               color = "black",
               size = I(2)) +
    geom_text_repel(data = plot_frame_first,
                    aes(x = .data$FirstV1,
                        y = .data$FirstV2,
                        label = .data$FirstObsLabel),
                    seed=0) +
    guides(color = FALSE, size = FALSE) +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 10)) +
    xlab(axis[1]) + ylab(axis[2]) +
    transition_manual(.data$percent) +
    labs(title = paste0('Fraction of Regularization: ', '{current_frame}', '%'))
}

#' @noRd
#' @importFrom rlang .data
#' @importFrom dplyr filter select summarize pull
#' @importFrom stats as.dendrogram cutree
#' @importFrom dendextend set color_branches as.ggdend
#' @importFrom ggplot2 ggplot geom_segment scale_x_continuous scale_y_continuous theme labs
#' @importFrom gganimate transition_time
carp_dynamic_dendro_plot <- function(x,
                             percent.seq,
                             ...){
  d <- as.dendrogram(x)

  adjust <- get_feature_paths(x, features = character()) %>% filter(.data$NCluster > 1) %>%
    select(.data$GammaPercent) %>%
    summarize(GammaPercent = max(.data$GammaPercent)) %>%
    pull

  segs_dynamic <- data.frame()
  lines_dynamic <- data.frame()
  for (per in percent.seq){
    k <- get_feature_paths(x, features = character()) %>% filter(.data$GammaPercent >= per) %>%
      select(.data$NCluster) %>%
      summarize(NCluster = max(.data$NCluster)) %>%
      pull

    # segments
    segs <- adjusted_dendrogram(d, rev = FALSE, k = k, cluster = TRUE, adjust = adjust)
    segs_dynamic <- rbind(segs_dynamic,cbind(segs, reg = per*100))

    # lines
    lines <- dendrogram_box(x, rev = FALSE, k = k, percent = per)
    lines_dynamic <- rbind(lines_dynamic,cbind(lines, reg = per*100))
  }

  dend <- as.ggdend(d)
  label <- dend$label$label

ggplot() +
  geom_segment(data = segs_dynamic, aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend, color = .data$col), show.legend = FALSE) +
  geom_segment(data = lines_dynamic, aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend, color = NA), show.legend = FALSE) +
  labs(y='Fraction of Regularization') +
  scale_x_continuous(breaks=seq_along(label),labels= levels(label)) +
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
#' @importFrom stats as.dendrogram cutree
#' @importFrom dendextend set color_branches as.ggdend
#' @importFrom ggplot2 ggplot geom_segment scale_x_continuous scale_y_continuous theme geom_tile scale_fill_gradient2 labs
#' @importFrom gganimate transition_manual
carp_dynamic_heatmap_plot <- function(x,
                              percent.seq,
                              refit = FALSE,
                              ...){
  adjust <- get_feature_paths(x, features = character()) %>% filter(.data$NCluster > 1) %>%
    select(.data$GammaPercent) %>%
    summarize(GammaPercent = max(.data$GammaPercent)) %>%
    pull

  data_mat_dynamic <- data.frame()
  segs_dynamic <- data.frame()
  lines_dynamic <- data.frame()
  for (per in percent.seq){
    U <- get_clustered_data(x, percent = per, refit = refit)
    k <- nlevels(get_cluster_labels(x, percent = per))
    h <- heatmapr(
      U,
      Rowv = as.hclust(x),
      dendrogram = "row",
      k_row = k
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
    d <- as.dendrogram(x)
    segs <- adjusted_dendrogram(d, rev = TRUE, k = k, cluster = TRUE, adjust = adjust)
    segs_dynamic <- rbind(segs_dynamic,cbind(segs, reg = per*100))

    # lines
    lines <- dendrogram_box(x, rev = TRUE, k = k, percent = per)
    lines_dynamic <- rbind(lines_dynamic,cbind(lines, reg = per*100))
  }

  r <- range(data_mat_dynamic$value)

ggplot(data =  data_mat_dynamic, aes(x = .data$x, y = .data$y)) +
    geom_tile(aes(fill = .data$value)) +
    scale_fill_gradient2(low = "#313695", mid = "#FFFFBF", high = "#A50026", midpoint = (floor(r[1])+ceiling(r[2]))/2, limits = c(floor(r[1]),ceiling(r[2]))) +
    scale_x_continuous(breaks=seq_along(cn),labels= cn) +
    scale_y_continuous(breaks=seq_along(rn),labels= rn) +
    geom_segment(data = segs_dynamic, aes(x = (.data$y/3+1)*length(cn)+0.5, xend = (.data$yend/3+1)*length(cn)+0.5, y = .data$x, yend = .data$xend, color = .data$col), show.legend = FALSE) +
    geom_segment(data = lines_dynamic, aes(x = (.data$y/3+1)*length(cn)+0.5, xend = (.data$yend/3+1)*length(cn)+0.5, y = .data$x, yend = .data$xend, color = NA), show.legend = FALSE) +
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
carp_heatmaply <- function(x,
                           ...,
                           percent,
                           k,
                           refit = FALSE){

  has_percent <- !missing(percent)
  has_k       <- !missing(k)

  n_args <- has_percent + has_k

  if(n_args >= 2){
    crv_error("At most one of ", sQuote("percent,"), " and ",
              sQuote("k"), " may be supplied.")
  }

  if(n_args == 0){
    U     <- get_clustered_data(x, percent = 0, refit = refit)

    heatmaply(U,
              Rowv = as.hclust(x),
              dendrogram = "row",
              ...)
  } else {
    U <- get_clustered_data(x, percent = percent, k = k, refit = refit)
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
#' @importFrom plotly add_markers add_paths add_text plot_ly hide_legend highlight animation_slider style
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

  if(n_args >= 1 & dynamic){
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
#' @importFrom dplyr filter select summarize pull desc
#' @importFrom stats as.dendrogram is.leaf cutree
#' @importFrom dendextend get_nodes_xy as.ggdend
#' @importFrom plotly add_segments add_markers add_text plot_ly hide_legend highlight animation_slider
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

  adjust <- get_feature_paths(x, features = character()) %>% filter(.data$NCluster > 1) %>%
    select(.data$GammaPercent) %>%
    summarize(GammaPercent = max(.data$GammaPercent)) %>%
    pull

  d <- as.dendrogram(x)
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

  tidy_segments <- adjusted_dendrogram (d, rev = FALSE, cluster = FALSE, adjust = adjust)

  allXY$y <- allXY$y*adjust

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

        k <- get_feature_paths(x, features = character()) %>% filter(.data$GammaPercent >= percent) %>%
          select(.data$NCluster) %>%
          summarize(NCluster = max(.data$NCluster)) %>%
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
          summarize(percent = mean(.data$GammaPercent)) %>%
          pull
      }

      segnm <- dendrogram_box(x, rev = FALSE, k = k, percent = percent, show_memnum = TRUE)
      seg <- segnm[[1]]
      m <- segnm[[2]]
      seg <- cbind(seg, perent=percent)

      dendro_static <- plot_ly(
        hoverinfo = "none") %>%
        add_segments(
          data = tidy_segments,
          x = ~x, y = ~y,
          xend = ~xend, yend = ~yend,
          color = I("black"),
          showlegend = FALSE)

      for (i in (1:k)){
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
      k <- get_feature_paths(x, features = character()) %>% filter(.data$GammaPercent >= per) %>%
        select(.data$NCluster) %>%
        summarize(NCluster = max(.data$NCluster)) %>%
        pull

      segnm <- dendrogram_box(x, rev = FALSE, k = k, percent = per, show_memnum = TRUE)
      seg <- segnm[[1]]
      m <- segnm[[2]]
      dynamic_seg <- rbind(dynamic_seg,cbind(seg, Regularization=per))

      allXY_per <- allXY[allXY$y > 0,]
      allXY_per$label <- seq_along(allXY_per$label)

      allTXT_per <- allTXT

      for (i in (1:k)){
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
carp_heatmaply_dynamic <- function(x,
                                   ...,
                                   percent.seq = percent.seq,
                                   slider_y = slider_y,
                                   refit = refit){
  adjust <- get_feature_paths(x, features = character()) %>% filter(.data$NCluster > 1) %>%
    select(.data$GammaPercent) %>%
    summarize(GammaPercent = max(.data$GammaPercent)) %>%
    pull

  data_mat_dynamic <- data.frame()
  seg_dynamic <- data.frame()
  for (per in percent.seq){
    U <- get_clustered_data(x, percent = per, refit = refit)
    k <- nlevels(get_cluster_labels(x, percent = per))
    h <- heatmapr(
      U,
      Rowv = as.hclust(x),
      dendrogram = "row",
      k_row = k
    )
    # data for heatmap
    data_mat <- h$matrix$data
    data_mat_dynamic <- rbind(data_mat_dynamic,cbind(value = as.vector(data_mat),
                                                     x = as.vector(col(data_mat)),
                                                     y = as.vector(row(data_mat)),
                                                     percent = as.vector(per)))

    # data for dendrogram
    dend <- h$rows
    segs <- as.ggdend(dend)$segments
    segs$col[is.na(segs$col)] <- "black" # default value for NA is "black"

    segs$y <- segs$y*adjust
    segs$yend <- segs$yend*adjust

    seg_dynamic <- rbind(seg_dynamic,cbind(seg = seq_along(segs$x), segs, percent = per))

    rn <- rownames(h$matrix$data)
    cn <- colnames(h$matrix$data)
  }

  colors <- sort(unique(seg_dynamic$col))

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
        range = c(0.5, length(rn) + 0.5),
        showticklabels = TRUE,
        showgrid = FALSE))

  for (i in seq_along(levels(factor(seg_dynamic$seg)))){
    p <- p %>%
      add_segments(data = seg_dynamic[seg_dynamic$seg==i,],
                   x = ~(y/3+1)*length(cn)+0.5, xend = ~(yend/3+1)*length(cn)+0.5, y = ~x, yend = ~xend, color = ~col,
                   showlegend = FALSE,
                   colors = colors,
                   hoverinfo = "none",
                   frame = ~percent*100)
  }

  p %>%
    animation_slider(y = slider_y,
                     currentvalue = list(prefix = "Regularization: ", suffix = "%"))
}

