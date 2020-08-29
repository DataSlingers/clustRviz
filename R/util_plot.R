# Adjust the height of dendrograms
#' @importFrom dendextend color_branches as.ggdend
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

# Generate box for dendrograms
#' @importFrom stats as.hclust cutree
dendrogram_box <- function(x, rev = FALSE, k, type, percent, show_memnum = FALSE){
  tree <- as.hclust(x, type = type)
  if (rev == TRUE) {
    tree <- rev(tree)
  }
  cluster <- cutree(tree, k = k)
  clustab <- table(cluster)[unique(cluster[tree$order])]
  clustsum <- cumsum(clustab)
  m <- c(0, clustsum) + 0.5
  line_x <- c(m[1],m[1],m)
  line_y <- c(0,percent,rep(percent,k+1))
  line_xend <- c(m[k+1],m[k+1],m)
  line_yend <- c(0,percent,rep(0,k+1))
  lines <- data.frame(x=line_x,y=line_y,xend=line_xend,yend=line_yend)

  if (show_memnum == FALSE){
    return(lines)
  } else {
    return(list(lines, m))
  }
}


# Manually "dodge" the labels for path plots
# from https://stackoverflow.com/questions/45065567/getting-coordinates-for-the-label-locations-from-ggrepel
#' @importFrom ggplot2 ggplot geom_point ggplot_build
#' @importFrom ggrepel geom_text_repel
#' @importFrom grid grid.force grid.get childNames forceGrob
get_ggrepel_coordinates <- function(x, y, label){
  df <- data.frame(x = x, y = y, label = label)
  p <- ggplot(data = df, aes(x = x, y = y)) +
    geom_text_repel(aes(label = label), size = 3) +
    geom_point()

  # Get x and y plot ranges
  xrg <- ggplot_build(p)$layout$panel_params[[1]]$x.range
  yrg <- ggplot_build(p)$layout$panel_params[[1]]$y.range

  forceGrob(p)
  grid.force()
  kids <- childNames(grid.get("textrepeltree", grep = TRUE))

  # Get positions of all ggrepel labels
  dts <- do.call(rbind, lapply(kids, get.xy.pos.labs, xrg = xrg, yrg = yrg))
  colnames(dts) <- c("x_adj", "y_adj")
  return(cbind(df,dts))
}

# get the x and y positions of a single ggrepel label
#' @importFrom grid grid.get convertX convertY
get.xy.pos.labs <- function(n, xrg, yrg) {
  grb <- grid.get(n)
  data.frame(
    x = xrg[1]+diff(xrg)*convertX(grb$x, "native", valueOnly = TRUE),
    y = yrg[1]+diff(yrg)*convertY(grb$y, "native", valueOnly = TRUE)
  )
}
