# These functions are to help generate the dendrogram

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
