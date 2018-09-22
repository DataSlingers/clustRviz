#####
# This function modifies the rect.hclust function of the R stats
# package.

# Modifications are the inclusion of the my.col.vec arguement to
# allow for colored/shaded dendrogram rectangles.
#
# Author of Changes: John Nagorski
# Date: 9-11-18

#' @importFrom graphics rect
#' @importFrom graphics par
#' @importFrom stats cutree
my.rect.hclust <- function(tree, k = NULL, which = NULL, x = NULL, h = NULL, border = 2,
                           cluster = NULL, my.col.vec = NULL, lwd = NULL) {
  if (length(h) > 1L | length(k) > 1L) {
    crv_error("'k' and 'h' must be a scalar")
  }
  if (!is.null(h)) {
    if (!is.null(k)) {
      crv_error("specify exactly one of 'k' and 'h'")
    }
    k <- min(which(rev(tree$height) < h))
    k <- max(k, 2)
  }
  else if (is.null(k)) {
    crv_error("specify exactly one of 'k' and 'h'")
  }
  if (k < 2 | k > length(tree$height)) {
    graphics::rect(
      xleft = graphics::par("usr")[1L],
      xright = graphics::par("usr")[2L],
      ybottom = graphics::par("usr")[3L],
      ytop = graphics::par("usr")[4L],
      border = "red",
      col = my.col.vec[1],
      lwd = lwd
    )
  } else {
    if (is.null(cluster)) {
      cluster <- stats::cutree(tree, k = k)
    }
    clustab <- table(cluster)[unique(cluster[tree$order])]
    m <- c(0, cumsum(clustab))
    if (!is.null(x)) {
      if (!is.null(which)) {
        crv_error("specify exactly one of 'which' and 'x'")
      }
      which <- x
      for (n in seq_along(x)) which[n] <- max(which(m < x[n]))
    }
    else if (is.null(which)) {
      which <- 1L:k
    }
    if (any(which > k)) {
      crv_error(gettextf("all elements of 'which' must be between 1 and %d",
                         k, domain = NA))
    }
    border <- rep_len(border, length(which))
    #### ADD
    my.col.vec <- rep(my.col.vec, length.out = length(which))
    ####
    retval <- list()
    for (n in seq_along(which)) {
      graphics::rect(m[which[n]] + 0.66, par("usr")[3L], m[which[n] +
                                                             1] + 0.33,
                     mean(rev(tree$height)[(k - 1):k]),
                     border = border[n],
                     col = my.col.vec[n],
                     lwd = lwd
      )
      retval[[n]] <- which(cluster == as.integer(names(clustab)[which[n]]))
    }
  }

  # invisible(retval)
}
