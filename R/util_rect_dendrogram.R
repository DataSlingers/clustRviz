#####
# This function modifies the rect.dendrogram function of the dendextend
# package: https://github.com/talgalili/dendextend/

# Modifications are the inclusion of the my.col.vec arguement to
# allow for colored/shaded dendrogram rectangles.
#
# Author of Changes: John Nagorski
# Date: 9-11-18
#' @importFrom dendextend is.dendrogram
#' @importFrom dendextend heights_per_k.dendrogram
#' @importFrom stats order.dendrogram
#' @importFrom stats cutree
#' @importFrom graphics par
#' @importFrom graphics grconvertX
#' @importFrom graphics grconvertY
#' @importFrom graphics rect
rect.dendrogram <- function(tree, k = NULL, which = NULL, x = NULL, h = NULL, border = 2,
                            cluster = NULL, horiz = FALSE, density = NULL, angle = 45,
                            text = NULL, text_cex = 1, text_col = 1, xpd = TRUE, lower_rect,
                            upper_rect = 0, prop_k_height = 0.5, stop_if_out = FALSE,
                            my.col.vec = NULL,
                            ...) {
  if (!dendextend::is.dendrogram(tree)) {
    crv_error("x is not a dendrogram object.")
  }
  if (length(h) > 1L | length(k) > 1L) {
    crv_error("'k' and 'h' must be a scalar(i.e.: of length 1)")
  }
  tree_heights <- dendextend::heights_per_k.dendrogram(tree)[-1]
  tree_order <- stats::order.dendrogram(tree)
  if (!is.null(h)) {
    if (!is.null(k)) {
      crv_error("specify exactly one of 'k' and 'h'")
    }
    ss_ks <- tree_heights < h
    k <- min(as.numeric(names(ss_ks))[ss_ks])
    k <- max(k, 2)
  }
  else if (is.null(k)) {
    crv_error("specify exactly one of 'k' and 'h'")
  }
  if (k < 2 | k > length(tree_heights)) {
    if (stop_if_out) {
      crv_error(gettextf("k must be between 2 and %d", length(tree_heights), domain = NA))
    }
    else {
      crv_warning(gettextf("k must be between 2 and %d", length(tree_heights), domain = NA))
    }
  }
  if (is.null(cluster)) {
    cluster <- stats::cutree(tree, k = k)
  }
  clustab <- table(cluster)[unique(cluster[tree_order])]
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
    crv_error(gettextf("all elements of 'which' must be between 1 and %d", k, domain = NA))
  }
  border <- rep_len(border, length(which))
  retval <- list()
  old_xpd <- graphics::par()["xpd"]
  graphics::par(xpd = xpd)
  my.col.vec <- rep(my.col.vec, length.out = length(which))
  for (n in seq_along(which)) {
    next_k_height <- tree_heights[names(tree_heights) ==
                                    k + 1]
    if (length(next_k_height) == 0) {
      next_k_height <- 0
      prop_k_height <- 1
    }
    if (!horiz) {
      xleft <- m[which[n]] + 0.66
      if (missing(lower_rect)) {
        lower_rect <- graphics::par("usr")[3L] - graphics::strheight("W") *
          (max(nchar(labels(tree))) + 1)
      }
      ybottom <- lower_rect
      xright <- m[which[n] + 1] + 0.33
      ytop <- tree_heights[names(tree_heights) == k] *
        prop_k_height + next_k_height * (1 - prop_k_height) +
        upper_rect
    }
    else {
      ybottom <- m[which[n]] + 0.66
      if (missing(lower_rect)) {
        lower_rect <- graphics::par("usr")[2L] + graphics::strwidth("X") *
          (max(nchar(labels(tree))) + 1)
      }
      xright <- lower_rect
      ytop <- m[which[n] + 1] + 0.33
      xleft <- tree_heights[names(tree_heights) == k] *
        prop_k_height + next_k_height * (1 - prop_k_height) +
        upper_rect
    }
    graphics::rect(xleft, ybottom, xright, ytop,
                   border = border[n],

                   density = density, angle = angle, my.col.vec[n], ...
    )
    if (!is.null(text)) {
      graphics::text((m[which[n]] + m[which[n] + 1] + 1) / 2, graphics::grconvertY(graphics::grconvertY(
        graphics::par("usr")[3L],
        "user", "ndc"
      ) + 0.02, "ndc", "user"), text[n],
      cex = text_cex, col = text_col
      )
    }
    retval[[n]] <- which(cluster == as.integer(names(clustab)[which[n]]))
  }
  graphics::par(xpd = old_xpd)
  invisible(retval)
}
