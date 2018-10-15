#####
# This function modifies the heatmap.2 function of the gplots
# package: https://github.com/cran/gplots

# Modifications allow the inclusion of of colored/shaded
# dendrogram rectangles; see util_rect_dendrogram.R
#
# Author of Changes: John Nagorski
# Date: 9-11-18
#' @importFrom stats dist hclust median reorder density sd as.dendrogram order.dendrogram
#' @importFrom graphics abline axis hist image layout lines mtext plot plot.new
#' @importFrom graphics strheight strwidth text title
#' @importFrom gtools invalid
my.heatmap.2 <- function(x, Rowv = TRUE,
                         Colv = if (symm) "Rowv" else TRUE,
                         distfun = stats::dist,
                         hclustfun = stats::hclust,
                         dendrogram = c("both", "row", "column", "none"),
                         symm = FALSE, scale = c("none", "row", "column"),
                         na.rm = TRUE, revC = identical(Colv, "Rowv"),
                         add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) || scale != "none",
                         col = "heat.colors", colsep, rowsep,
                         sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote,
                         notecex = 1, notecol = "cyan", na.color = par("bg"),
                         trace = c("column", "row", "both", "none"),
                         tracecol = "cyan", hline = stats::median(breaks),
                         vline = stats::median(breaks), linecol = tracecol,
                         margins = c(5, 5),
                         ColSideColors, RowSideColors, cexRow = 0.2 + 1 / log10(nr),
                         cexCol = 0.2 + 1 / log10(nc), labRow = NULL, labCol = NULL,
                         key = TRUE, keysize = 1.5,
                         density.info = c("histogram", "density", "none"),
                         denscol = tracecol, symkey = min(x < 0, na.rm = TRUE) || symbreaks,
                         densadj = 0.25, main = NULL, xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL,
                         k.row = 3, k.col = 4, border = "red",
                         Row.hclust = NULL,
                         Col.hclust = NULL, my.col.vec = NULL, lwd = 2, ...) {
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low) / (high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale)) {
    "none"
  } else {
    match.arg(scale)
  }
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col)) {
    col <- get(col, mode = "function")
  }
  if (!missing(breaks) && (scale != "none")) {
    crv_warning(
      "Using scale=\"row\" or scale=\"column\" when breaks are",
      "specified can produce unpredictable results.",
      "Please consider using only one or the other."
    )
  }
  if (is.null(Rowv) || is.na(Rowv)) {
    Rowv <- FALSE
  }
  if (is.null(Colv) || is.na(Colv)) {
    Colv <- FALSE
  } else if (Colv == "Rowv" && !isTRUE(Rowv)) {
    Colv <- FALSE
  }
  if (length(di <- dim(x)) != 2 || !is.numeric(x)) {
    crv_error("`x' must be a numeric matrix")
  }
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1) {
    crv_error("`x' must have at least 2 rows and 2 columns")
  }
  if (!is.numeric(margins) || length(margins) != 2) {
    crv_error("`margins' must be a numeric vector of length 2")
  }
  if (missing(cellnote)) {
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  }
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv)) {
        dendrogram <- "column"
      } else {
        dedrogram <- "none"
      }
      crv_warning(
        "Discrepancy: Rowv is FALSE, while dendrogram is `",
        dendrogram, "'. Omitting row dendogram."
      )
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv)) {
        dendrogram <- "row"
      } else {
        dendrogram <- "none"
      }
      crv_warning(
        "Discrepancy: Colv is FALSE, while dendrogram is `",
        dendrogram, "'. Omitting column dendogram."
      )
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- stats::as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd)) {
      crv_error("row dendrogram ordering gave index of wrong length")
    }
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- stats::as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd)) {
      crv_error("row dendrogram ordering gave index of wrong length")
    }
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc) {
      crv_error("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    }
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else {
      colInd <- rowInd
    }
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm) {
      x
    } else {
      t(x)
    } ))
    ddc <- stats::as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd)) {
      crv_error("column dendrogram ordering gave index of wrong length")
    }
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm) {
      x
    } else {
      t(x)
    } ))
    ddc <- stats::as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd)) {
      crv_error("column dendrogram ordering gave index of wrong length")
    }
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow)) {
    labRow <- if (is.null(rownames(x))) {
      (1:nr)[rowInd]
    } else {
      rownames(x)
    }
  } else {
    labRow <- labRow[rowInd]
  }
  if (is.null(labCol)) {
    labCol <- if (is.null(colnames(x))) {
      (1:nc)[colInd]
    } else {
      colnames(x)
    }
  } else {
    labCol <- labCol[colInd]
  }
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, stats::sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, stats::sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) <
      1) {
    if (missing(col) || is.function(col)) {
      breaks <- 16
    } else {
      breaks <- length(col) + 1
    }
  }
  if (length(breaks) == 1) {
    if (!symbreaks) {
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks
      )
    } else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function") {
    col <- col(ncol)
  }
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei)) {
    lhei <- c(keysize, 4)
  }
  if (missing(lwid) || is.null(lwid)) {
    lwid <- c(keysize, 4)
  }
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    if (!missing(ColSideColors)) {
      if (!is.character(ColSideColors) || length(ColSideColors) !=
          nc) {
        crv_error("'ColSideColors' must be a character vector of length ncol(x)")
      }
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] +
                      1)
      lhei <- c(lhei[1], 0.2, lhei[2])
    }
    if (!missing(RowSideColors)) {
      if (!is.character(RowSideColors) || length(RowSideColors) !=
          nr) {
        crv_error("'RowSideColors' must be a character vector of length nrow(x)")
      }
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) -
                                           1), 1), lmat[, 2] + 1)
      lwid <- c(lwid[1], 0.2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  if (length(lhei) != nrow(lmat)) {
    crv_error("lhei must have length = nrow(lmat) = ", nrow(lmat))
  }
  if (length(lwid) != ncol(lmat)) {
    crv_error("lwid must have length = ncol(lmat) =", ncol(lmat))
  }
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  graphics::layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  if (!missing(RowSideColors)) {
    par(mar = c(margins[1], 0, 0, 0.5))
    graphics::image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
  }
  if (!missing(ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2]))
    graphics::image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
  }
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr")) {
      ddr <- rev(ddr)
    }
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else {
    iy <- 1:nr
  }
  graphics::image(1:nc, 1:nr, x,
                  xlim = 0.5 + c(0, nc), ylim = 0.5 +
                    c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col,
                  breaks = breaks, ...
  )
  retval$carpet <- x
  if (exists("ddr")) {
    retval$rowDendrogram <- ddr
  }
  if (exists("ddc")) {
    retval$colDendrogram <- ddc
  }
  retval$breaks <- breaks
  retval$col <- col
  if (!gtools::invalid(na.color) & any(is.na(x))) {
    mmat <- ifelse(is.na(x), 1, NA)
    graphics::image(1:nc, 1:nr, mmat,
                    axes = FALSE, xlab = "", ylab = "",
                    col = na.color, add = TRUE
    )
  }
  graphics::axis(1, 1:nc,
                 labels = labCol, las = 2, line = -0.5, tick = 0,
                 cex.axis = cexCol
  )
  if (!is.null(xlab)) {
    graphics::mtext(xlab, side = 1, line = margins[1] - 1.25)
  }
  graphics::axis(4, iy,
                 labels = labRow, las = 2, line = -0.5, tick = 0,
                 cex.axis = cexRow
  )
  if (!is.null(ylab)) {
    graphics::mtext(ylab, side = 4, line = margins[2] - 1.25)
  }
  if (!missing(add.expr)) {
    eval(substitute(add.expr))
  }
  if (!missing(colsep)) {
    for (csep in colsep) rect(
      xleft = csep + 0.5, ybottom = 0,
      xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) +
        1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor
    )
  }
  if (!missing(rowsep)) {
    for (rsep in rowsep) rect(
      xleft = 0, ybottom = (ncol(x) +
                              1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) +
                                                                               1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1,
      col = sepcolor, border = sepcolor
    )
  }
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        graphics::abline(
          v = i - 0.5 + vline.vals, col = linecol,
          lty = 2
        )
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      graphics::lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        graphics::abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      graphics::lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote)) {
    graphics::text(
      x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
      col = notecol, cex = notecex
    )
  }
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    graphics::plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    tree <- Row.hclust
    cluster <- stats::cutree(tree, k = k.row)
    cluster
    clustab <- table(cluster)[unique(cluster[tree$order])]
    clustab
    m <- c(0, cumsum(clustab))
    m
    which <- 1:k.row
    which
    retval <- list()
    my.col.vec.row <- rep(my.col.vec, length.out = length(which))
    if (k.row == 1) {
      rect(
        xleft = par("usr")[1L],
        xright = par("usr")[2L],
        ybottom = par("usr")[3L],
        ytop = par("usr")[4L],
        border = "red",
        col = my.col.vec[1],
        lwd = lwd
      )
    } else {
      for (n in seq_along(which)) {
        rect(par("usr")[2L], m[which[n]] + 0.66,
             mean(rev(tree$height)[(k.row - 1):k.row]),
             m[which[n] + 1] + 0.33,
             border = "red",
             col = my.col.vec.row[n],
             lwd = lwd
        )
      }
    }
  }
  else {
    graphics::plot.new()
  }
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    graphics::plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    tree <- Col.hclust
    cluster <- stats::cutree(tree, k = k.col)
    cluster
    clustab <- table(cluster)[unique(cluster[tree$order])]
    clustab
    m <- c(0, cumsum(clustab))
    m
    which <- 1:k.col
    retval <- list()
    my.col.vec.col <- rep(my.col.vec, length.out = length(which))
    if (k.col == 1) {
      rect(
        xleft = par("usr")[1L],
        xright = par("usr")[2L],
        ybottom = par("usr")[3L],
        ytop = par("usr")[4L], border = "red", col = my.col.vec[1], lwd = lwd
      )
    } else {
      for (n in seq_along(which)) {
        rect(m[which[n]] + 0.66, par("usr")[3L], m[which[n] +
                                                     1] + 0.33,
             mean(rev(tree$height)[(k.col - 1):k.col]),
             border = "red",
             col = my.col.vec.col[n], lwd = lwd
        )
      }
    }
  }
  else {
    graphics::plot.new()
  }
  if (!is.null(main)) {
    graphics::title(main, cex.main = 1.5 * op[["cex.main"]])
  }
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    z <- seq(min.raw, max.raw, length = length(col))
    graphics::image(
      z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
      xaxt = "n", yaxt = "n"
    )
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    graphics::axis(1, at = xv, labels = lv)
    if (scale == "row") {
      graphics::mtext(side = 1, "Row Z-Score", line = 2)
    } else if (scale == "column") {
      graphics::mtext(side = 1, "Column Z-Score", line = 2)
    } else {
      graphics::mtext(side = 1, "Value", line = 2)
    }
    if (density.info == "density") {
      dens <- stats::density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      graphics::lines(dens$x, dens$y / max(dens$y) * 0.95,
                      col = denscol,
                      lwd = 1
      )
      graphics::axis(2,
                     at = pretty(dens$y) / max(dens$y) * 0.95,
                     pretty(dens$y)
      )
      graphics::title("Color Key\nand Density Plot")
      par(cex = 0.5)
      graphics::mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- graphics::hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      graphics::lines(hx, hy / max(hy) * 0.95,
                      lwd = 1, type = "s",
                      col = denscol
      )
      graphics::axis(2, at = pretty(hy) / max(hy) * 0.95, pretty(hy))
      graphics::title("Color Key\nand Histogram")
      par(cex = 0.5)
      graphics::mtext(side = 2, "Count", line = 2)
    }
    else {
      graphics::title("Color Key")
    }
  }
  else {
    graphics::plot.new()
  }
  retval$colorTable <- data.frame(
    low = retval$breaks[-length(retval$breaks)],
    high = retval$breaks[-1], color = retval$col
  )
  invisible(retval)
}
