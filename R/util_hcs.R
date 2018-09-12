#####
# These function modify the iorder and hc functions of the
# hclust_in_R package: https://github.com/bwlewis/hclust_in_R

# Modifications are custom naming of the heirarchical clustering
# method.
#
# Author of Changes: John Nagorski
# Date: 9-11-18
iorder <- function(m) {
  N <- nrow(m) + 1
  iorder <- rep(0, N)
  iorder[1] <- m[N - 1, 1]
  iorder[2] <- m[N - 1, 2]
  loc <- 2
  for (i in seq(N - 2, 1))
  {
    for (j in seq(1, loc))
    {
      if (iorder[j] == i) {
        iorder[j] <- m[i, 1]
        if (j == loc) {
          loc <- loc + 1
          iorder[loc] <- m[i, 2]
        } else {
          loc <- loc + 1
          for (k in seq(loc, j + 2)) iorder[k] <- iorder[k - 1]
          iorder[j + 1] <- m[i, 2]
        }
      }
    }
  }
  -iorder
}
cvxhc <- function(clust.path, gamma.path, labels) {
  h <- sort(gamma.path, decreasing = FALSE)
  # h = log(h + 1)
  h <- h[-1]
  N <- length(clust.path)
  n <- -(1:N) # Tracks group membership
  m <- matrix(0, nrow = N - 1, ncol = 2) # hclust merge output
  for (j in seq(1, N - 1))
  {

    # affected observation indicies
    old.clustering <- clust.path[[j]]
    new.clustering <- clust.path[[j + 1]]
    i <- new.clustering[!(new.clustering %in% old.clustering)][[1]]

    i
    p <- n[i]
    p
    n.neg <- sum(p < 0)
    n.pos <- sum(p > 0)
    n.neg
    n.pos
    # if n.neg > 0 and n.pos ==0, joining two singletons
    if ((n.neg > 0) & (n.pos == 0)) {
      i
      p <- n[i]
      p
    } else if ((n.neg > 0) & (n.pos > 0)) {
      i <- c(i[p < 0][1], i[p > 0][1])
      i
      p <- n[i]
    } else if ((n.neg == 0) & (n.pos > 0)) {
      tmp.clusts <- unique(p)
      i <- c(i[p == tmp.clusts[1]][1], i[p == tmp.clusts[2]][1])
      p <- n[i]
    }
    p

    # R's convention is to order each m[j,] pair as follows:
    p <- p[order(p)]

    m[j, ] <- p
    # Agglomerate this pair and all previous groups they belong to
    # into the current jth group:
    grp <- c(i, which(n %in% n[i[n[i] > 0]]))
    n[grp] <- j
  }
  structure(list(
    merge = m, height = h, order = iorder(m),
    labels = labels, method = "CVX",
    call = match.call(), dist.method = "euclidean"
  ),
  class = "hclust"
  )
}
