context("Test Clustering Assignments")

test_that("Star graph clustering works", {
  get_cluster_assignments <- clustRviz:::get_cluster_assignments
  n <- 5
  E <- matrix(c(1, 2, ## This checks fusing beginnings to beginnings
                1, 3,
                1, 4,
                1, 5), ncol=2, byrow=TRUE)

  ## Fully clustered
  edge_indicator <- matrix(c(1, 1, 1, 1), nrow=1)
  cluster_assignment <- get_cluster_assignments(E, edge_indicator, n)[[1]]
  expect_equal(cluster_assignment$no, 1)
  expect_equal(cluster_assignment$csize, 5)
  expect_equal(cluster_assignment$membership, rep(1, 5))

  ## No clustering
  edge_indicator <- matrix(c(0, 0, 0, 0), nrow=1)
  cluster_assignment <- get_cluster_assignments(E, edge_indicator, n)[[1]]
  expect_equal(cluster_assignment$no, 5)
  expect_equal(cluster_assignment$csize, rep(1, 5))
  expect_equal(cluster_assignment$membership, seq(1, 5))

  ## Partial clustering
  edge_indicator <- matrix(c(1, 1, 0, 0), nrow=1)
  cluster_assignment <- get_cluster_assignments(E, edge_indicator, n)[[1]]
  expect_equal(cluster_assignment$no, 3)
  expect_equal(cluster_assignment$csize, c(3, 1, 1))
  expect_equal(cluster_assignment$membership, c(1, 1, 1, 2, 3))
})

test_that("Chain graph clustering works", {
  get_cluster_assignments <- clustRviz:::get_cluster_assignments
  n <- 5
  E <- matrix(c(1, 2, ## This checks fusing beginnings to ends
                2, 3,
                3, 4,
                4, 5), ncol=2, byrow=TRUE)

  ## Fully clustered
  edge_indicator <- matrix(c(1, 1, 1, 1), nrow=1)
  cluster_assignment <- get_cluster_assignments(E, edge_indicator, n)[[1]]
  expect_equal(cluster_assignment$no, 1)
  expect_equal(cluster_assignment$csize, 5)
  expect_equal(cluster_assignment$membership, rep(1, 5))

  ## No clustering
  edge_indicator <- matrix(c(0, 0, 0, 0), nrow=1)
  cluster_assignment <- get_cluster_assignments(E, edge_indicator, n)[[1]]
  expect_equal(cluster_assignment$no, 5)
  expect_equal(cluster_assignment$csize, rep(1, 5))
  expect_equal(cluster_assignment$membership, seq(1, 5))

  ## Partial clustering
  edge_indicator <- matrix(c(1, 0, 0, 1), nrow=1)
  cluster_assignment <- get_cluster_assignments(E, edge_indicator, n)[[1]]
  expect_equal(cluster_assignment$no, 3)
  expect_equal(cluster_assignment$csize, c(2, 1, 2))
  expect_equal(cluster_assignment$membership, c(1, 1, 2, 3, 3))
})

test_that("Disconnected graph doesn't cluster", {
  # This shouldn't happen in clustRviz, but our code can handle it
  get_cluster_assignments <- clustRviz:::get_cluster_assignments
  n <- 6
  E <- matrix(c(1, 2,
                2, 3,
                4, 5,
                5, 6), ncol=2, byrow=TRUE)

  edge_indicator <- matrix(c(1, 1, 1, 1), nrow=1)
  cluster_assignment <- get_cluster_assignments(E, edge_indicator, n)[[1]]
  expect_equal(cluster_assignment$no, 2)
  expect_equal(cluster_assignment$csize, c(3, 3))
  expect_equal(cluster_assignment$membership, c(1, 1, 1, 2, 2, 2))
})

test_that("Clusters get merged when only ends are shared", {
  # This shouldn't happen in clustRviz, but our code can handle it
  get_cluster_assignments <- clustRviz:::get_cluster_assignments
  n <- 5
  E <- matrix(c(1, 5,
                2, 5,
                3, 4), ncol=2, byrow=TRUE)

  edge_indicator <- matrix(c(1, 1, 1), nrow=1)
  cluster_assignment <- get_cluster_assignments(E, edge_indicator, n)[[1]]
  expect_equal(cluster_assignment$no, 2)
  expect_equal(cluster_assignment$csize, c(3, 2))
  expect_equal(cluster_assignment$membership, c(1, 1, 2, 2, 1))

  ## Add some singletons
  n <- 7
  cluster_assignment <- get_cluster_assignments(E, edge_indicator, n)[[1]]
  expect_equal(cluster_assignment$no, 4)
  expect_equal(cluster_assignment$csize, c(3, 2, 1, 1))
  expect_equal(cluster_assignment$membership, c(1, 1, 2, 2, 1, 3, 4))
})

test_that("Functionality works row-wise", {
  set.seed(25)
  get_cluster_assignments <- clustRviz:::get_cluster_assignments

  n <- 50
  k <- 20
  E <- matrix(sample(n, 2 * k, replace = TRUE), ncol=2)
  E <- unique(E)

  edge_indicator <- matrix(rbinom(NROW(E) * 5, size=1, prob=0.5), ncol=NROW(E))

  expect_equal(get_cluster_assignments(E, edge_indicator, n),
               lapply(1:5, function(i) get_cluster_assignments(E, edge_indicator[i, , drop=FALSE], n)[[1]]))
})

test_that("Cluster labels are increasing", {
  ## For stability, we require the cluster labels to be assigned in increasing order
  ##
  ## That is, if we get the smallest index for each cluster, those indices should be
  ## increasing
  set.seed(125)
  get_cluster_assignments <- clustRviz:::get_cluster_assignments

  n <- 40
  k <- 30
  E <- matrix(sample(n, 2 * k, replace = TRUE), ncol=2)
  E <- unique(E)

  edge_indicator <- matrix(rbinom(NROW(E) * 50, size=1, prob=0.5), ncol=NROW(E))

  cluster_assignments <- get_cluster_assignments(E, edge_indicator, n)

  is_increasing <- function(x) all(x == cummax(x))
  unique_increasing <- function(x) is_increasing(unique(x))
  expect_true(all(vapply(cluster_assignments, function(x) unique_increasing(x$membership), logical(1))))
})

test_that("Cluster labels are sequential (no gaps)", {
  set.seed(500)
  get_cluster_assignments <- clustRviz:::get_cluster_assignments

  n <- 400
  k <- 200
  E <- matrix(sample(n, 2 * k, replace = TRUE), ncol=2)
  E <- unique(E)

  edge_indicator <- matrix(rbinom(NROW(E) * 50, size=1, prob=0.5), ncol=NROW(E))

  cluster_assignments <- get_cluster_assignments(E, edge_indicator, n)

  has_unique_labels <- function(x){
    num_unique <- length(unique(x))
    all(seq_len(num_unique) %in% x) && all(x %in% seq_len(num_unique))
  }
  expect_true(all(vapply(cluster_assignments, function(x) has_unique_labels(x$membership), logical(1))))
})

test_that("Number of clusters is correctly calculated", {
  set.seed(750)
  get_cluster_assignments <- clustRviz:::get_cluster_assignments

  n <- 400
  k <- 200
  E <- matrix(sample(n, 2 * k, replace = TRUE), ncol=2)
  E <- unique(E)

  edge_indicator <- matrix(rbinom(NROW(E) * 50, size=1, prob=0.75), ncol=NROW(E))

  cluster_assignments <- get_cluster_assignments(E, edge_indicator, n)

  expect_true(all(vapply(cluster_assignments,
                         function(x) (x$no == length(unique(x$membership))) && (x$no == length(x$csize)),
                         logical(1))))
})
