context("Test convex_clustering() Solver")

test_that("convex_clustering() errors early with incorrect input", {
  # Pre-processing parameters must be boolean flags
  expect_error(convex_clustering(presidential_speech, X.center = NA, lambda_grid = 1:5))
  expect_error(convex_clustering(presidential_speech, X.center = c(TRUE, FALSE), lambda_grid = 1:5))

  expect_error(convex_clustering(presidential_speech, X.scale = NA, lambda_grid = 1:5))
  expect_error(convex_clustering(presidential_speech, X.scale = c(TRUE, FALSE), lambda_grid = 1:5))

  # Fail on unknown flags
  expect_error(convex_clustering(presidential_speech, flag="unknown", lambda_grid = 1:5), regexp = "flag")
  expect_error(convex_clustering(presidential_speech, "value", lambda_grid = 1:5), regexp = "Unknown")

  ps <- presidential_speech
  ps[1,1] <- Inf; expect_error(convex_clustering(ps, lambda_grid = 1:5))

  ## lambda_grid must be supplied, strictly positive, and ordered (only gets a warning)
  expect_error(convex_clustering(presidential_speech))
  expect_error(convex_clustering(presidential_speech), lambda_grid = numeric())
  expect_error(convex_clustering(presidential_speech, lambda_grid = 0))
  expect_error(convex_clustering(presidential_speech, lambda_grid = c(0, 3)))
  expect_warning(convex_clustering(presidential_speech, lambda_grid = c(3, 2, 1)))
})

test_that("convex_clustering() matches cvxclustr", {
  skip_on_cran()
  skip_if_not_installed("cvxclustr")
  data(mammals, package = "cvxclustr")

  ## Example modified from help pages of cvxclustr
  X <- as.matrix(mammals[,-1])
  clust_fit <- convex_clustering(X, X.center = FALSE, X.scale = FALSE,
                                 lambda_grid = seq(4, 40, length.out = 10))

  ## Calculate matching `cvxclustr` solution
  Xt <- t(X)

  ## Match CARP() selected weights
  w <- clustRviz:::weight_mat_to_vec(clust_fit$weights)

  ## Perform clustering
  suppressWarnings(cvxclust_fit <- cvxclustr::cvxclust(Xt, w, tol = 1e-7,
                                                       gamma = clust_fit$lambda_grid))

  ## cvxclustr seems to use a pretty loose stopping tolerance, so this is a loose check...
  for(i in seq_along(clust_fit$lambda_grid)){
    expect_equal(clust_fit$U[,,i], t(cvxclust_fit$U[[i]]),
                 check.attributes = FALSE, tolerance = 1e-4)
  }
})
