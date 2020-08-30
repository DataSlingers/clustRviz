context("Test convex_biclustering() Solver")

test_that("convex_biclustering() errors early with incorrect input", {
  # Pre-processing parameters must be boolean flags
  expect_error(convex_biclustering(presidential_speech, X.center.global = NA, lambda_grid = 1:5))
  expect_error(convex_biclustering(presidential_speech, X.center.global = c(TRUE, FALSE), lambda_grid = 1:5))

  # Fail on unknown flags
  expect_error(convex_biclustering(presidential_speech, flag="unknown", lambda_grid = 1:5), regexp = "flag")
  expect_error(convex_biclustering(presidential_speech, "value", lambda_grid = 1:5), regexp = "Unknown")

  ps <- presidential_speech
  ps[1,1] <- Inf; expect_error(convex_biclustering(ps, lambda_grid = 1:5))

  ## lambda_grid must be supplied, strictly positive, and ordered (only gets a warning)
  expect_error(convex_biclustering(presidential_speech))
  expect_error(convex_biclustering(presidential_speech), lambda_grid = numeric())
  expect_error(convex_biclustering(presidential_speech, lambda_grid = 0))
  expect_error(convex_biclustering(presidential_speech, lambda_grid = c(0, 3)))
  expect_warning(convex_biclustering(presidential_speech, lambda_grid = c(3, 2, 1)))
})

test_that("convex_biclustering() matches cvxclustr", {
  skip_on_cran()
  skip_if_not_installed("cvxclustr")
  skip_if_not_installed("cvxbiclustr")
  skip_if_not_installed("Matrix")
  data(mammals, package = "cvxclustr")
  library(Matrix)

  ## Example modified from help pages of cvxclustr
  X <- as.matrix(mammals[,-1])
  biclust_fit <- convex_biclustering(X, X.center.global = FALSE, lambda_grid = seq(4, 40, length.out = 10))

  ## Calculate matching `cvxbiclustr` solution

  ## Match CBASS() selected edges and weights
  D_row <- biclust_fit$D_row
  D_col <- t(biclust_fit$D_col)

  w_row <- biclust_fit$row_weights; w_row <- clustRviz:::weight_mat_to_vec(w_row); w_row <- w_row[w_row != 0]
  w_col <- biclust_fit$col_weights; w_col <- clustRviz:::weight_mat_to_vec(w_col); w_col <- w_col[w_col != 0]

  w_row <- w_row / sum(w_row); w_row <- w_row / sqrt(NROW(X))
  w_col <- w_col / sum(w_col); w_col <- w_col / sqrt(NCOL(X))

  ## Perform clustering
  capture.output(cobra_fit <- cvxbiclustr::cobra(X,
                                                 E_row = Matrix(D_row, sparse = TRUE),
                                                 E_col = Matrix(D_col, sparse = TRUE),
                                                 w_row = w_row,
                                                 w_col = w_col,
                                                 gamma = biclust_fit$lambda_grid))

  obj <- function(U, lambda, l1 = FALSE){
    DU <- D_row %*% U; DTU <- D_col %*% t(U)
    if(l1){
      sum((X - U)^2)/2 + lambda * sum(w_row * abs(DU)) + lambda * sum(w_col * abs(t(DTU)))
    } else {
      sum((X - U)^2)/2 + lambda * sum(w_row * sqrt(rowSums(DU^2))) + lambda * sum(w_col * sqrt(rowSums(DTU^2)))
    }
  }

  ## Check that the objectives are of (roughly) equal quality...
  for(ix in seq_along(biclust_fit$lambda_grid)){
    lambda <- biclust_fit$lambda_grid[ix]

    my_U <- biclust_fit$U[,,ix]
    ec_U <- cobra_fit$U[[ix]]

    expect_true(obj(my_U, lambda) <= 1.001 * obj(ec_U, lambda))
  }
})
