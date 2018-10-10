context("Test Built-In Weight Functions")

test_that("Dense RBF works with fixed phi", {
    weight_func <- dense_rbf_kernel_weights(phi = 1)
    weight_results <- weight_func(presidential_speech)

    weight_mat_manual <- exp(-as.matrix(dist(presidential_speech))^2)
    diag(weight_mat_manual) <- 0

    expect_equal(weight_results$weight_mat, weight_mat_manual)
})

test_that("Dense RBF works with learned phi", {
  weight_func <- dense_rbf_kernel_weights()
  weight_results <- weight_func(presidential_speech)
  phi <- weight_results$type$phi

  weight_mat_manual <- exp(-phi * as.matrix(dist(presidential_speech))^2)
  diag(weight_mat_manual) <- 0

  expect_equal(weight_results$weight_mat, weight_mat_manual)
})

test_that("Sparse RBF with full k is a no-op", {
  sparse_weight_func <- sparse_rbf_kernel_weights(k = NROW(presidential_speech) - 1)
  dense_weight_func  <- dense_rbf_kernel_weights()

  expect_equal(sparse_weight_func(presidential_speech)$weight_mat,
               dense_weight_func(presidential_speech)$weight_mat)
})

test_that("Sparse RBF with learned k is same as if k were known a priori", {
  weight_func    <- sparse_rbf_kernel_weights()
  weight_results <- weight_func(presidential_speech)
  k <- weight_results$type$k

  weight_func2    <- sparse_rbf_kernel_weights(k = k)
  weight_results2 <- weight_func2(presidential_speech)

  expect_equal(weight_results$weight_mat,
               weight_results2$weight_mat)
})

test_that("Dense RBF works with Manhattan distance", {
  weight_func <- dense_rbf_kernel_weights(dist.method = "manhattan")
  weight_results <- weight_func(presidential_speech)
  phi <- weight_results$type$phi

  weight_mat_manual <- exp(-phi * as.matrix(dist(presidential_speech, method = "manhattan"))^2)
  diag(weight_mat_manual) <- 0

  expect_equal(weight_results$weight_mat, weight_mat_manual)
})

test_that("Dense RBF checks inputs", {
  expect_error(dense_rbf_kernel_weights(dist.method = "London"))
  expect_error(dense_rbf_kernel_weights(dist.method = NA))
  expect_error(dense_rbf_kernel_weights(p = 0))
  expect_error(dense_rbf_kernel_weights(p = -3))
  expect_error(dense_rbf_kernel_weights(p = NA))
  expect_error(dense_rbf_kernel_weights(p = c(1, 2)))
  expect_error(dense_rbf_kernel_weights(phi = 0)(presidential_speech))
})

test_that("Sparse RBF checks inputs", {
  expect_error(sparse_rbf_kernel_weights(dist.method = "London"))
  expect_error(sparse_rbf_kernel_weights(dist.method = NA))
  expect_error(sparse_rbf_kernel_weights(p = 0))
  expect_error(sparse_rbf_kernel_weights(p = -3))
  expect_error(sparse_rbf_kernel_weights(p = NA))
  expect_error(sparse_rbf_kernel_weights(p = c(1, 2)))
  expect_error(sparse_rbf_kernel_weights(phi = 0)(presidential_speech))
  expect_error(sparse_rbf_kernel_weights(k = 0)(presidential_speech))
  expect_error(sparse_rbf_kernel_weights(k = -1)(presidential_speech))
})

test_that("Print method works - Dense RBF", {
  weight_func <- dense_rbf_kernel_weights(phi = 1)
  weight_fit_obj <- weight_func(presidential_speech)$type
  weight_print <- capture_print(weight_fit_obj)

  expect_str_contains(weight_print, "Radial Basis Function Kernel Weights")
  expect_str_contains(weight_print, "Distance Metric: Euclidean")
  expect_str_contains(weight_print, stringr::fixed("Scale parameter (phi): 1 [User-Supplied]"))


  weight_func <- dense_rbf_kernel_weights()
  weight_fit_obj <- weight_func(presidential_speech)$type
  weight_print <- capture_print(weight_fit_obj)

  expect_str_contains(weight_print, "Radial Basis Function Kernel Weights")
  expect_str_contains(weight_print, "Distance Metric: Euclidean")
  expect_str_contains(weight_print, stringr::fixed("Scale parameter (phi): 0.01 [Data-Driven]"))
})

test_that("Print method works - Sparse RBF", {
  weight_func <- sparse_rbf_kernel_weights(k = 10)
  weight_fit_obj <- weight_func(presidential_speech)$type
  weight_print <- capture_print(weight_fit_obj)
  expect_str_contains(weight_print, stringr::fixed("Sparsified: 10 Nearest Neighbors [User-Supplied]"))

  weight_func <- sparse_rbf_kernel_weights()
  weight_fit_obj <- weight_func(presidential_speech)$type
  weight_print <- capture_print(weight_fit_obj)
  expect_str_contains(weight_print, stringr::fixed("Sparsified: 4 Nearest Neighbors [Data-Driven]"))
})

test_that("Print method works - User-Function", {
  expect_str_contains(capture_print(clustRviz:::UserFunction()),
                      "Source: User-Provided Function")
})


test_that("Print method works - User-Matrix", {
  expect_str_contains(capture_print(clustRviz:::UserMatrix()),
                      "Source: User-Provided Matrix")
})

test_that("Dense weights connectedness check works", {
  check_weight_matrix <- clustRviz:::check_weight_matrix
  eye <- function(n) diag(1, nrow = n, ncol = n)

  W <- eye(3)
  expect_error(check_weight_matrix(W))

  W <- matrix(1, 3, 3)
  expect_no_error(check_weight_matrix(W))

  W <- matrix(1, 3, 3); W[1,2] <- W[2,1] <- W[1,3] <- W[3,1] <- 0
  ## No check for first ob since we are only checking lower triangle
  expect_error(check_weight_matrix(W), regexp = "observation 2")

  W <- matrix(1, 3, 4)
  expect_error(check_weight_matrix(W), regexp = "square")
})
