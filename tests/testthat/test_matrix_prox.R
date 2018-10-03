context("Test C++ matrix prox")

test_that("L1 matrix prox works", {
  set.seed(125)
  n <- 25
  p <- 50

  X <- matrix(rnorm(n * p), nrow = n, ncol = p)

  MatrixProx <- clustRviz:::MatrixProx
  weights <- rep(1, n)

  expect_equal(X, MatrixProx(X, lambda = 0, weights = weights, l1 = TRUE))
  expect_equal(abs(X) + 4, MatrixProx(abs(X) + 5, lambda = 1, weights = weights, l1 = TRUE))
  expect_equal(-abs(X) - 4, MatrixProx(-abs(X) - 5, lambda = 1, weights = weights, l1 = TRUE))

  ## Now we check that weights work
  X <- matrix(1:25, nrow = 25, ncol = 1)
  weights <- 1:25
  expect_equal(matrix(0, nrow = 25, ncol = 1),
               MatrixProx(X, lambda = 1, weights = weights, l1 = TRUE))

  X <- matrix(5, nrow = 6, ncol = 1)
  weights <- seq(0, 5)
  expect_equal(matrix(5 - weights, nrow = 6, ncol = 1),
               MatrixProx(X, lambda = 1, weights = weights, l1 = TRUE))
})

test_that("L2 prox works", {
  set.seed(125)
  MatrixProx <- clustRviz:::MatrixProx
  num_unique_cols <- clustRviz:::num_unique_cols
  n <- 25

  ## If X has a single column, same as L1 prox
  X <- matrix(rnorm(n, sd = 3), ncol = 1)
  weights <- rexp(n)

  expect_equal(MatrixProx(X, lambda = 1, weights = weights, l1 = TRUE),
               MatrixProx(X, lambda = 1, weights = weights, l1 = FALSE))

  p <- 5
  X <- matrix(1, nrow = n, ncol = p)
  weights <- seq(0, 5, length.out = 25)

  expect_equal(1, num_unique_cols(MatrixProx(X, lambda = 1, weights = weights, l1 = FALSE)))

  y <- matrix(c(3, 4), nrow = 1)

  expect_equal(MatrixProx(y, 1, 1, l1 = FALSE), y * (1 - 1/5))
  expect_equal(MatrixProx(y, 1, 3, l1 = FALSE), y * (1 - 3/5))
  expect_equal(MatrixProx(y, 2, 1, l1 = FALSE), y * (1 - 2/5))
  expect_equal(MatrixProx(y, 2, 3, l1 = FALSE), y * 0)

  y <- -1 * y
  expect_equal(MatrixProx(y, 1, 1, l1 = FALSE), y * (1 - 1/5))
  expect_equal(MatrixProx(y, 1, 3, l1 = FALSE), y * (1 - 3/5))
  expect_equal(MatrixProx(y, 2, 1, l1 = FALSE), y * (1 - 2/5))
  expect_equal(MatrixProx(y, 2, 3, l1 = FALSE), y * 0)
})
