context("Test tensor PCA projection")

test_that("tensor_projection works", {
  tensor_projection <- clustRviz:::tensor_projection
  n <- 50
  p <- 30
  q <- 100

  X <- array(rnorm(n * p * q), dim = c(n, p, q))

  k <- 10
  P <- matrix(rnorm(p * k), nrow = p, ncol = k)

  XP <- array(NA, dim = c(n, k, q))
  for(i in seq(1, q)){
    XP[,,i] <- X[,,i] %*% P
  }

  expect_equal(tensor_projection(X, P), XP)
})
