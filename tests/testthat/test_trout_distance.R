context("Test trout_distance")

test_that("Trout distance doesn't depend on signs", {
  trout_dist <- clustRviz:::trout_dist

  X <- matrix(c(1, 1, -1, -1), ncol = 2)
  expect_equal(as.vector(trout_dist(X)), 0)

  X <- matrix(c(1, -1, 1, -1), ncol = 2)
  expect_equal(as.vector(trout_dist(X)), 0)

  X <- matrix(c(1 + 1i, -1 -1i, 1 + 1i, -1 -1i), ncol = 2)
  expect_equal(as.vector(trout_dist(X)), 0)
})

test_that("Trout distance minimizes distance", {
  trout_dist <- clustRviz:::trout_dist
  set.seed(125)

  X1 <- rnorm(25) + (0 + 1i) * rnorm(25)
  X2 <- rnorm(25) + (0 + 1i) * rnorm(25)
  X <- rbind(X1, X2)

  theta_grid <- seq(0, 2 * pi, length.out = 501)
  d <- Vectorize(function(theta) sum(Mod(X1 - exp((0 + 1i) * theta) * X2)^2))
  min_d <- sqrt(min(d(theta_grid)))

  td <- as.vector(trout_dist(X))
  expect_lte(td, min_d)
})
