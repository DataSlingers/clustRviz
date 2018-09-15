context("test clustering.CARP")

test_that("No additional arguments gives the whole path", {
  carp_fit <- CARP(presidential_speech)

  clustering_result <- clustering(carp_fit)

  ## The number of unique elements in clustering.assignment should increase uniformly
  expect_equal(apply(clustering_result$clustering.assignment, 1, num_unique),
               1:NROW(presidential_speech))

  ## The number of unique means in cluster.means should increase uniformly
  expect_equal(vapply(clustering_result$cluster.means, NCOL, integer(1)),
               1:NROW(presidential_speech))
})

test_that("Works with user `k`", {
  carp_fit <- CARP(presidential_speech)

  clustering_result <- clustering(carp_fit, k = 5)

  expect_equal(length(table(clustering_result$clustering.assignment)), 5)
  expect_equal(NCOL(clustering_result$cluster.means), 5)

  clustering_result <- clustering(carp_fit, k = NROW(presidential_speech))
  expect_equal(clustering_result$cluster.means,
               t(presidential_speech),
               check.attributes=FALSE)

  clustering_result <- clustering(carp_fit, k = 1)
  expect_equal(clustering_result$cluster.means,
               matrix(colMeans(presidential_speech)),
               check.attributes=FALSE)

  expect_error(clustering(carp_fit, k = 0))
  expect_error(clustering(carp_fit, k = NROW(presidential_speech) + 1))
})

test_that("Works with user `percent`", {
  carp_fit <- CARP(presidential_speech)

  clustering_result <- clustering(carp_fit, percent = 0.5)

  n_clusters <- length(table(clustering_result$clustering.assignment))
  expect_equal(NCOL(clustering_result$cluster.means), n_clusters)

  clustering_result <- clustering(carp_fit, percent = 0)
  expect_equal(clustering_result$cluster.means,
               t(presidential_speech),
               check.attributes=FALSE)

  clustering_result <- clustering(carp_fit, percent = 1)
  expect_equal(clustering_result$cluster.means,
               matrix(colMeans(presidential_speech)),
               check.attributes=FALSE)

  expect_error(clustering(carp_fit, percent = +1.2))
  expect_error(clustering(carp_fit, percent = -0.2))
})

test_that("Supplying percent and k gives error", {
  carp_fit <- CARP(presidential_speech)

  expect_error(clustering(carp_fit, k = 5, percent = 0.25))
})
