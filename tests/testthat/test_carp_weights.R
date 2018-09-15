context("Test CARP Weight Handling")

test_that("CARP works with user weight function", {
    # FIXME: This should work - see GitHub #8
    ## uniform_weight_func <- function(X) matrix(1, nrow=NROW(X), ncol=NROW(X))
    ## CARP(presidential_speech, weights = uniform_weight_func)

    # By manual testing, this is a good weight function / matrix
    my_weight_func <- function(X) {
      mat <- exp(-0.01 * as.matrix(dist(X))^2)
      diag(mat) <- 0
      mat[mat < quantile(mat, 0.73)] <- 0
      mat
    }

    expect_true(clustRviz:::is_connected_adj_mat(my_weight_func(presidential_speech)))
    expect_no_error(CARP(presidential_speech, weights = my_weight_func))
})

test_that("CARP works with user weight matrix", {
    mat <- exp(-0.01 * as.matrix(dist(presidential_speech))^2)
    diag(mat) <- 0
    mat[mat < quantile(mat, 0.73)] <- 0
    expect_no_error(CARP(presidential_speech, weights = mat))
})

test_that("CARP errors with negative weights", {
  mat <- -1 * exp(-0.01 * as.matrix(dist(presidential_speech))^2)
  diag(mat) <- 0
  mat[mat < quantile(mat, 0.73)] <- 0
  expect_error(CARP(presidential_speech, weights = mat))
})

test_that("CARP errors with unconnected graphs", {
  mat <- exp(-0.01 * as.matrix(dist(presidential_speech))^2)
  diag(mat) <- 0
  mat[mat < quantile(mat, 0.95)] <- 0
  expect_error(CARP(presidential_speech, weights = mat))
})
