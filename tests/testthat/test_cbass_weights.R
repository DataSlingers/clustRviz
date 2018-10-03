context("Test CBASS Weight Handling")

## Tests for row weights
test_that("CBASS works with user weight function for row weights", {
  # FIXME: This should work - see GitHub #8
  ## uniform_weight_func <- function(X) matrix(1, nrow=NROW(X), ncol=NROW(X))
  ## CBASS(presidential_speech, row_weights = uniform_weight_func)

  # By manual testing, this is a good weight function / matrix
  my_weight_func <- function(X) {
    mat <- exp(-0.01 * as.matrix(dist(X))^2)
    diag(mat) <- 0
    mat[mat < quantile(mat, 0.73)] <- 0
    mat
  }

  expect_true(clustRviz:::is_connected_adj_mat(my_weight_func(presidential_speech)))
  expect_no_error(CBASS(presidential_speech, row_weights = my_weight_func))
})

test_that("CBASS works with user weight matrix for row weights", {
  mat <- exp(-0.01 * as.matrix(dist(presidential_speech))^2)
  diag(mat) <- 0
  mat[mat < quantile(mat, 0.73)] <- 0
  expect_no_error(CBASS(presidential_speech, row_weights = mat))
})

test_that("CBASS errors with negative weights for row weights", {
  mat <- -1 * exp(-0.01 * as.matrix(dist(presidential_speech))^2)
  diag(mat) <- 0
  mat[mat < quantile(mat, 0.73)] <- 0
  expect_error(CBASS(presidential_speech, row_weights = mat))
})

test_that("CBASS errors with unconnected graphs for row weights", {
  mat <- exp(-0.01 * as.matrix(dist(presidential_speech))^2)
  diag(mat) <- 0
  mat[mat < quantile(mat, 0.95)] <- 0
  expect_error(CBASS(presidential_speech, row_weights = mat))
})

## Tests for column / feature weights
test_that("CBASS works with user weight function for column weights", {
  # FIXME: This should work - see GitHub #8
  ## uniform_weight_func <- function(X) matrix(1, nrow=NROW(X), ncol=NROW(X))
  ## CBASS(presidential_speech, col_weights = uniform_weight_func)

  # By manual testing, this is a good weight function / matrix
  my_weight_func <- function(X) {
    mat <- exp(-0.01 * as.matrix(dist(X))^2)
    diag(mat) <- 0
    mat[mat < quantile(mat, 0.73)] <- 0
    mat
  }

  expect_true(clustRviz:::is_connected_adj_mat(my_weight_func(t(presidential_speech))))
  expect_no_error(CBASS(presidential_speech, col_weights = my_weight_func))
})

test_that("CBASS works with user weight matrix for column weights", {
  mat <- exp(-0.01 * as.matrix(dist(t(presidential_speech)))^2)
  diag(mat) <- 0
  mat[mat < quantile(mat, 0.73)] <- 0
  expect_no_error(CBASS(presidential_speech, col_weights = mat))
})

test_that("CBASS errors with negative weights for column weights", {
  mat <- -1 * exp(-0.01 * as.matrix(dist(t(presidential_speech)))^2)
  diag(mat) <- 0
  mat[mat < quantile(mat, 0.73)] <- 0
  expect_error(CBASS(presidential_speech, col_weights = mat))
})

test_that("CBASS errors with unconnected graphs for column weights", {
  mat <- exp(-0.01 * as.matrix(dist(t(presidential_speech)))^2)
  diag(mat) <- 0
  mat[mat < quantile(mat, 0.95)] <- 0
  expect_error(CBASS(presidential_speech, col_weights = mat))
})
