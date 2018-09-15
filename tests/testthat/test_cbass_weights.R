context("Test CBASS Weight Handling")

## Tests for observation weights
test_that("CBASS works with user weight function for observation weights", {
  # FIXME: This should work - see GitHub #8
  ## uniform_weight_func <- function(X) matrix(1, nrow=NROW(X), ncol=NROW(X))
  ## CBASS(presidential_speech, obs_weights = uniform_weight_func)

  # By manual testing, this is a good weight function / matrix
  my_weight_func <- function(X) {
    mat <- exp(-0.01 * as.matrix(dist(X))^2)
    diag(mat) <- 0
    mat[mat < quantile(mat, 0.73)] <- 0
    mat
  }

  expect_true(clustRviz:::is_connected_adj_mat(my_weight_func(presidential_speech)))
  expect_no_error(CBASS(presidential_speech, obs_weights = my_weight_func))
})

test_that("CBASS works with user weight matrix for observation weights", {
  mat <- exp(-0.01 * as.matrix(dist(presidential_speech))^2)
  diag(mat) <- 0
  mat[mat < quantile(mat, 0.73)] <- 0
  expect_no_error(CBASS(presidential_speech, obs_weights = mat))
})

test_that("CBASS errors with negative weights for observation weights", {
  mat <- -1 * exp(-0.01 * as.matrix(dist(presidential_speech))^2)
  diag(mat) <- 0
  mat[mat < quantile(mat, 0.73)] <- 0
  expect_error(CBASS(presidential_speech, obs_weights = mat))
})

test_that("CBASS errors with unconnected graphs for observation weights", {
  mat <- exp(-0.01 * as.matrix(dist(presidential_speech))^2)
  diag(mat) <- 0
  mat[mat < quantile(mat, 0.95)] <- 0
  expect_error(CBASS(presidential_speech, obs_weights = mat))
})

## Tests for variable / feature weights
test_that("CBASS works with user weight function for variable weights", {
  # FIXME: This should work - see GitHub #8
  ## uniform_weight_func <- function(X) matrix(1, nrow=NROW(X), ncol=NROW(X))
  ## CBASS(presidential_speech, var_weights = uniform_weight_func)

  # By manual testing, this is a good weight function / matrix
  my_weight_func <- function(X) {
    mat <- exp(-0.01 * as.matrix(dist(X))^2)
    diag(mat) <- 0
    mat[mat < quantile(mat, 0.73)] <- 0
    mat
  }

  expect_true(clustRviz:::is_connected_adj_mat(my_weight_func(t(presidential_speech))))
  expect_no_error(CBASS(presidential_speech, var_weights = my_weight_func))
})

test_that("CBASS works with user weight matrix for variable weights", {
  mat <- exp(-0.01 * as.matrix(dist(t(presidential_speech)))^2)
  diag(mat) <- 0
  mat[mat < quantile(mat, 0.73)] <- 0
  expect_no_error(CBASS(presidential_speech, var_weights = mat))
})

test_that("CBASS errors with negative weights for variable weights", {
  mat <- -1 * exp(-0.01 * as.matrix(dist(t(presidential_speech)))^2)
  diag(mat) <- 0
  mat[mat < quantile(mat, 0.73)] <- 0
  expect_error(CBASS(presidential_speech, var_weights = mat))
})

test_that("CBASS errors with unconnected graphs for variable weights", {
  mat <- exp(-0.01 * as.matrix(dist(t(presidential_speech)))^2)
  diag(mat) <- 0
  mat[mat < quantile(mat, 0.95)] <- 0
  expect_error(CBASS(presidential_speech, var_weights = mat))
})
