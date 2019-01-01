context("Test CBASS Sparsity Invariants")

test_that("CBASS begins with no fusions", {
  clustRviz_options(keep_debug_info = TRUE)

  cbass_fit <- CBASS(presidential_speech, back_track = FALSE)
  expect_zeros(cbass_fit$debug$row$v_zero_indices[, 1])
  expect_zeros(cbass_fit$debug$col$v_zero_indices[, 1])

  clustRviz_reset_options()
})

test_that("CBASS ends with full fusions", {
  clustRviz_options(keep_debug_info = TRUE)

  cbass_fit <- CBASS(presidential_speech, back_track = FALSE)
  expect_ones(cbass_fit$debug$row$v_zero_indices[, NCOL(cbass_fit$debug$row$v_zero_indices)])
  expect_ones(cbass_fit$debug$col$v_zero_indices[, NCOL(cbass_fit$debug$col$v_zero_indices)])

  clustRviz_reset_options()
})

test_that("CBASS-VIZ begins with no fusions", {
  clustRviz_options(keep_debug_info = TRUE)

  cbass_fit <- CBASS(presidential_speech, back_track = TRUE)
  expect_zeros(cbass_fit$debug$row$v_zero_indices[, 1])
  expect_zeros(cbass_fit$debug$col$v_zero_indices[, 1])

  clustRviz_reset_options()
})

test_that("CBASS-VIZ ends with full fusions", {
  clustRviz_options(keep_debug_info = TRUE)

  cbass_fit <- CBASS(presidential_speech, back_track = TRUE)
  expect_ones(cbass_fit$debug$row$v_zero_indices[, NCOL(cbass_fit$debug$row$v_zero_indices)])
  expect_ones(cbass_fit$debug$col$v_zero_indices[, NCOL(cbass_fit$debug$col$v_zero_indices)])

  clustRviz_reset_options()
})

test_that("CBASS begins with no fusions (uniform weights)", {
  clustRviz_options(keep_debug_info = TRUE)

  uniform_weights <- function(X) matrix(1, nrow = NROW(X), ncol = NROW(X))
  cbass_fit <- CBASS(presidential_speech, back_track = FALSE, row_weights = uniform_weights, col_weights = uniform_weights)
  expect_zeros(cbass_fit$debug$row$v_zero_indices[, 1])
  expect_zeros(cbass_fit$debug$col$v_zero_indices[, 1])

  clustRviz_reset_options()
})

test_that("CBASS ends with full fusions (uniform weights)", {
  clustRviz_options(keep_debug_info = TRUE)

  uniform_weights <- function(X) matrix(1, nrow = NROW(X), ncol = NROW(X))
  cbass_fit <- CBASS(presidential_speech, back_track = FALSE, row_weights = uniform_weights, col_weights = uniform_weights)
  expect_ones(cbass_fit$debug$row$v_zero_indices[, NCOL(cbass_fit$debug$row$v_zero_indices)])
  expect_ones(cbass_fit$debug$col$v_zero_indices[, NCOL(cbass_fit$debug$col$v_zero_indices)])

  clustRviz_reset_options()
})

test_that("CBASS-VIZ begins with no fusions (uniform weights)", {
  clustRviz_options(keep_debug_info = TRUE)

  uniform_weights <- function(X) matrix(1, nrow = NROW(X), ncol = NROW(X))
  cbass_fit <- CBASS(presidential_speech, back_track = TRUE, row_weights = uniform_weights, col_weights = uniform_weights)
  expect_zeros(cbass_fit$debug$row$v_zero_indices[, 1])
  expect_zeros(cbass_fit$debug$col$v_zero_indices[, 1])

  clustRviz_reset_options()
})

test_that("CBASS-VIZ ends with full fusions (uniform weights)", {
  clustRviz_options(keep_debug_info = TRUE)

  uniform_weights <- function(X) matrix(1, nrow = NROW(X), ncol = NROW(X))
  cbass_fit <- CBASS(presidential_speech, back_track = TRUE, row_weights = uniform_weights, col_weights = uniform_weights)
  expect_ones(cbass_fit$debug$row$v_zero_indices[, NCOL(cbass_fit$debug$row$v_zero_indices)])
  expect_ones(cbass_fit$debug$col$v_zero_indices[, NCOL(cbass_fit$debug$col$v_zero_indices)])

  clustRviz_reset_options()
})
