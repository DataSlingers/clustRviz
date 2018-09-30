context("Test CBASS Sparsity Invariants")

test_that("CBASS begins with no fusions", {
  cbass_fit <- CBASS(presidential_speech, alg.type = 'cbass')
  expect_zeros(cbass_fit$cbass.sol.path$v.col.zero.inds[, 1])
  expect_zeros(cbass_fit$cbass.sol.path$v.row.zero.inds[, 1])
})

test_that("CBASS ends with full fusions", {
  cbass_fit <- CBASS(presidential_speech, alg.type = 'cbass')
  expect_ones(cbass_fit$cbass.sol.path$v.col.zero.inds[, NCOL(cbass_fit$cbass.sol.path$v.col.zero.inds)])
  expect_ones(cbass_fit$cbass.sol.path$v.row.zero.inds[, NCOL(cbass_fit$cbass.sol.path$v.row.zero.inds)])
})

test_that("CBASS-VIZ begins with no fusions", {
  cbass_fit <- CBASS(presidential_speech, alg.type = 'cbassviz')
  expect_zeros(cbass_fit$cbass.sol.path$v.col.zero.inds[, 1])
  expect_zeros(cbass_fit$cbass.sol.path$v.row.zero.inds[, 1])
})

test_that("CBASS-VIZ ends with full fusions", {
  cbass_fit <- CBASS(presidential_speech, alg.type = 'cbassviz')
  expect_ones(cbass_fit$cbass.sol.path$v.col.zero.inds[, NCOL(cbass_fit$cbass.sol.path$v.col.zero.inds)])
  expect_ones(cbass_fit$cbass.sol.path$v.row.zero.inds[, NCOL(cbass_fit$cbass.sol.path$v.row.zero.inds)])
})

test_that("CBASS begins with no fusions (uniform weights)", {
  uniform_weights <- function(X) matrix(1, nrow = NROW(X), ncol = NROW(X))
  cbass_fit <- CBASS(presidential_speech, alg.type = 'cbass', row_weights = uniform_weights, col_weights = uniform_weights)
  expect_zeros(cbass_fit$cbass.sol.path$v.col.zero.inds[, 1])
  expect_zeros(cbass_fit$cbass.sol.path$v.row.zero.inds[, 1])
})

test_that("CBASS ends with full fusions (uniform weights)", {
  uniform_weights <- function(X) matrix(1, nrow = NROW(X), ncol = NROW(X))
  cbass_fit <- CBASS(presidential_speech, alg.type = 'cbass', row_weights = uniform_weights, col_weights = uniform_weights)
  expect_ones(cbass_fit$cbass.sol.path$v.col.zero.inds[, NCOL(cbass_fit$cbass.sol.path$v.col.zero.inds)])
  expect_ones(cbass_fit$cbass.sol.path$v.row.zero.inds[, NCOL(cbass_fit$cbass.sol.path$v.row.zero.inds)])
})

test_that("CBASS-VIZ begins with no fusions (uniform weights)", {
  uniform_weights <- function(X) matrix(1, nrow = NROW(X), ncol = NROW(X))
  cbass_fit <- CBASS(presidential_speech, alg.type = 'cbassviz', row_weights = uniform_weights, col_weights = uniform_weights)
  expect_zeros(cbass_fit$cbass.sol.path$v.col.zero.inds[, 1])
  expect_zeros(cbass_fit$cbass.sol.path$v.row.zero.inds[, 1])
})

test_that("CBASS-VIZ ends with full fusions (uniform weights)", {
  uniform_weights <- function(X) matrix(1, nrow = NROW(X), ncol = NROW(X))
  cbass_fit <- CBASS(presidential_speech, alg.type = 'cbassviz', row_weights = uniform_weights, col_weights = uniform_weights)
  expect_ones(cbass_fit$cbass.sol.path$v.col.zero.inds[, NCOL(cbass_fit$cbass.sol.path$v.col.zero.inds)])
  expect_ones(cbass_fit$cbass.sol.path$v.row.zero.inds[, NCOL(cbass_fit$cbass.sol.path$v.row.zero.inds)])
})
