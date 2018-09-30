context("Test CARP Sparsity Invariants")

test_that("CARP begins with no fusions", {
  carp_fit <- CARP(presidential_speech, alg.type = 'carp')
  expect_zeros(carp_fit$carp.sol.path$v.zero.inds[, 1])
})

test_that("CARP ends with all fusions", {
  carp_fit <- CARP(presidential_speech, alg.type = 'carp')
  expect_ones(carp_fit$carp.sol.path$v.zero.inds[, NCOL(carp_fit$carp.sol.path$v.zero.inds)])
})

test_that("CARP-VIZ begins with no fusions", {
  carp_fit <- CARP(presidential_speech, alg.type = 'carpviz')
  expect_zeros(carp_fit$carp.sol.path$v.zero.inds[, 1])
})

test_that("CARP-VIZ ends with all fusions", {
  carp_fit <- CARP(presidential_speech,alg.type='carpviz')
  expect_ones(carp_fit$carp.sol.path$v.zero.inds[, NCOL(carp_fit$carp.sol.path$v.zero.inds)])
})

test_that("CARP begins with no fusions (uniform weights)", {
  weight_mat <- matrix(1, nrow=NROW(presidential_speech), ncol=NROW(presidential_speech))
  carp_fit   <- CARP(presidential_speech, weights = weight_mat, alg.type = 'carp')
  expect_zeros(carp_fit$carp.sol.path$v.zero.inds[, 1])
})

test_that("CARP ends with all fusions (uniform weights)", {
  weight_mat <- matrix(1, nrow=NROW(presidential_speech), ncol=NROW(presidential_speech))
  carp_fit   <- CARP(presidential_speech, weights = weight_mat, alg.type = 'carp')
  expect_ones(carp_fit$carp.sol.path$v.zero.inds[, NCOL(carp_fit$carp.sol.path$v.zero.inds)])
})

test_that("CARP-VIZ begins with no fusions (uniform weights)", {
  weight_mat <- matrix(1, nrow=NROW(presidential_speech), ncol=NROW(presidential_speech))
  carp_fit   <- CARP(presidential_speech, weights = weight_mat, alg.type = 'carpviz')
  expect_zeros(carp_fit$carp.sol.path$v.zero.inds[, 1])
})

test_that("CARP-VIZ ends with all fusions (uniform weights)", {
  weight_mat <- matrix(1, nrow=NROW(presidential_speech), ncol=NROW(presidential_speech))
  carp_fit   <- CARP(presidential_speech, weights = weight_mat, alg.type = 'carpviz')
  expect_ones(carp_fit$carp.sol.path$v.zero.inds[, NCOL(carp_fit$carp.sol.path$v.zero.inds)])
})
