context("Test CARP Sparsity Invariants")

test_that("CARP begins with no fusions", {
  clustRviz_options(keep_debug_info = TRUE)

  carp_fit <- CARP(presidential_speech, alg.type = 'carp')
  expect_zeros(carp_fit$debug$v_zero_indices[, 1])

  clustRviz_reset_options()
})

test_that("CARP ends with all fusions", {
  clustRviz_options(keep_debug_info = TRUE)

  carp_fit <- CARP(presidential_speech, alg.type = 'carp')
  expect_ones(carp_fit$debug$v_zero_indices[, NCOL(carp_fit$debug$v_zero_indices)])

  clustRviz_reset_options()
})

test_that("CARP-VIZ begins with no fusions", {
  clustRviz_options(keep_debug_info = TRUE)

  carp_fit <- CARP(presidential_speech, alg.type = 'carpviz')
  expect_zeros(carp_fit$debug$v_zero_indices[, 1])

  clustRviz_reset_options()
})

test_that("CARP-VIZ ends with all fusions", {
  clustRviz_options(keep_debug_info = TRUE)

  carp_fit <- CARP(presidential_speech,alg.type='carpviz')
  expect_ones(carp_fit$debug$v_zero_indices[, NCOL(carp_fit$debug$v_zero_indices)])

  clustRviz_reset_options()
})

test_that("CARP begins with no fusions (uniform weights)", {
  clustRviz_options(keep_debug_info = TRUE)

  weight_mat <- matrix(1, nrow=NROW(presidential_speech), ncol=NROW(presidential_speech))
  carp_fit   <- CARP(presidential_speech, weights = weight_mat, alg.type = 'carp')
  expect_zeros(carp_fit$debug$v_zero_indices[, 1])

  clustRviz_reset_options()
})

test_that("CARP ends with all fusions (uniform weights)", {
  clustRviz_options(keep_debug_info = TRUE)

  weight_mat <- matrix(1, nrow=NROW(presidential_speech), ncol=NROW(presidential_speech))
  carp_fit   <- CARP(presidential_speech, weights = weight_mat, alg.type = 'carp')
  expect_ones(carp_fit$debug$v_zero_indices[, NCOL(carp_fit$debug$v_zero_indices)])

  clustRviz_reset_options()
})

test_that("CARP-VIZ begins with no fusions (uniform weights)", {
  clustRviz_options(keep_debug_info = TRUE)

  weight_mat <- matrix(1, nrow=NROW(presidential_speech), ncol=NROW(presidential_speech))
  carp_fit   <- CARP(presidential_speech, weights = weight_mat, alg.type = 'carpviz')
  expect_zeros(carp_fit$debug$v_zero_indices[, 1])

  clustRviz_reset_options()
})

test_that("CARP-VIZ ends with all fusions (uniform weights)", {
  clustRviz_options(keep_debug_info = TRUE)

  weight_mat <- matrix(1, nrow=NROW(presidential_speech), ncol=NROW(presidential_speech))
  carp_fit   <- CARP(presidential_speech, weights = weight_mat, alg.type = 'carpviz')
  expect_ones(carp_fit$debug$v_zero_indices[, NCOL(carp_fit$debug$v_zero_indices)])

  clustRviz_reset_options()
})
