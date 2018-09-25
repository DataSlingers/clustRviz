context("test clustering.CARP")


test_that("CARP-returned sparsity pattern begins with no fusions",{
  carp.fit <- CARP(presidential_speech, alg.type = 'carp')
  expect_equal(
    sum(carp.fit$carp.sol.path$v.zero.inds[, 1]),
    0
  )
})
test_that("CARP-returned sparsity pattern ends with full fusions",{
  carp.fit <- CARP(presidential_speech, alg.type = 'carp')
  expect_true(
    all(carp.fit$carp.sol.path$v.zero.inds[, ncol(carp.fit$carp.sol.path$v.zero.inds)] == 1)
  )
})

test_that("CARPVIZ-returned sparsity pattern begins with no fusions",{
  carp.fit <- CARP(presidential_speech, alg.type = 'carpviz')
  expect_equal(
    sum(carp.fit$carp.sol.path$v.zero.inds[, 1]),
    0
  )
})
test_that("CARPVIZ-returned sparsity pattern ends with full fusions",{
  carp.fit <- CARP(presidential_speech,alg.type='carpviz')
  expect_true(
    all(carp.fit$carp.sol.path$v.zero.inds[, ncol(carp.fit$carp.sol.path$v.zero.inds)] == 1)
  )
})





test_that("CARP-returned sparsity pattern begins with no fusions (uniform weights)",{
  weight_mat <-matrix(1, nrow=nrow(presidential_speech), ncol=nrow(presidential_speech))
  carp.fit <- CARP(presidential_speech, weights = weight_mat, alg.type = 'carp')
  expect_equal(
    sum(carp.fit$carp.sol.path$v.zero.inds[, 1]),
    0
  )
})
test_that("CARP-returned sparsity pattern ends with full fusions (uniform weights)",{
  weight_mat <-matrix(1, nrow=nrow(presidential_speech), ncol=nrow(presidential_speech))
  carp.fit <- CARP(presidential_speech, weights = weight_mat, alg.type='carp')
  expect_true(
    all(carp.fit$carp.sol.path$v.zero.inds[, ncol(carp.fit$carp.sol.path$v.zero.inds)] == 1)
  )
})

test_that("CARPVIZ-returned sparsity pattern begins with no fusions (uniform weights)",{
  weight_mat <-matrix(1, nrow=nrow(presidential_speech), ncol=nrow(presidential_speech))
  carp.fit <- CARP(presidential_speech, alg.type = 'carpviz', weights = weight_mat)
  expect_equal(
    sum(carp.fit$carp.sol.path$v.zero.inds[, 1]),
    0
  )
})
test_that("CARPVIZ-returned sparsity pattern ends with full fusions (uniform weights)",{
  weight_mat <-matrix(1, nrow=nrow(presidential_speech), ncol=nrow(presidential_speech))
  carp.fit <- CARP(presidential_speech, alg.type='carpviz', weights = weight_mat)
  expect_true(
    all(carp.fit$carp.sol.path$v.zero.inds[, ncol(carp.fit$carp.sol.path$v.zero.inds)] == 1)
  )
})
