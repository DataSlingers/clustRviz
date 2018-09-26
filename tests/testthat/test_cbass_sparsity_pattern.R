context("test CBASS sparsity pattern")



test_that("CBASS-returned observation sparsity pattern begins with no fusions",{
  cbass.fit <- CBASS(presidential_speech, alg.type = 'cbass')
  expect_equal(
    sum(cbass.fit$cbass.sol.path$v.col.zero.inds[, 1]),
    0
  )
})
test_that("CBASS-returned observation sparsity pattern ends with full fusions",{
  cbass.fit <- CBASS(presidential_speech, alg.type = 'cbass')
  expect_true(
    all(cbass.fit$cbass.sol.path$v.col.zero.inds[, ncol(cbass.fit$cbass.sol.path$v.col.zero.inds)] == 1)
  )
})
test_that("CBASS-returned variable sparsity pattern begins with no fusions",{
  cbass.fit <- CBASS(presidential_speech, alg.type = 'cbass')
  expect_equal(
    sum(cbass.fit$cbass.sol.path$v.row.zero.inds[, 1]),
    0
  )
})
test_that("CBASS-returned variable sparsity pattern ends with full fusions",{
  cbass.fit <- CBASS(presidential_speech, alg.type='cbass')
  expect_true(
    all(cbass.fit$cbass.sol.path$v.row.zero.inds[, ncol(cbass.fit$cbass.sol.path$v.row.zero.inds)] == 1)
  )
})






test_that("CBASSVIZ-returned observation sparsity pattern begins with no fusions",{
  cbass.fit <- CBASS(presidential_speech, alg.type = 'cbassviz')
  expect_equal(
    sum(cbass.fit$cbass.sol.path$v.col.zero.inds[, 1]),
    0
  )
})
test_that("CBASSVIZ-returned observation sparsity pattern ends with full fusions",{
  cbass.fit <- CBASS(presidential_speech, alg.type = 'cbassviz')
  expect_true(
    all(cbass.fit$cbass.sol.path$v.col.zero.inds[, ncol(cbass.fit$cbass.sol.path$v.col.zero.inds)]==1)
  )
})
test_that("CBASSVIZ-returned variable sparsity pattern begins with no fusions",{
  cbass.fit <- CBASS(presidential_speech, alg.type = 'cbassviz')
  expect_equal(
    sum(cbass.fit$cbass.sol.path$v.row.zero.inds[, 1]),
    0
  )
})
test_that("CBASSVIZ-returned variable sparsity pattern ends with full fusions",{
  cbass.fit <- CBASS(presidential_speech, alg.type='cbassviz')
  expect_true(
    all(cbass.fit$cbass.sol.path$v.row.zero.inds[, ncol(cbass.fit$cbass.sol.path$v.row.zero.inds)] == 1)
  )
})





test_that("CBASS-returned observation sparsity pattern begins with no fusions (uniform weights)",{
  weight_mat_obs <-matrix(1, nrow=nrow(presidential_speech), ncol=nrow(presidential_speech))
  weight_mat_var <-matrix(1, nrow=ncol(presidential_speech), ncol=ncol(presidential_speech))
  cbass.fit <- CBASS(presidential_speech, alg.type = 'cbass', obs_weights = weight_mat_obs, var_weights = weight_mat_var)
  expect_equal(
    sum(cbass.fit$cbass.sol.path$v.col.zero.inds[, 1]),
    0
  )
})
test_that("CBASS-returned observation sparsity pattern ends with full fusions (uniform weights)",{
  weight_mat_obs <-matrix(1, nrow=nrow(presidential_speech), ncol=nrow(presidential_speech))
  weight_mat_var <-matrix(1, nrow=ncol(presidential_speech), ncol=ncol(presidential_speech))
  cbass.fit <- CBASS(presidential_speech, alg.type = 'cbass', obs_weights = weight_mat_obs, var_weights = weight_mat_var)
  expect_true(
    all(cbass.fit$cbass.sol.path$v.col.zero.inds[, ncol(cbass.fit$cbass.sol.path$v.col.zero.inds)] == 1)
  )
})
test_that("CBASS-returned variable sparsity pattern begins with no fusions (uniform weights)",{
  weight_mat_obs <-matrix(1, nrow=nrow(presidential_speech), ncol=nrow(presidential_speech))
  weight_mat_var <-matrix(1, nrow=ncol(presidential_speech), ncol=ncol(presidential_speech))
  cbass.fit <- CBASS(presidential_speech, alg.type = 'cbass', obs_weights = weight_mat_obs, var_weights = weight_mat_var)
  expect_equal(
    sum(cbass.fit$cbass.sol.path$v.row.zero.inds[, 1]),
    0
  )
})
test_that("CBASS-returned variable sparsity pattern ends with full fusions (uniform weights)",{
  weight_mat_obs <-matrix(1, nrow=nrow(presidential_speech), ncol=nrow(presidential_speech))
  weight_mat_var <-matrix(1, nrow=ncol(presidential_speech), ncol=ncol(presidential_speech))
  cbass.fit <- CBASS(presidential_speech, alg.type = 'cbass', obs_weights = weight_mat_obs, var_weights = weight_mat_var)
  expect_true(
    all(cbass.fit$cbass.sol.path$v.row.zero.inds[, ncol(cbass.fit$cbass.sol.path$v.row.zero.inds)]==1)
  )
})




test_that("CBASSVIZ-returned observation sparsity pattern begins with no fusions (uniform weights)",{
  weight_mat_obs <-matrix(1, nrow=nrow(presidential_speech), ncol=nrow(presidential_speech))
  weight_mat_var <-matrix(1, nrow=ncol(presidential_speech), ncol=ncol(presidential_speech))
  cbass.fit <- CBASS(presidential_speech, alg.type = 'cbassviz', obs_weights = weight_mat_obs, var_weights = weight_mat_var)
  expect_equal(
    sum(cbass.fit$cbass.sol.path$v.col.zero.inds[, 1]),
    0
  )
})
test_that("CBASSVIZ-returned observation sparsity pattern ends with full fusions (uniform weights)",{
  weight_mat_obs <-matrix(1, nrow=nrow(presidential_speech), ncol=nrow(presidential_speech))
  weight_mat_var <-matrix(1, nrow=ncol(presidential_speech), ncol=ncol(presidential_speech))
  cbass.fit <- CBASS(presidential_speech, alg.type = 'cbassviz', obs_weights = weight_mat_obs, var_weights = weight_mat_var)
  expect_true(
    all(cbass.fit$cbass.sol.path$v.col.zero.inds[, ncol(cbass.fit$cbass.sol.path$v.col.zero.inds)] == 1)
  )
})
test_that("CBASSVIZ-returned variable sparsity pattern begins with no fusions (uniform weights)",{
  weight_mat_obs <-matrix(1, nrow=nrow(presidential_speech), ncol=nrow(presidential_speech))
  weight_mat_var <-matrix(1, nrow=ncol(presidential_speech), ncol=ncol(presidential_speech))
  cbass.fit <- CBASS(presidential_speech, alg.type = 'cbassviz', obs_weights = weight_mat_obs, var_weights = weight_mat_var)
  expect_equal(
    sum(cbass.fit$cbass.sol.path$v.row.zero.inds[, 1]),
    0
  )
})
test_that("CBASSVIZ-returned variable sparsity pattern ends with full fusions (uniform weights)",{
  weight_mat_obs <-matrix(1, nrow=nrow(presidential_speech), ncol=nrow(presidential_speech))
  weight_mat_var <-matrix(1, nrow=ncol(presidential_speech), ncol=ncol(presidential_speech))
  cbass.fit <- CBASS(presidential_speech, alg.type = 'cbassviz', obs_weights = weight_mat_obs, var_weights = weight_mat_var)
  expect_true(
    all(cbass.fit$cbass.sol.path$v.row.zero.inds[, ncol(cbass.fit$cbass.sol.path$v.row.zero.inds)]==1)
  )
})

