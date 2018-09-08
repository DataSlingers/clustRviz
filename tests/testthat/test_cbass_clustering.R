context("test clustering.CBASS")

test_that("Works with user `k.obs`", {
  cbass_fit <- CBASS(presidential_speech)

  clustering_result <- clustering(cbass_fit, k.obs = 5)

  expect_equal(length(table(clustering_result$clustering.assignment.obs)), 5)
  expect_equal(num_unique_rows(clustering_result$cluster.mean.matrix), 5)

  expect_error(clustering(cbass_fit, k.obs = 0))
  expect_error(clustering(cbass_fit, k.obs = NROW(presidential_speech) + 1))

  expect_equal(length(clustering_result$cluster.assignment.obs, NROW(presidential_speech)))
  expect_equal(length(clustering_result$cluster.assignment.var, NCOL(presidential_speech)))
})

test_that("Works with user `k.var`", {
  cbass_fit <- CBASS(presidential_speech)

  clustering_result <- clustering(cbass_fit, k.var = 5)

  expect_equal(length(table(clustering_result$clustering.assignment.var)), 5)
  expect_equal(num_unique_cols(clustering_result$cluster.mean.matrix), 5)

  expect_error(clustering(cbass_fit, k.var = 0))
  expect_error(clustering(cbass_fit, k.var = NCOL(presidential_speech) + 1))

  expect_equal(length(clustering_result$cluster.assignment.obs, NROW(presidential_speech)))
  expect_equal(length(clustering_result$cluster.assignment.var, NCOL(presidential_speech)))
})

test_that("Works with user `percent`", {
  cbass_fit <- CBASS(presidential_speech)

  clustering_result <- clustering(cbass_fit, percent = 0)
  # FIXME: This might be a bug - see GitHub #10
  expect_equal(clustering_result$cluster.mean.matrix,
               presidential_speech - mean(presidential_speech))

  clustering_result <- clustering(cbass_fit, percent = 1)
  # FIXME: This might be a bug - see GitHub #10
  expect_equal(clustering_result$cluster.mean.matrix,
               0 * presidential_speech)

  expect_error(clustering(cbass_fit, percent = +1.2))
  expect_error(clustering(cbass_fit, percent = -0.2))
})

test_that("Exactly one argument required", {
  cbass_fit <- CBASS(presidential_speech)

  expect_error(clustering(cbass_fit))
  expect_error(clustering(cbass_fit, k.obs = 5, k.var = 5))
  expect_error(clustering(cbass_fit, k.obs = 5, percent = 0.5))
  expect_error(clustering(cbass_fit, k.var = 5, percent = 0.5))
  expect_error(clustering(cbass_fit, k.obs = 5, k.var = 5, percent = 0.5))
})

test_that("Observation labels returned by clustering.CBASS have correct length (=NROW(X))", {
  cbass_fit <- CBASS(presidential_speech)

  # Percent
  expect_equal(NROW(presidential_speech),
               length(clustering(cbass_fit, percent = 0.2)$clustering.assignment.obs))

  # k.obs
  expect_equal(NROW(presidential_speech),
               length(clustering(cbass_fit, k.obs = 10)$clustering.assignment.obs))

  # k.var
  expect_equal(NROW(presidential_speech),
               length(clustering(cbass_fit, k.var = 10)$clustering.assignment.obs))
})

test_that("Variable labels returned by clustering.CBASS have correct length (=NCOL(X))", {
  cbass_fit <- CBASS(presidential_speech)

  # Percent
  expect_equal(NCOL(presidential_speech),
               length(clustering(cbass_fit, percent = 0.2)$clustering.assignment.var))

  # k.obs
  expect_equal(NCOL(presidential_speech),
               length(clustering(cbass_fit, k.obs = 10)$clustering.assignment.var))

  # k.var
  expect_equal(NCOL(presidential_speech),
               length(clustering(cbass_fit, k.var = 10)$clustering.assignment.var))

})

test_that("Clustered data matrix returned by clustering.CBASS has correct dimensions (=dim(X))",{
  cbass_fit <- CBASS(presidential_speech)

  # Percent
  expect_equal(dim(presidential_speech),
               dim(clustering(cbass_fit, percent = .2)$cluster.mean.matrix))

  # k.obs
  expect_equal(dim(presidential_speech),
               dim(clustering(cbass_fit, k.obs = 10)$cluster.mean.matrix))

  # k.var
  expect_equal(dim(presidential_speech),
               dim(clustering(cbass_fit, k.var = 10)$cluster.mean.matrix))
})

test_that("CBASS clustering (w/ X.center.global == TRUE) returns zero matrix at full regularization", {
  cbass_fit <- CBASS(presidential_speech, X.center.global = TRUE)
  expect_equal(clustering(cbass_fit, percent = 1)$cluster.mean.matrix,
               0 * presidential_speech)
})

test_that("CBASS clustering (w/ X.center.global == FALSE) returns global mean at full regularization", {
  cbass_fit <- CBASS(presidential_speech, X.center.global = FALSE)
  expect_equal(clustering(cbass_fit, percent = 1)$cluster.mean.matrix,
               0 * presidential_speech + mean(presidential_speech))
})

test_that("CBASS clustering (w/ X.center.global == TRUE) returns (centered) original data at 0 regularization", {
  cbass_fit <- CBASS(presidential_speech, X.center.global = TRUE)
  expect_equal(clustering(cbass_fit, percent = 0)$cluster.mean.matrix,
               presidential_speech - mean(presidential_speech))
})

test_that("CBASS clustering (w/ X.center.global == FALSE) returns original data at 0 regularization", {
  cbass_fit <- CBASS(presidential_speech, X.center.global = FALSE)
  expect_equal(clustering(cbass_fit, percent = 0)$cluster.mean.matrix,
               presidential_speech)
})

test_that("CBASS clustering cluster.mean.matrix respects obsveration/row labels", {
  cbass_fit <- CBASS(presidential_speech)
  expect_equal(rownames(presidential_speech),
               rownames(clustering(cbass_fit, percent = .3)$cluster.mean.matrix))
})

test_that("CBASS clustering cluster.mean.matrix respects variable/column labels",{
  cbass_fit <- CBASS(presidential_speech)
  expect_equal(colnames(presidential_speech),
               colnames(clustering(cbass_fit, percent = .3)$cluster.mean.matrix))
})
