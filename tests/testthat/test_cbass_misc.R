context("Test Miscellaneous Features of CBASS")

test_that("CBASS creates unique names if needed", {
  ## Give "Grover the Good" due credit
  pres2 <- rbind(presidential_speech,
                 jitter(presidential_speech["Grover Cleveland", , drop=FALSE]))

  cbass_fit <- CBASS(pres2, X.center.global = FALSE)

  expect_equal(rownames(cbass_fit$X)[45], "Grover Cleveland_1")

  expect_equal(head(get_clustered_data(cbass_fit, percent = 0), -1),
               presidential_speech)
})

test_that("CBASS supports factor labels", {
  X <- matrix(rnorm(9), 3, 3)

  cbass_fit <- CBASS(X, row_labels = factor(c("a", "b", "c")),
                        col_labels = factor(c("d", "e", "f")))

  expect_equal(rownames(cbass_fit$X), c("a", "b", "c"))
  expect_equal(colnames(cbass_fit$X), c("d", "e", "f"))
})

test_that("CBASS stores mean of original data", {
  cbass_fit <- CBASS(presidential_speech, X.center.global = TRUE)
  expect_equal(mean(presidential_speech), cbass_fit$mean_adjust)

  cbass_fit <- CBASS(presidential_speech, X.center.global = FALSE)
  expect_equal(0, cbass_fit$mean_adjust)
})
