context("Test Miscellaneous Features of CARP")

test_that("CARP creates unique names if needed", {
  ## Give "Grover the Good" due credit
  pres2 <- rbind(presidential_speech,
                 jitter(presidential_speech["Grover Cleveland", , drop=FALSE]))

  carp_fit <- CARP(pres2)

  expect_equal(rownames(carp_fit$X)[45], "Grover Cleveland_1")
})

test_that("CARP supports factor labels", {
  X <- matrix(rnorm(9), 3, 3)
  colnames(X) <- c("a", "b", "c")

  carp_fit <- CARP(X, labels = factor(c("a", "b", "c")))

  expect_equal(rownames(carp_fit$X), c("a", "b", "c"))
})

test_that("CARP stores scale factors", {
  carp_fit_std    <- CARP(presidential_speech, X.center = TRUE,  X.scale = TRUE)

  expect_equal(carp_fit_std$center_vector, colMeans(presidential_speech))
  expect_equal(carp_fit_std$scale_vector, apply(presidential_speech, 2, sd))

  carp_fit_no_std <- CARP(presidential_speech, X.center = FALSE, X.scale = FALSE)

  expect_equal(carp_fit_no_std$scale_vector, rep(1, NCOL(presidential_speech)))
  expect_equal(carp_fit_no_std$center_vector, rep(0, NCOL(presidential_speech)))
})
