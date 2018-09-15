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
