context("CBASS() Error Handling")

test_that("CBASS() fails with non-finite numerical input", {
  ps <- presidential_speech

  ps[1,1] <- NA; expect_error(CBASS(ps))
  ps[1,1] <- NaN; expect_error(CBASS(ps))
  ps[1,1] <- Inf; expect_error(CBASS(ps))
})

test_that("CBASS() errors early with incorrect input", {
  # ADMM Relaxation parameter must be positive
  expect_error(CBASS(presidential_speech, rho = 0))
  expect_error(CBASS(presidential_speech, rho = -1))

  # Pre-processing parameters must be boolean flags
  expect_error(CBASS(presidential_speech, X.center.global = NA))
  expect_error(CBASS(presidential_speech, X.center.global = c(TRUE, FALSE)))

  # Must use at least one core
  expect_error(CBASS(presidential_speech, ncores = NA))
  expect_error(CBASS(presidential_speech, ncores = c(1, 5)))
  expect_error(CBASS(presidential_speech, ncores = 0L))
  expect_error(CBASS(presidential_speech, ncores = -1L))

  # Must use a known algorithm
  expect_error(CBASS(presidential_speech, alg.type = "unknown"))
  expect_error(CBASS(presidential_speech, alg.type = NA))
  expect_error(CBASS(presidential_speech, alg.type = "CBASS"))

  # Must use a t > 1
  expect_error(CBASS(presidential_speech, t = 1))
  expect_error(CBASS(presidential_speech, t = 0))
  expect_error(CBASS(presidential_speech, t = -3))
  expect_error(CBASS(presidential_speech, t = NA))
  expect_error(CBASS(presidential_speech, t = c(1.3, 1.2)))

  # Check max.iter
  expect_error(CBASS(presidential_speech, max.iter = 1))
  expect_error(CBASS(presidential_speech, max.iter = 75.5))
  expect_error(CBASS(presidential_speech, max.iter = 0))
  expect_error(CBASS(presidential_speech, max.iter = -3))
  expect_error(CBASS(presidential_speech, max.iter = NA))
  expect_error(CBASS(presidential_speech, max.iter = c(5, 10)))

  # Check burn.in
  expect_error(CBASS(presidential_speech, burn.in = 75.5))
  expect_error(CBASS(presidential_speech, burn.in = 0))
  expect_error(CBASS(presidential_speech, burn.in = -3))
  expect_error(CBASS(presidential_speech, burn.in = NA))
  expect_error(CBASS(presidential_speech, burn.in = c(5, 10)))
  expect_error(CBASS(presidential_speech, burn.in = 50, max.iter = 25))

  # Fail on unknown flags
  expect_error(CBASS(presidential_speech, flag="unknown"), regexp = "flag")
  expect_error(CBASS(presidential_speech, "value"), regexp = "Unknown")
})
