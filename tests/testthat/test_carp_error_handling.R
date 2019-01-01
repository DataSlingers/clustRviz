context("CARP() Error Handling")

test_that("CARP() fails with non-finite numerical input", {
  ps <- presidential_speech

  ps[1,1] <- NA; expect_error(CARP(ps))
  ps[1,1] <- NaN; expect_error(CARP(ps))
  ps[1,1] <- Inf; expect_error(CARP(ps))
})

test_that("CARP() errors early with incorrect input", {
  # Pre-processing parameters must be boolean flags
  expect_error(CARP(presidential_speech, X.center = NA))
  expect_error(CARP(presidential_speech, X.center = c(TRUE, FALSE)))

  expect_error(CARP(presidential_speech, X.scale = NA))
  expect_error(CARP(presidential_speech, X.scale = c(TRUE, FALSE)))

  # Must use at least one core
  expect_error(CARP(presidential_speech, ncores = NA))
  expect_error(CARP(presidential_speech, ncores = c(1, 5)))
  expect_error(CARP(presidential_speech, ncores = 0L))
  expect_error(CARP(presidential_speech, ncores = -1L))

  # Check `exact` argument
  expect_error(CARP(presidential_speech, exact = "unknown"))
  expect_error(CARP(presidential_speech, exact = NA))
  expect_error(CARP(presidential_speech, exact = c(TRUE, FALSE)))
  expect_error(CARP(presidential_speech, exact = 1L))

  # Check `back_track` argument
  expect_error(CARP(presidential_speech, back_track = "unknown"))
  expect_error(CARP(presidential_speech, back_track = NA))
  expect_error(CARP(presidential_speech, back_track = c(TRUE, FALSE)))
  expect_error(CARP(presidential_speech, back_track = 1L))

  # Must use a t > 1
  expect_error(CARP(presidential_speech, t = 1))
  expect_error(CARP(presidential_speech, t = 0))
  expect_error(CARP(presidential_speech, t = -3))
  expect_error(CARP(presidential_speech, t = NA))
  expect_error(CARP(presidential_speech, t = c(1.3, 1.2)))

  # Fail on unknown flags
  expect_error(CARP(presidential_speech, flag="unknown"), regexp = "flag")
  expect_error(CARP(presidential_speech, "value"), regexp = "Unknown")
})
