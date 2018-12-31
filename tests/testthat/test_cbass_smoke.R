context("Smoke tests for CBASS")

test_that("CBASS-VIZ [L2] works", {
  expect_no_error(cbass_fit <- CBASS(presidential_speech, back_track = TRUE))
  expect_no_error(print(cbass_fit))
})

test_that("CBASS-VIZ [L1] works", {
  expect_no_error(cbass_fit <- CBASS(presidential_speech, back_track = TRUE, norm = 1))
  expect_no_error(print(cbass_fit))
})

test_that("CBASS [L2] works", {
  expect_no_error(cbass_fit <- CBASS(presidential_speech, back_track = FALSE, t = 1.2))
  expect_no_error(CBASS(presidential_speech, back_track = FALSE, t = 1.1))
  expect_no_error(CBASS(presidential_speech, back_track = FALSE, t = 1.05))
  expect_no_error(print(cbass_fit))
})

test_that("CBASS [L1] works", {
  expect_no_error(cbass_fit <- CBASS(presidential_speech, back_track = FALSE, t = 1.2, norm = 1))
  expect_no_error(CBASS(presidential_speech, back_track = FALSE, t = 1.1, norm = 1))
  expect_no_error(CBASS(presidential_speech, back_track = FALSE, t = 1.05, norm = 1))
  expect_no_error(print(cbass_fit))
})

test_that("CBASS runs with X.center.global = FALSE & progress printing", {
  expect_no_error(CBASS(presidential_speech, X.center.global = FALSE, status = TRUE))
})
