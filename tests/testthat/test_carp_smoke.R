context("Smoke tests for CARP")

test_that("CARP-VIZ [L2] works", {
  ## Also smoke test status printing here
  expect_no_error(carp_fit <- CARP(presidential_speech, back_track = TRUE, status = TRUE))
  expect_no_error(print(carp_fit))
})

test_that("CARP-VIZ [L1] works", {
  expect_no_error(carp_fit <- CARP(presidential_speech, back_track = TRUE, norm = 1))
  expect_no_error(print(carp_fit))
})

test_that("CARP [L2] works", {
  expect_no_error(carp_fit <- CARP(presidential_speech, back_track = FALSE, t = 1.2))
  expect_no_error(CARP(presidential_speech, back_track = FALSE, t = 1.1))
  expect_no_error(CARP(presidential_speech, back_track = FALSE, t = 1.05))
  expect_no_error(print(carp_fit))
})

test_that("CARP [L1] works", {
  expect_no_error(carp_fit <- CARP(presidential_speech, back_track = FALSE, t = 1.2, norm = 1))
  expect_no_error(CARP(presidential_speech, back_track = FALSE, t =1.1, norm = 1))
  expect_no_error(CARP(presidential_speech, back_track = FALSE, t = 1.05, norm = 1))
  expect_no_error(print(carp_fit))
})
