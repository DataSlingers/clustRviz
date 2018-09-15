context("Smoke tests for CARP")

test_that("CARP-VIZ [L2] works", {
  expect_no_error(carp_fit <- CARP(presidential_speech, alg.type="carpviz"))
  expect_no_error(print(carp_fit))
})

test_that("CARP-VIZ [L1] works", {
  expect_no_error(carp_fit <- CARP(presidential_speech, alg.type="carpvizl1"))
  expect_no_error(print(carp_fit))
})

test_that("CARP [L2] works", {
  expect_no_error(carp_fit <- CARP(presidential_speech, alg.type="carp", t=1.2))
  expect_no_error(CARP(presidential_speech, alg.type="carp", t=1.1))
  expect_no_error(CARP(presidential_speech, alg.type="carp", t=1.05))
  expect_no_error(print(carp_fit))
})

test_that("CARP [L1] works", {
  expect_no_error(carp_fit <- CARP(presidential_speech, alg.type="carpl1", t=1.2))
  expect_no_error(CARP(presidential_speech, alg.type="carpl1", t=1.1))
  expect_no_error(CARP(presidential_speech, alg.type="carpl1", t=1.05))
  expect_no_error(print(carp_fit))
})
