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

test_that("CBASS interactive dendrogram plot works", {
  cbass_fit <- CBASS(presidential_speech)

  ## static row
  expect_no_error(plot(cbass_fit, type = "row.dendrogram", interactive = T, dynamic = F))
  expect_no_error(plot(cbass_fit, type = "row.dendrogram", interactive = T, dynamic = F, k.row = 3))
  expect_no_error(plot(cbass_fit, type = "row.dendrogram", interactive = T, dynamic = F, percent = 0.25))

  ## static col
  expect_no_error(plot(cbass_fit, type = "col.dendrogram", interactive = T, dynamic = F))
  expect_no_error(plot(cbass_fit, type = "col.dendrogram", interactive = T, dynamic = F, k.col = 3))
  expect_no_error(plot(cbass_fit, type = "col.dendrogram", interactive = T, dynamic = F, percent = 0.25))

  ## dynamic row
  expect_no_error(plot(cbass_fit, type = "row.dendrogram", interactive = T, dynamic = T))

  ## dynamic col
  expect_no_error(plot(cbass_fit, type = "col.dendrogram", interactive = T, dynamic = T))
})

test_that("CBASS interactive path plot works", {
  cbass_fit <- CBASS(presidential_speech)

  ## static row
  expect_no_error(plot(cbass_fit, type = "row.path", interactive = T, dynamic = F))
  expect_no_error(plot(cbass_fit, type = "row.path", interactive = T, dynamic = F, k.row = 3))
  expect_no_error(plot(cbass_fit, type = "row.path", interactive = T, dynamic = F, percent = 0.25))

  ## static col
  expect_no_error(plot(cbass_fit, type = "col.path", interactive = T, dynamic = F))
  expect_no_error(plot(cbass_fit, type = "col.path", interactive = T, dynamic = F, k.col = 3))
  expect_no_error(plot(cbass_fit, type = "col.path", interactive = T, dynamic = F, percent = 0.25))

  ## dynamic row
  expect_no_error(plot(cbass_fit, type = "row.path", interactive = T, dynamic = T))

  ## dynamic col
  expect_no_error(plot(cbass_fit, type = "col.path", interactive = T, dynamic = T))
})

test_that("CBASS interactive heatmap works", {
  cbass_fit <- CBASS(presidential_speech[1:5,1:5])

  ## dynamic
  expect_no_error(plot(cbass_fit, type = "heatmap", interactive = T, dynamic = T))
})

