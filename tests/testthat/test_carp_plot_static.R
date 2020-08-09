context("Test Static CARP Plots")

## Add tests for static CARP plots
##
## Right now these are just "smoke" tests (i.e., runs in the default way without error)
## but we will add actual tests later (see GH #44)

test_that("CARP path plot works", {
  carp_fit <- CARP(presidential_speech)

  ## Main settings work
  expect_no_error(plot(carp_fit, type = "path"))
  expect_is(plot(carp_fit, type = "path"), "gg")
  expect_no_error(plot(carp_fit, type = "path", k = 3, axis = c("amount", "appropri")))
  expect_no_error(plot(carp_fit, type = "path", k = 3, axis = c("PC1", "PC3")))

  ## Must give at most one of `percent` or `k`
  expect_error(plot(carp_fit, type = "path", percent = 0.5, k = 3))

  ## Error checking on `percent` and `k`
  expect_error(plot(carp_fit, type = "path", percent = 1.5))
  expect_error(plot(carp_fit, type = "path", percent = -0.5))
  expect_error(plot(carp_fit, type = "path", percent = NA))
  expect_error(plot(carp_fit, type = "path", percent = c(0.25, 0.75)))

  expect_error(plot(carp_fit, type = "path", k = 3.5))
  expect_error(plot(carp_fit, type = "path", k = 0))
  expect_error(plot(carp_fit, type = "path", k = -1))
  expect_error(plot(carp_fit, type = "path", k = NROW(presidential_speech) + 1))

  ## Error on unknown arguments
  expect_error(plot(carp_fit, type = "path", 5))
  expect_error(plot(carp_fit, type = "path", a = 5))
})

test_that("CARP dendrogram plot works", {
  carp_fit <- CARP(presidential_speech)

  ## Main settings work
  expect_no_error(plot(carp_fit, type = "dendrogram"))
  expect_equal(plot(carp_fit, type = "dendrogram"), invisible(carp_fit))

  ## Must give at most one of `percent` or `k`
  expect_error(plot(carp_fit, type = "dendrogram", percent = 0.5, k = 3))

  ## Error checking on `percent` and `k`
  expect_error(plot(carp_fit, type = "dendrogram", percent = 1.5))
  expect_error(plot(carp_fit, type = "dendrogram", percent = -0.5))
  expect_error(plot(carp_fit, type = "dendrogram", percent = NA))
  expect_error(plot(carp_fit, type = "dendrogram", percent = c(0.25, 0.75)))

  expect_error(plot(carp_fit, type = "dendrogram", k = 3.5))
  expect_error(plot(carp_fit, type = "dendrogram", k = 0))
  expect_error(plot(carp_fit, type = "dendrogram", k = -1))
  expect_error(plot(carp_fit, type = "dendrogram", k = NROW(presidential_speech) + 1))
})

test_that("CARP heatmaply (javascript) plot works", {
  carp_fit <- CARP(presidential_speech)

  ## Main settings work
  expect_no_error(plot(carp_fit, type = "heatmap", interactive = TRUE))
  expect_is(plot(carp_fit, type = "heatmap", interactive = TRUE), "plotly")
  expect_no_error(plot(carp_fit, type = "heatmap", interactive = TRUE, k = 3))
  expect_no_error(plot(carp_fit, type = "heatmap", interactive = TRUE, percent = 0.5))

  ## Must give at most one of `percent` or `k`
  expect_error(plot(carp_fit, type = "heatmap", interactive = TRUE, percent = 0.5, k = 3))

  ## Error checking on `percent` and `k`
  expect_error(plot(carp_fit, type = "heatmap", interactive = TRUE, percent = 1.5))
  expect_error(plot(carp_fit, type = "heatmap", interactive = TRUE, percent = -0.5))
  expect_error(plot(carp_fit, type = "heatmap", interactive = TRUE, percent = NA))
  expect_error(plot(carp_fit, type = "heatmap", interactive = TRUE, percent = c(0.25, 0.75)))

  expect_error(plot(carp_fit, type = "heatmap", interactive = TRUE, k = 3.5))
  expect_error(plot(carp_fit, type = "heatmap", interactive = TRUE, k = 0))
  expect_error(plot(carp_fit, type = "heatmap", interactive = TRUE, k = -1))
  expect_error(plot(carp_fit, type = "heatmap", interactive = TRUE, k = NROW(presidential_speech) + 1))
})
