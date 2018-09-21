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

  ## Can only plot pre-calculated PCs right now
  ## FIXME -- See GH #24
  expect_error(plot(carp_fit, type = "path", k = 3, axis = c("amount", "appropri")))
  expect_error(plot(carp_fit, type = "path", k = 3, axis = c("PC1", "PC5")))

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

  ## Error checking on `show_clusters`
  expect_error(plot(carp_fit, type = "path", show_clusters = NA))
  expect_error(plot(carp_fit, type = "path", show_clusters = c(TRUE, FALSE)))

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

  ## Error checking on `show_clusters`
  expect_error(plot(carp_fit, type = "dendrogram", show_clusters = NA))
  expect_error(plot(carp_fit, type = "dendrogram", show_clusters = c(TRUE, FALSE)))

  ## Error checking on two specially handled arguments
  expect_error(plot(carp_fit, type = "dendrogram", k = 3, dend.branch.width = 0))
  expect_error(plot(carp_fit, type = "dendrogram", k = 3, dend.branch.width = -2))
  expect_error(plot(carp_fit, type = "dendrogram", k = 3, dend.labels.cex   = 0))
  expect_error(plot(carp_fit, type = "dendrogram", k = 3, dend.labels.cex   = -2))

  ## Error if no regularization given, but show_clusters = TRUE
  expect_error(plot(carp_fit, type = "dendrogram", show_clusters = TRUE))
})
