context("Test Static CBASS Plots")

## Add tests for static CBASS plots
##
## Right now these are just "smoke" tests (i.e., runs in the default way without error)
## but we will add actual tests later (see GH #44)

test_that("CBASS path plot works for observations", {
  cbass_fit <- CBASS(presidential_speech)

  ## Main settings work
  expect_no_error(plot(cbass_fit, type = "obs.path"))
  expect_is(plot(cbass_fit, type = "obs.path"), "gg")

  ## Can only plot pre-calculated PCs right now
  ## FIXME -- See GH #24
  expect_error(plot(cbass_fit, type = "obs.path", k.obs  = 3, axis = c("amount", "appropri")))
  expect_error(plot(cbass_fit, type = "obs.path", k.obs  = 3, axis = c("PC1", "PC5")))

  ## Must give at most one of `percent`, `k.obs`, `k.var`
  expect_error(plot(cbass_fit, type = "obs.path", percent = 0.5, k.obs = 3))
  expect_error(plot(cbass_fit, type = "obs.path", percent = 0.5, k.var = 3))
  expect_error(plot(cbass_fit, type = "obs.path", k.obs = 3, k.var = 3))

  ## Error checking on `percent` and `k`
  expect_error(plot(cbass_fit, type = "obs.path", percent = 1.5))
  expect_error(plot(cbass_fit, type = "obs.path", percent = -0.5))
  expect_error(plot(cbass_fit, type = "obs.path", percent = NA))
  expect_error(plot(cbass_fit, type = "obs.path", percent = c(0.25, 0.75)))

  expect_error(plot(cbass_fit, type = "obs.path", k.obs = 3.5))
  expect_error(plot(cbass_fit, type = "obs.path", k.obs = 0))
  expect_error(plot(cbass_fit, type = "obs.path", k.obs = -1))
  expect_error(plot(cbass_fit, type = "obs.path", k.obs = NROW(presidential_speech) + 1))

  expect_error(plot(cbass_fit, type = "obs.path", k.var = 3.5))
  expect_error(plot(cbass_fit, type = "obs.path", k.var = 0))
  expect_error(plot(cbass_fit, type = "obs.path", k.var = -1))
  expect_error(plot(cbass_fit, type = "obs.path", k.var = NCOL(presidential_speech) + 1))

  ## Error checking on `show_clusters`
  expect_error(plot(cbass_fit, type = "obs.path", show_clusters = NA))
  expect_error(plot(cbass_fit, type = "obs.path", show_clusters = c(TRUE, FALSE)))

  ## Error on unknown arguments
  expect_error(plot(cbass_fit, type = "obs.path", 5))
  expect_error(plot(cbass_fit, type = "obs.path", a = 5))

  ## Error if `show_clusters = TRUE` but no regularization given
  expect_error(plot(cbass_fit, type = "obs.path", show_clusters = TRUE))
})

test_that("CBASS path plot works for variables", {
  cbass_fit <- CBASS(presidential_speech)

  ## Main settings work
  expect_no_error(plot(cbass_fit, type = "var.path", k.obs = 3))
  expect_is(plot(cbass_fit, type = "var.path", k.obs  = 3), "gg")

  ## Can only plot pre-calculated PCs right now
  ## FIXME -- See GH #24
  expect_error(plot(cbass_fit, type = "var.path", k.obs  = 3, axis = c("amount", "appropri")))
  expect_error(plot(cbass_fit, type = "var.path", k.obs  = 3, axis = c("PC1", "PC5")))

  ## Must give at most one of `percent`, `k.obs`, `k.var`
  expect_error(plot(cbass_fit, type = "var.path", percent = 0.5, k.obs = 3))
  expect_error(plot(cbass_fit, type = "var.path", percent = 0.5, k.var = 3))
  expect_error(plot(cbass_fit, type = "var.path", k.obs = 3, k.var = 3))

  ## Error checking on `percent` and `k`
  expect_error(plot(cbass_fit, type = "var.path", percent = 1.5))
  expect_error(plot(cbass_fit, type = "var.path", percent = -0.5))
  expect_error(plot(cbass_fit, type = "var.path", percent = NA))
  expect_error(plot(cbass_fit, type = "var.path", percent = c(0.25, 0.75)))

  expect_error(plot(cbass_fit, type = "var.path", k.obs = 3.5))
  expect_error(plot(cbass_fit, type = "var.path", k.obs = 0))
  expect_error(plot(cbass_fit, type = "var.path", k.obs = -1))
  expect_error(plot(cbass_fit, type = "var.path", k.obs = NROW(presidential_speech) + 1))

  expect_error(plot(cbass_fit, type = "var.path", k.var = 3.5))
  expect_error(plot(cbass_fit, type = "var.path", k.var = 0))
  expect_error(plot(cbass_fit, type = "var.path", k.var = -1))
  expect_error(plot(cbass_fit, type = "var.path", k.var = NCOL(presidential_speech) + 1))

  ## Error checking on `show_clusters`
  expect_error(plot(cbass_fit, type = "var.path", show_clusters = NA))
  expect_error(plot(cbass_fit, type = "var.path", show_clusters = c(TRUE, FALSE)))

  ## Error on unknown arguments
  expect_error(plot(cbass_fit, type = "var.path", 5))
  expect_error(plot(cbass_fit, type = "var.path", a = 5))

  ## Error if `show_clusters = TRUE` but no regularization given
  expect_error(plot(cbass_fit, type = "var.path", show_clusters = TRUE))
})

test_that("CBASS dendrogram plot works for observations", {
  cbass_fit <- CBASS(presidential_speech)

  ## Main settings work
  expect_no_error(plot(cbass_fit, type = "obs.dendrogram", k.obs = 3))
  expect_equal(plot(cbass_fit, type = "obs.dendrogram", k.obs = 3), invisible(cbass_fit))

  ## Must give at most one of `percent`, `k.obs`, `k.var`
  expect_error(plot(cbass_fit, type = "obs.dendrogram", percent = 0.5, k.obs = 3))
  expect_error(plot(cbass_fit, type = "obs.dendrogram", percent = 0.5, k.var = 3))
  expect_error(plot(cbass_fit, type = "obs.dendrogram", k.obs = 3, k.var = 3))

  ## Error checking on `percent`, `k.obs`, `k.var`
  expect_error(plot(cbass_fit, type = "obs.dendrogram", percent = 1.5))
  expect_error(plot(cbass_fit, type = "obs.dendrogram", percent = -0.5))
  expect_error(plot(cbass_fit, type = "obs.dendrogram", percent = NA))
  expect_error(plot(cbass_fit, type = "obs.dendrogram", percent = c(0.25, 0.75)))

  expect_error(plot(cbass_fit, type = "obs.dendrogram", k.obs = 3.5))
  expect_error(plot(cbass_fit, type = "obs.dendrogram", k.obs = 0))
  expect_error(plot(cbass_fit, type = "obs.dendrogram", k.obs = -1))
  expect_error(plot(cbass_fit, type = "obs.dendrogram", k.obs = NROW(presidential_speech) + 1))

  expect_error(plot(cbass_fit, type = "obs.dendrogram", k.var = 3.5))
  expect_error(plot(cbass_fit, type = "obs.dendrogram", k.var = 0))
  expect_error(plot(cbass_fit, type = "obs.dendrogram", k.var = -1))
  expect_error(plot(cbass_fit, type = "obs.dendrogram", k.var = NCOL(presidential_speech) + 1))

  ## Error checking on `show_clusters`
  expect_error(plot(cbass_fit, type = "obs.dendrogram", show_clusters = NA))
  expect_error(plot(cbass_fit, type = "obs.dendrogram", show_clusters = c(TRUE, FALSE)))

  ## Error checking on two specially handled arguments
  expect_error(plot(cbass_fit, type = "obs.dendrogram", k.obs = 3, dend.branch.width = 0))
  expect_error(plot(cbass_fit, type = "obs.dendrogram", k.obs = 3, dend.branch.width = -2))
  expect_error(plot(cbass_fit, type = "obs.dendrogram", k.obs = 3, dend.labels.cex   = 0))
  expect_error(plot(cbass_fit, type = "obs.dendrogram", k.obs = 3, dend.labels.cex   = -2))

  ## No cluster specification required if `show_clusters` is FALSE
  expect_no_error(plot(cbass_fit, type = "obs.dendrogram", show_clusters = FALSE))

  ## Error if `show_clusters = TRUE` but no regularization given
  expect_error(plot(cbass_fit, type = "obs.dendrogram", show_clusters = TRUE))
})

test_that("CBASS dendrogram plot works for observations", {
  cbass_fit <- CBASS(presidential_speech)

  ## Main settings work
  expect_no_error(plot(cbass_fit, type = "var.dendrogram", k.obs = 3))
  expect_equal(plot(cbass_fit, type = "var.dendrogram", k.obs = 3), invisible(cbass_fit))

  ## Must give at most one of `percent`, `k.obs`, `k.var`
  expect_error(plot(cbass_fit, type = "var.dendrogram", percent = 0.5, k.obs = 3))
  expect_error(plot(cbass_fit, type = "var.dendrogram", percent = 0.5, k.var = 3))
  expect_error(plot(cbass_fit, type = "var.dendrogram", k.obs = 3, k.var = 3))

  ## Error checking on `percent`, `k.obs`, `k.var`
  expect_error(plot(cbass_fit, type = "var.dendrogram", percent = 1.5))
  expect_error(plot(cbass_fit, type = "var.dendrogram", percent = -0.5))
  expect_error(plot(cbass_fit, type = "var.dendrogram", percent = NA))
  expect_error(plot(cbass_fit, type = "var.dendrogram", percent = c(0.25, 0.75)))

  expect_error(plot(cbass_fit, type = "var.dendrogram", k.obs = 3.5))
  expect_error(plot(cbass_fit, type = "var.dendrogram", k.obs = 0))
  expect_error(plot(cbass_fit, type = "var.dendrogram", k.obs = -1))
  expect_error(plot(cbass_fit, type = "var.dendrogram", k.obs = NROW(presidential_speech) + 1))

  expect_error(plot(cbass_fit, type = "var.dendrogram", k.var = 3.5))
  expect_error(plot(cbass_fit, type = "var.dendrogram", k.var = 0))
  expect_error(plot(cbass_fit, type = "var.dendrogram", k.var = -1))
  expect_error(plot(cbass_fit, type = "var.dendrogram", k.var = NCOL(presidential_speech) + 1))

  ## Error checking on `show_clusters`
  expect_error(plot(cbass_fit, type = "var.dendrogram", show_clusters = NA))
  expect_error(plot(cbass_fit, type = "var.dendrogram", show_clusters = c(TRUE, FALSE)))

  ## Error checking on two specially handled arguments
  expect_error(plot(cbass_fit, type = "var.dendrogram", k.obs = 3, dend.branch.width = 0))
  expect_error(plot(cbass_fit, type = "var.dendrogram", k.obs = 3, dend.branch.width = -2))
  expect_error(plot(cbass_fit, type = "var.dendrogram", k.obs = 3, dend.labels.cex   = 0))
  expect_error(plot(cbass_fit, type = "var.dendrogram", k.obs = 3, dend.labels.cex   = -2))

  ## No cluster specification required if `show_clusters` is FALSE
  expect_no_error(plot(cbass_fit, type = "var.dendrogram", show_clusters = FALSE))

  ## Error if `show_clusters = TRUE` but no regularization given
  expect_error(plot(cbass_fit, type = "var.dendrogram", show_clusters = TRUE))
})
