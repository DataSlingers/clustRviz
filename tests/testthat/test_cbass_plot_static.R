context("Test Static CBASS Plots")

## Add tests for static CBASS plots
##
## Right now these are just "smoke" tests (i.e., runs in the default way without error)
## but we will add actual tests later (see GH #44)

test_that("CBASS path plot works for rows", {
  cbass_fit <- CBASS(presidential_speech)

  ## Main settings work
  expect_no_error(plot(cbass_fit, type = "row.path"))
  expect_is(plot(cbass_fit, type = "row.path"), "gg")
  expect_no_error(plot(cbass_fit, type = "row.path", k.row  = 3, axis = c("amount", "appropri")))
  expect_no_error(plot(cbass_fit, type = "row.path", k.row  = 3, axis = c("PC1", "PC3")))

  ## Must give at most one of `percent`, `k.row`, `k.col`
  expect_error(plot(cbass_fit, type = "row.path", percent = 0.5, k.row = 3))
  expect_error(plot(cbass_fit, type = "row.path", percent = 0.5, k.col = 3))
  expect_error(plot(cbass_fit, type = "row.path", k.row = 3, k.col = 3))

  ## Error checking on `percent` and `k`
  expect_error(plot(cbass_fit, type = "row.path", percent = 1.5))
  expect_error(plot(cbass_fit, type = "row.path", percent = -0.5))
  expect_error(plot(cbass_fit, type = "row.path", percent = NA))
  expect_error(plot(cbass_fit, type = "row.path", percent = c(0.25, 0.75)))

  expect_error(plot(cbass_fit, type = "row.path", k.row = 3.5))
  expect_error(plot(cbass_fit, type = "row.path", k.row = 0))
  expect_error(plot(cbass_fit, type = "row.path", k.row = -1))
  expect_error(plot(cbass_fit, type = "row.path", k.row = NROW(presidential_speech) + 1))

  expect_error(plot(cbass_fit, type = "row.path", k.col = 3.5))
  expect_error(plot(cbass_fit, type = "row.path", k.col = 0))
  expect_error(plot(cbass_fit, type = "row.path", k.col = -1))
  expect_error(plot(cbass_fit, type = "row.path", k.col = NCOL(presidential_speech) + 1))

  ## Error on unknown arguments
  expect_error(plot(cbass_fit, type = "row.path", 5))
  expect_error(plot(cbass_fit, type = "row.path", a = 5))
})

test_that("CBASS path plot works for columns", {
  cbass_fit <- CBASS(presidential_speech)

  ## Main settings work
  expect_no_error(plot(cbass_fit, type = "col.path", k.row = 3))
  expect_is(plot(cbass_fit, type = "col.path", k.row  = 3), "gg")
  expect_no_error(plot(cbass_fit, type = "col.path", k.row  = 3, axis = c("Barack Obama", "George Washington")))
  expect_no_error(plot(cbass_fit, type = "col.path", k.row  = 3, axis = c("PC1", "PC3")))

  ## Must give at most one of `percent`, `k.row`, `k.col`
  expect_error(plot(cbass_fit, type = "col.path", percent = 0.5, k.row = 3))
  expect_error(plot(cbass_fit, type = "col.path", percent = 0.5, k.col = 3))
  expect_error(plot(cbass_fit, type = "col.path", k.row = 3, k.col = 3))

  ## Error checking on `percent` and `k`
  expect_error(plot(cbass_fit, type = "col.path", percent = 1.5))
  expect_error(plot(cbass_fit, type = "col.path", percent = -0.5))
  expect_error(plot(cbass_fit, type = "col.path", percent = NA))
  expect_error(plot(cbass_fit, type = "col.path", percent = c(0.25, 0.75)))

  expect_error(plot(cbass_fit, type = "col.path", k.row = 3.5))
  expect_error(plot(cbass_fit, type = "col.path", k.row = 0))
  expect_error(plot(cbass_fit, type = "col.path", k.row = -1))
  expect_error(plot(cbass_fit, type = "col.path", k.row = NROW(presidential_speech) + 1))

  expect_error(plot(cbass_fit, type = "col.path", k.col = 3.5))
  expect_error(plot(cbass_fit, type = "col.path", k.col = 0))
  expect_error(plot(cbass_fit, type = "col.path", k.col = -1))
  expect_error(plot(cbass_fit, type = "col.path", k.col = NCOL(presidential_speech) + 1))

  ## Error on unknown arguments
  expect_error(plot(cbass_fit, type = "col.path", 5))
  expect_error(plot(cbass_fit, type = "col.path", a = 5))
})

test_that("CBASS dendrogram plot works for rows", {
  cbass_fit <- CBASS(presidential_speech)

  ## Main settings work
  expect_no_error(plot(cbass_fit, type = "row.dendrogram", k.row = 3))
  expect_s3_class(plot(cbass_fit, type = "row.dendrogram", k.row = 3), 'ggplot')

  ## Must give at most one of `percent`, `k.row`, `k.col`
  expect_error(plot(cbass_fit, type = "row.dendrogram", percent = 0.5, k.row = 3))
  expect_error(plot(cbass_fit, type = "row.dendrogram", percent = 0.5, k.col = 3))
  expect_error(plot(cbass_fit, type = "row.dendrogram", k.row = 3, k.col = 3))

  ## Error checking on `percent`, `k.row`, `k.col`
  expect_error(plot(cbass_fit, type = "row.dendrogram", percent = 1.5))
  expect_error(plot(cbass_fit, type = "row.dendrogram", percent = -0.5))
  expect_error(plot(cbass_fit, type = "row.dendrogram", percent = NA))
  expect_error(plot(cbass_fit, type = "row.dendrogram", percent = c(0.25, 0.75)))

  expect_error(plot(cbass_fit, type = "row.dendrogram", k.row = 3.5))
  expect_error(plot(cbass_fit, type = "row.dendrogram", k.row = 0))
  expect_error(plot(cbass_fit, type = "row.dendrogram", k.row = -1))
  expect_error(plot(cbass_fit, type = "row.dendrogram", k.row = NROW(presidential_speech) + 1))

  expect_error(plot(cbass_fit, type = "row.dendrogram", k.col = 3.5))
  expect_error(plot(cbass_fit, type = "row.dendrogram", k.col = 0))
  expect_error(plot(cbass_fit, type = "row.dendrogram", k.col = -1))
  expect_error(plot(cbass_fit, type = "row.dendrogram", k.col = NCOL(presidential_speech) + 1))
})

test_that("CBASS dendrogram plot works for columns", {
  cbass_fit <- CBASS(presidential_speech)

  ## Main settings work
  expect_no_error(plot(cbass_fit, type = "col.dendrogram", k.col = 3))
  expect_s3_class(plot(cbass_fit, type = "col.dendrogram", k.col = 3), 'ggplot')

  ## Must give at most one of `percent`, `k.row`, `k.col`
  expect_error(plot(cbass_fit, type = "col.dendrogram", percent = 0.5, k.row = 3))
  expect_error(plot(cbass_fit, type = "col.dendrogram", percent = 0.5, k.col = 3))
  expect_error(plot(cbass_fit, type = "col.dendrogram", k.row = 3, k.col = 3))

  ## Error checking on `percent`, `k.row`, `k.col`
  expect_error(plot(cbass_fit, type = "col.dendrogram", percent = 1.5))
  expect_error(plot(cbass_fit, type = "col.dendrogram", percent = -0.5))
  expect_error(plot(cbass_fit, type = "col.dendrogram", percent = NA))
  expect_error(plot(cbass_fit, type = "col.dendrogram", percent = c(0.25, 0.75)))

  expect_error(plot(cbass_fit, type = "col.dendrogram", k.row = 3.5))
  expect_error(plot(cbass_fit, type = "col.dendrogram", k.row = 0))
  expect_error(plot(cbass_fit, type = "col.dendrogram", k.row = -1))
  expect_error(plot(cbass_fit, type = "col.dendrogram", k.row = NROW(presidential_speech) + 1))

  expect_error(plot(cbass_fit, type = "col.dendrogram", k.col = 3.5))
  expect_error(plot(cbass_fit, type = "col.dendrogram", k.col = 0))
  expect_error(plot(cbass_fit, type = "col.dendrogram", k.col = -1))
  expect_error(plot(cbass_fit, type = "col.dendrogram", k.col = NCOL(presidential_speech) + 1))
})

test_that("CBASS heatmap plot works", {
  cbass_fit <- CBASS(presidential_speech)

  ## Main settings work
  expect_no_error(plot(cbass_fit, type = "heatmap"))
  expect_s3_class(plot(cbass_fit, type = "heatmap"), 'ggplot')
  expect_s3_class(plot(cbass_fit, type = "heatmap", percent = 0.5), 'ggplot')
  expect_s3_class(plot(cbass_fit, type = "heatmap", k.row = 5), 'ggplot')
  expect_s3_class(plot(cbass_fit, type = "heatmap", k.col = 5), 'ggplot')

  ## Must give at most one of `percent`, `k.row`, `k.col`
  expect_error(plot(cbass_fit, type = "heatmap", percent = 0.5, k.row = 3))
  expect_error(plot(cbass_fit, type = "heatmap", percent = 0.5, k.col = 3))
  expect_error(plot(cbass_fit, type = "heatmap", k.row = 3, k.col = 3))
  expect_error(plot(cbass_fit, type = "heatmap", k.row = 3, k.col = 3, percent = 0.5))

  ## Error checking on `percent`, `k.row`, `k.col`
  expect_error(plot(cbass_fit, type = "heatmap", percent = 1.5))
  expect_error(plot(cbass_fit, type = "heatmap", percent = -0.5))
  expect_error(plot(cbass_fit, type = "heatmap", percent = NA))
  expect_error(plot(cbass_fit, type = "heatmap", percent = c(0.25, 0.75)))

  expect_error(plot(cbass_fit, type = "heatmap", k.row = 3.5))
  expect_error(plot(cbass_fit, type = "heatmap", k.row = 0))
  expect_error(plot(cbass_fit, type = "heatmap", k.row = -1))
  expect_error(plot(cbass_fit, type = "heatmap", k.row = NROW(presidential_speech) + 1))

  expect_error(plot(cbass_fit, type = "heatmap", k.col = 3.5))
  expect_error(plot(cbass_fit, type = "heatmap", k.col = 0))
  expect_error(plot(cbass_fit, type = "heatmap", k.col = -1))
  expect_error(plot(cbass_fit, type = "heatmap", k.col = NCOL(presidential_speech) + 1))
})

test_that("CBASS heatmaply (javascript) plot works", {
  cbass_fit <- CBASS(presidential_speech)

  ## Main settings work
  expect_no_error(plot(cbass_fit, type = "heatmap", interactive = TRUE))
  expect_is(plot(cbass_fit, type = "heatmap", interactive = TRUE), "plotly")
  expect_no_error(plot(cbass_fit, type = "heatmap", interactive = TRUE, k.row = 3))
  expect_no_error(plot(cbass_fit, type = "heatmap", interactive = TRUE, k.col = 3))
  expect_no_error(plot(cbass_fit, type = "heatmap", interactive = TRUE, percent = 0.3))

  ## Must give at most one of `percent`, `k.row`, `k.col`
  expect_error(plot(cbass_fit, type = "heatmap", interactive = TRUE, percent = 0.5, k.row = 3))
  expect_error(plot(cbass_fit, type = "heatmap", interactive = TRUE, percent = 0.5, k.col = 3))
  expect_error(plot(cbass_fit, type = "heatmap", interactive = TRUE, k.row = 3, k.col = 3))

  ## Error checking on `percent` and `k`
  expect_error(plot(cbass_fit, type = "heatmap", interactive = TRUE, percent = 1.5))
  expect_error(plot(cbass_fit, type = "heatmap", interactive = TRUE, percent = -0.5))
  expect_error(plot(cbass_fit, type = "heatmap", interactive = TRUE, percent = NA))
  expect_error(plot(cbass_fit, type = "heatmap", interactive = TRUE, percent = c(0.25, 0.75)))

  expect_error(plot(cbass_fit, type = "heatmap", interactive = TRUE, k.row = 3.5))
  expect_error(plot(cbass_fit, type = "heatmap", interactive = TRUE, k.row = 0))
  expect_error(plot(cbass_fit, type = "heatmap", interactive = TRUE, k.row = -1))
  expect_error(plot(cbass_fit, type = "heatmap", interactive = TRUE, k.row = NROW(presidential_speech) + 1))

  expect_error(plot(cbass_fit, type = "heatmap", interactive = TRUE, k.col = 3.5))
  expect_error(plot(cbass_fit, type = "heatmap", interactive = TRUE, k.col = 0))
  expect_error(plot(cbass_fit, type = "heatmap", interactive = TRUE, k.col = -1))
  expect_error(plot(cbass_fit, type = "heatmap", interactive = TRUE, k.col = NCOL(presidential_speech) + 1))
})
