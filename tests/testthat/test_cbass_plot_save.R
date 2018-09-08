context("Test saveviz.CBASS")

## FIXME - Make these rigorous tests: right now, we're just checking for no errors on the default path

test_that("saveviz.CBASS can save a dynamic heatmap", {
  skip_on_cran()

  cbass_fit <- CBASS(presidential_speech)
  expect_no_error(saveviz(cbass_fit, file.name = tempfile(), plot.type = "heatmap", image.type = "dynamic"))
})

test_that("saveviz.CBASS can save a dynamic observation dendrogram", {
  skip_on_cran()

  cbass_fit <- CBASS(presidential_speech)
  expect_no_error(saveviz(cbass_fit, file.name = tempfile(), plot.type = "obs.dendrogram", image.type = "dynamic"))
})

test_that("saveviz.CBASS can save a dynamic variable dendrogram", {
  skip_on_cran()

  cbass_fit <- CBASS(presidential_speech)
  expect_no_error(saveviz(cbass_fit, file.name = tempfile(), plot.type = "var.dendrogram", image.type = "dynamic"))
})

test_that("saveviz.CBASS can save a static (fixed k.obs) heatmap", {
  skip_on_cran()

  cbass_fit <- CBASS(presidential_speech)
  expect_no_error(saveviz(cbass_fit, file.name = tempfile(), plot.type = "heatmap", image.type = "static", k.obs = 5))
})

test_that("saveviz.CBASS can save a static (fixed k.obs) observation dendrogram", {
  skip_on_cran()

  cbass_fit <- CBASS(presidential_speech)
  expect_no_error(saveviz(cbass_fit, file.name = tempfile(), plot.type = "obs.dendrogram", image.type = "static", k.obs = 5))
})

test_that("saveviz.CBASS can save a static (fixed k.obs) variable dendrogram", {
  skip_on_cran()

  cbass_fit <- CBASS(presidential_speech)
  expect_no_error(saveviz(cbass_fit, file.name = tempfile(), plot.type = "var.dendrogram", image.type = "static", k.obs = 5))
})

test_that("saveviz.CBASS can save a static (fixed k.var) heatmap", {
  skip_on_cran()

  cbass_fit <- CBASS(presidential_speech)
  expect_no_error(saveviz(cbass_fit, file.name = tempfile(), plot.type = "heatmap", image.type = "static", k.var = 5))
})

test_that("saveviz.CBASS can save a static (fixed k.var) observation dendrogram", {
  skip_on_cran()

  cbass_fit <- CBASS(presidential_speech)
  expect_no_error(saveviz(cbass_fit, file.name = tempfile(), plot.type = "obs.dendrogram", image.type = "static", k.var = 5))
})

test_that("saveviz.CBASS can save a static (fixed k.var) variable dendrogram", {
  skip_on_cran()

  cbass_fit <- CBASS(presidential_speech)
  expect_no_error(saveviz(cbass_fit, file.name = tempfile(), plot.type = "var.dendrogram", image.type = "static", k.var = 5))
})

test_that("saveviz.CBASS can save a static (fixed percent) heatmap", {
  skip_on_cran()

  cbass_fit <- CBASS(presidential_speech)
  expect_no_error(saveviz(cbass_fit, file.name = tempfile(), plot.type = "heatmap", image.type = "static", percent = 0.5))
})

test_that("saveviz.CBASS can save a static (fixed percent) observation dendrogram", {
  skip_on_cran()

  cbass_fit <- CBASS(presidential_speech)
  expect_no_error(saveviz(cbass_fit, file.name = tempfile(), plot.type = "obs.dendrogram", image.type = "static", percent = 0.5))
})

test_that("saveviz.CBASS can save a static (fixed percent) variable dendrogram", {
  skip_on_cran()

  cbass_fit <- CBASS(presidential_speech)
  expect_no_error(saveviz(cbass_fit, file.name = tempfile(), plot.type = "var.dendrogram", image.type = "static", percent = 0.5))
})
