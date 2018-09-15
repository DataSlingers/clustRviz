context("Test saveviz.CARP")

## FIXME - Make these rigorous tests: right now, we're just checking for no errors on the default path

test_that("saveviz.CARP can save a dynamic path", {
  skip_on_cran(); skip_on_travis()

  carp_fit <- CARP(presidential_speech)
  expect_no_error(saveviz(carp_fit, file.name = tempfile(), plot.type = "path", image.type = "dynamic"))
})

test_that("saveviz.CARP can save a dynamic dendrogram", {
  skip_on_cran(); skip_on_travis()

  carp_fit <- CARP(presidential_speech)
  expect_no_error(saveviz(carp_fit, file.name = tempfile(), plot.type = "dendrogram", image.type = "dynamic"))
})


test_that("saveviz.CARP can save a static (fixed k) path", {
  skip_on_cran(); skip_on_travis()

  carp_fit <- CARP(presidential_speech)
  expect_no_error(saveviz(carp_fit, file.name = tempfile(), plot.type = "path", image.type = "static", k = 5))
})

test_that("saveviz.CARP can save a static (fixed k) dendrogram", {
  skip_on_cran(); skip_on_travis()

  carp_fit <- CARP(presidential_speech)
  expect_no_error(saveviz(carp_fit, file.name = tempfile(), plot.type = "dendrogram", image.type = "static", k = 5))
})

test_that("saveviz.CARP can save a static (fixed percent) path", {
  skip_on_cran(); skip_on_travis()

  carp_fit <- CARP(presidential_speech)
  expect_no_error(saveviz(carp_fit, file.name = tempfile(), plot.type = "path", image.type = "static", percent = 0.5))
})

test_that("saveviz.CARP can save a static (fixed percent) dendrogram", {
  skip_on_cran(); skip_on_travis()

  carp_fit <- CARP(presidential_speech)
  expect_no_error(saveviz(carp_fit, file.name = tempfile(), plot.type = "dendrogram", image.type = "static", percent = 0.5))
})
