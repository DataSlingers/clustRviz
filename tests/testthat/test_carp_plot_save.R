context("Test saveviz.CARP")

## FIXME - Make these rigorous tests: right now, we're just checking for no errors on the default path

test_that("saveviz.CARP can save a dynamic path", {
  skip("Need to rework dynamic visualizations")

  carp_fit <- CARP(presidential_speech)
  expect_no_error(saveviz(carp_fit, file.name = tempfile(), plot.type = "path", image.type = "dynamic"))
})

test_that("saveviz.CARP can save a dynamic dendrogram", {
  skip("Need to rework dynamic visualizations")

  carp_fit <- CARP(presidential_speech)
  expect_no_error(saveviz(carp_fit, file.name = tempfile(), plot.type = "dendrogram", image.type = "dynamic"))
})


test_that("saveviz.CARP can save a static path as a PNG", {
  skip_on_cran()

  carp_fit <- CARP(presidential_speech)
  temp_file <- file.path(tempdir(), "tester.png")
  expect_equal(invisible(temp_file), saveviz(carp_fit, file.name = temp_file, type = "path", image.type = "static"))
  expect_true(file.exists(temp_file))
})

test_that("saveviz.CARP can save a static dendrogram as a JPG", {
  skip_on_cran()

  carp_fit <- CARP(presidential_speech)
  temp_file <- file.path(tempdir(), "tester.jpg")
  expect_equal(invisible(temp_file), saveviz(carp_fit, file.name = temp_file, type = "dendrogram", image.type = "static"))
  expect_true(file.exists(temp_file))
})
