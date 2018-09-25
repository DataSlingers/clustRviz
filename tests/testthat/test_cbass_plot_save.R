context("Test saveviz.CBASS")

## FIXME - Make these rigorous tests: right now, we're just checking for no errors on the default path

test_that("saveviz.CBASS can save a dynamic heatmap", {
  skip("Need to rework dynamic visualizations")

  cbass_fit <- CBASS(presidential_speech)
  expect_no_error(saveviz(cbass_fit, file.name = tempfile(), plot.type = "heatmap", dynamic = TRUE))
})

test_that("saveviz.CBASS can save a dynamic observation dendrogram", {
  skip("Need to rework dynamic visualizations")

  cbass_fit <- CBASS(presidential_speech)
  expect_no_error(saveviz(cbass_fit, file.name = tempfile(), plot.type = "obs.dendrogram", dynamic = TRUE))
})

test_that("saveviz.CBASS can save a dynamic variable dendrogram", {
  skip("Need to rework dynamic visualizations")

  cbass_fit <- CBASS(presidential_speech)
  expect_no_error(saveviz(cbass_fit, file.name = tempfile(), plot.type = "var.dendrogram", dynamic = TRUE))
})
test_that("saveviz.CBASS can save a static heatmap as a PNG", {
  skip_on_cran()

  cbass_fit <- CBASS(presidential_speech)
  temp_file <- file.path(tempdir(), "tester.png")
  expect_equal(invisible(temp_file), saveviz(cbass_fit, file.name = temp_file, type = "heatmap", dynamic = FALSE))
  expect_true(file.exists(temp_file))
})

test_that("saveviz.CBASS can save a static observation dendrogram as a JPG", {
  skip_on_cran()

  cbass_fit <- CBASS(presidential_speech)
  temp_file <- file.path(tempdir(), "tester.jpg")
  expect_equal(invisible(temp_file), saveviz(cbass_fit, file.name = temp_file, type = "obs.dendrogram", dynamic = FALSE))
  expect_true(file.exists(temp_file))
})


test_that("saveviz.CBASS can save a static variable dendrogram as a PDF", {
  skip_on_cran()

  cbass_fit <- CBASS(presidential_speech)
  temp_file <- file.path(tempdir(), "tester.pdf")
  expect_equal(invisible(temp_file), saveviz(cbass_fit, file.name = temp_file, type = "var.dendrogram", dynamic = FALSE))
  expect_true(file.exists(temp_file))
})

test_that("saveviz.CBASS can save a dynamic observation dendrogram as a GIF", {
  skip_on_cran()

  cbass_fit <- CBASS(presidential_speech)
  temp_file <- file.path(tempdir(), "tester.gif")
  expect_no_warning(sv_res <- saveviz(cbass_fit, file.name = temp_file, type = "obs.dendrogram", dynamic = TRUE))
  expect_equal(sv_res, invisible(temp_file))
  expect_true(file.exists(temp_file))
})

test_that("saveviz.CBASS can save a dynamic variable dendrogram as a GIF", {
  skip_on_cran()

  cbass_fit <- CBASS(presidential_speech)
  temp_file <- file.path(tempdir(), "tester.gif")
  expect_no_warning(sv_res <- saveviz(cbass_fit, file.name = temp_file, type = "var.dendrogram", dynamic = TRUE))
  expect_equal(sv_res, invisible(temp_file))
  expect_true(file.exists(temp_file))
})

test_that("saveviz.CBASS can save a dynamic heatmap as a GIF", {
  skip_on_cran()

  cbass_fit <- CBASS(presidential_speech)
  temp_file <- file.path(tempdir(), "tester.gif")
  expect_no_warning(sv_res <- saveviz(cbass_fit, file.name = temp_file, type = "heatmap", dynamic = TRUE))
  expect_equal(sv_res, invisible(temp_file))
  expect_true(file.exists(temp_file))
})
