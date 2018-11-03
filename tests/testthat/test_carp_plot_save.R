context("Test saveviz.CARP")

## FIXME - Make these rigorous tests: right now, we're just checking for no errors on the default path

test_that("saveviz.CARP can save a static path as a PNG", {
  skip_on_cran()

  carp_fit <- CARP(presidential_speech)
  temp_file <- file.path(tempdir(), "tester.png")
  expect_equal(invisible(temp_file), saveviz(carp_fit, file.name = temp_file, type = "path", dynamic = FALSE))
  expect_true(file.exists(temp_file))
})

test_that("saveviz.CARP can save a static dendrogram as a JPG", {
  skip_on_cran()

  carp_fit <- CARP(presidential_speech)
  temp_file <- file.path(tempdir(), "tester.jpg")
  expect_equal(invisible(temp_file), saveviz(carp_fit, file.name = temp_file, type = "dendrogram", dynamic = FALSE))
  expect_true(file.exists(temp_file))
})

test_that("saveviz.CARP can save a dynamic dendrogram as a GIF", {
  skip_on_cran()

  carp_fit <- CARP(presidential_speech)
  temp_file <- file.path(tempdir(), "tester.gif")
  expect_no_warning(sv_res <- saveviz(carp_fit, file.name = temp_file, type = "dendrogram", dynamic = TRUE))
  expect_equal(sv_res, invisible(temp_file))
  expect_true(file.exists(temp_file))
})

test_that("saveviz.CARP can save a dynamic path as a GIF", {
  skip_on_cran()

  carp_fit <- CARP(presidential_speech)
  temp_file <- file.path(tempdir(), "tester.gif")
  expect_no_warning(sv_res <- saveviz(carp_fit, file.name = temp_file, type = "path", dynamic = TRUE))
  expect_equal(sv_res, invisible(temp_file))
  expect_true(file.exists(temp_file))
})

test_that("saveviz.CARP can export JS to a file", {
  skip_on_cran()

  carp_fit <- CARP(presidential_speech)
  temp_file <- file.path(tempdir(), "tester.html")
  expect_warning(saveviz(carp_fit, file.name = temp_file, type = "js"))
  expect_no_warning(sv_res <- saveviz(carp_fit, file.name = temp_file, type = "js", dynamic = FALSE))
  expect_equal(sv_res, invisible(temp_file))
  expect_true(file.exists(temp_file))
})
