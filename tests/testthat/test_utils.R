context("Test Utility Functions")

capitalize_string <- function(x){
  x <- gsub("_", " ", x)
  vapply(strsplit(x, " "),
         function(x) paste(paste0(toupper(substring(x, 1, 1)), substring(x, 2)), collapse = " "),
         character(1))
}

test_that("Validators work", {
  is_logical_scalar <- clustRviz:::is_logical_scalar
  is_numeric_scalar <- clustRviz:::is_numeric_scalar
  is_integer_scalar <- clustRviz:::is_integer_scalar
  is_percent_scalar <- clustRviz:::is_percent_scalar

  expect_true(is_logical_scalar(TRUE))
  expect_true(is_logical_scalar(FALSE))
  expect_false(is_logical_scalar(NA))
  expect_false(is_logical_scalar(0))
  expect_false(is_logical_scalar("a"))
  expect_false(is_logical_scalar(c(TRUE, TRUE)))

  expect_true(is_numeric_scalar(3))
  expect_true(is_numeric_scalar(3.5))
  expect_true(is_numeric_scalar(0))
  expect_true(is_numeric_scalar(-4))
  expect_false(is_numeric_scalar(NA))
  expect_false(is_numeric_scalar(c(2, 5)))
  expect_false(is_numeric_scalar("a"))

  expect_true(is_integer_scalar(3))
  expect_false(is_integer_scalar(3.5))
  expect_true(is_integer_scalar(0))
  expect_true(is_integer_scalar(-4))
  expect_false(is_integer_scalar(NA))
  expect_false(is_integer_scalar(c(2, 5)))
  expect_false(is_integer_scalar("a"))

  expect_true(is_percent_scalar(0.3))
  expect_true(is_percent_scalar(1))
  expect_true(is_percent_scalar(0))
  expect_false(is_percent_scalar(1.5))
  expect_false(is_percent_scalar(-1.5))
  expect_false(is_percent_scalar(NA))
  expect_false(is_percent_scalar(c(0.2, 0.5)))
  expect_false(is_percent_scalar("a"))

  is_square <- clustRviz:::is_square
  expect_true(is_square(matrix(1, 5, 5)))
  expect_false(is_square(matrix(1, 5, 3)))
})

test_that("Capitalization works", {
  capitalize_string <- clustRviz:::capitalize_string

  expect_equal("A", capitalize_string("a"))
  expect_equal("Abc", capitalize_string("abc"))
  expect_equal("ABc", capitalize_string("ABc"))
  expect_equal("A Fantastic Cow", capitalize_string("a fantastic cow"))
})
