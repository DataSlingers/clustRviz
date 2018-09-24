context("Logging Tests")

test_that("Logging controls work", {
    expect_error(clustRviz_logger_level("BAD LEVEL"))

    clustRviz_logger_level("INFO")
    expect_equal("INFO", clustRviz_logger_level())

    clustRviz_logger_level("MESSAGE")
})

test_that("INFO and DEBUG message print as expected", {
    clustRviz_logger_level("MESSAGE")

    expect_silent(clustRviz:::crv_info("A message"))
    expect_silent(clustRviz:::crv_debug("A message"))

    clustRviz_logger_level("DEBUG")

    expect_output(clustRviz:::crv_info("A message"), "[INFO]")
    expect_output(clustRviz:::crv_info("A message"), "A message")
    expect_output(clustRviz:::crv_debug("The message"), "[DEBUG]")
    expect_output(clustRviz:::crv_debug("The message"), "The message")

    clustRviz_logger_level("MESSAGE")
})

test_that("Supressing messages works", {
    # At INFO level, everything is shown in R
    clustRviz_logger_level("INFO")

    expect_error(clustRviz:::crv_error("ERROR"))
    expect_warning(clustRviz:::crv_warning("WARNING"))
    expect_message(clustRviz:::crv_message("MESSAGE"))

    # At MESSAGE level, everything is shown in R
    clustRviz_logger_level("MESSAGE")

    expect_error(clustRviz:::crv_error("ERROR"))
    expect_warning(clustRviz:::crv_warning("WARNING"))
    expect_message(clustRviz:::crv_message("MESSAGE"))

    # At WARNING level, we don't get a message
    clustRviz_logger_level("WARNING")

    expect_error(clustRviz:::crv_error("ERROR"))
    expect_warning(clustRviz:::crv_warning("WARNING"))
    expect_no_message(clustRviz:::crv_message("MESSAGE"))

    # At ERROR level, we don't get a message or warning
    clustRviz_logger_level("ERROR")

    expect_error(clustRviz:::crv_error("ERROR"))
    expect_no_warning(clustRviz:::crv_warning("WARNING"))
    expect_no_message(clustRviz:::crv_message("MESSAGE"))

    clustRviz_logger_level("MESSAGE")
})

test_that("No extra newlines", {
    clustRviz_logger_level("DEBUG")

    e <- tryCatch(clustRviz:::crv_error("MY ERROR"), error=identity)
    expect_equal(str_count(e$message, "\n"), 1)

    e <- tryCatch(clustRviz:::crv_warning("MY WARNING"), warning=identity)
    expect_equal(str_count(e$message, "\n"), 1)

    e <- tryCatch(clustRviz:::crv_message("MY MESSAGE"), message=identity)
    expect_equal(str_count(e$message, "\n"), 1)

    e <- tryCatch(clustRviz:::crv_error("MY ERROR\nON TWO LINES"), error=identity)
    expect_equal(str_count(e$message, "\n"), 2)

    clustRviz_logger_level("MESSAGE")
})

test_that("Function capture works at R level", {
    clustRviz_logger_level("MESSAGE")

    f <- function(x){clustRviz:::crv_error("ERROR MESSAGE")}

    e <- tryCatch(f(), error=identity)

    expect_str_contains(e$message, "ERROR MESSAGE")
    expect_str_contains(e$message, "(Called from f)")
    expect_true(is.null(e$call))
    expect_true(is.null(e$cppstack))

    f <- function(x){clustRviz:::crv_error("ERROR MESSAGE", call=FALSE)}
    e <- tryCatch(f(), error=identity)

    expect_false(grepl("\\(Called from f\\)", e$message))

    f <- function(x){clustRviz:::crv_error("ERROR MESSAGE", call="my func")}
    e <- tryCatch(f(), error=identity)

    expect_true(grepl("\\(Called from my func\\)", e$message))

    f <- function(x){clustRviz:::crv_warning("WARNING MESSAGE", call=FALSE)}
    e <- tryCatch(f(), warning=identity)

    expect_false(grepl("\\(Called from f\\)", e$message))

    f <- function(x){clustRviz:::crv_warning("WARNING MESSAGE", call="my func")}
    e <- tryCatch(f(), warning=identity)

    expect_true(grepl("\\(Called from my func\\)", e$message))
})
