library(stringr)

expect_no_error <- function(object, ..., all=FALSE, info=NULL, label=NULL){
  expect_error(object, regexp=NA, ..., all=all, info=info, label=label)
}

expect_no_warning <- function(object, ..., all=FALSE, info=NULL, label=NULL){
  expect_warning(object, regexp=NA, ..., all=all, info=info, label=label)
}

expect_no_message <- function(object, ..., all=FALSE, info=NULL, label=NULL){
  expect_message(object, regexp=NA, ..., all=all, info=info, label=label)
}

expect_str_contains <- function(object, expected, info=NULL, label=NULL){
  if(!is.character(object)) object <- as.character(object)
  if(!is.character(expected)) expected <- as.character(expected)

  expect_true(all(str_detect(object, expected)),
              info=info, label=label)
}

expect_zero <- function(object, ..., info=NULL, label=NULL, expected.label=NULL){
  expect_equal(object, 0, ..., info=info, label=label, expected.label=expected.label)
}

expect_zeros <- function(object, ..., info=NULL, label=NULL, expected.label=NULL){
  expect_equal(object, rep(0, length(object)), ..., info=info, label=label, expected.label=expected.label)
}

expect_ones <- function(object, ..., info=NULL, label=NULL, expected.label=NULL){
  expect_equal(object, rep(1, length(object)), ..., info=info, label=label, expected.label=expected.label)
}

capture_print <- function(x, ...){
  paste(capture.output(print(x, ...)), collapse="\n")
}

list_all_equal <- function(x) {
  all(vapply(seq_len(length(x) - 1), function(i) isTRUE(all.equal(x[[i]], x[[i + 1]])), logical(1)))
}

num_unique      <- clustRviz:::num_unique
num_unique_rows <- clustRviz:::num_unique_rows
num_unique_cols <- clustRviz:::num_unique_cols
