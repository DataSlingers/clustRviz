context("CARP() Missing Data Imputation")

test_that("CARP Imputes Missing Data", {
  test_dat <- matrix(c(NA, 1, 2,
                       3,  1, 2,
                       100, 100, 100), byrow = TRUE, ncol = 3)

  carp_fit <- suppressWarnings(CARP(test_dat)) ## Warns that missForest can't impute well...

  # We should impute 3 for the missing value when we refit
  expect_equal(get_clustered_data(carp_fit, k = 2)[1,1], 3, check.attributes = FALSE)

  # Check we impute correctly in "U-space" as well.
  expect_equal(get_clustered_data(carp_fit, k = 2, refit = FALSE)[1,1],
               get_clustered_data(carp_fit, k = 2, refit = FALSE)[2,1],
               check.attributes = FALSE)

  # If we want to refit, can't get an estimated centroid from all NA data
  expect_true(is.nan(get_clustered_data(carp_fit, k = 3, refit = TRUE)[1,1]))
  # But we can impute if we don't refit (note that it's a pretty bad imputation...)
  expect_false(is.nan(get_clustered_data(carp_fit, k = 3, refit = FALSE)[1,1]))

  # We should get an error when imputation fails to remove all NAs
  expect_error(CARP(test_dat, impute_func = identity))

  # We should also get an error when imputation puts in non-finite or NaN values
  expect_error(CARP(test_dat, impute_func = function(X){X[is.na(X)] <- NaN}))
  expect_error(CARP(test_dat, impute_func = function(X){X[is.na(X)] <- Inf}))

  # We should get an error when we can't impute (too many NAs)
  test_dat[,1] <- NA
  expect_error(CARP(test_dat, impute_func = function(X){X[,1] <- mean(X[,1], na.rm = TRUE); X}))
})
