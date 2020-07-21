context("CBASS() Missing Data Imputation")

test_that("CBASS Imputes Missing Data", {
  test_dat <- matrix(c(1, 2, 3,
                       4, NA, 6,
                       7, 8, 9), byrow = TRUE, ncol = 3)

  cbass_fit <- CBASS(test_dat)

  # We should impute 5 for the missing value regardless of refit
  expect_equal(get_clustered_data(cbass_fit, k.row = 1)[2,2], 5)
  expect_equal(get_clustered_data(cbass_fit, k.row = 1, refit = FALSE)[2,2], 5)

  # If we want to refit, can't get an estimated centroid from all NA data
  expect_true(is.nan(get_clustered_data(cbass_fit, k.row = 3, refit = TRUE)[2,2]))
  # But we can impute if we don't refit (note that it's a pretty bad imputation...)
  expect_false(is.nan(get_clustered_data(cbass_fit, k.row = 3, refit = FALSE)[2,2]))

  # We should get an error when we can't succesfully impute
  expect_error(CBASS(presidential_speech * NA))
})
