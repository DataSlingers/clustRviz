context("test print.CBASS")

test_that("print.CBASS works for default settings", {
  cbass_fit <- CBASS(presidential_speech)
  cbass_print <- capture_print(cbass_fit)

  expect_str_contains(cbass_print, "Algorithm:[ ]+CBASS-VIZ")
  expect_str_contains(cbass_print, "Number of Observations:[ ]+44")
  expect_str_contains(cbass_print, "Number of Variables:[ ]+75")
  expect_str_contains(cbass_print, "Columnwise centering:[ ]+TRUE")
  expect_str_contains(cbass_print, "Columnwise scaling:[ ]+FALSE")
  expect_str_contains(cbass_print, "Source: Radial Basis Function Kernel Weights")
  expect_str_contains(cbass_print, "Distance Metric: Euclidean")
  expect_str_contains(cbass_print, stringr::fixed("Scale parameter (phi): 0.01 [Data-Driven]"))
  expect_str_contains(cbass_print, stringr::fixed("Sparsified: 4 Nearest Neighbors [Data-Driven]"))
})

test_that("print.CBASS works for other algorithms", {
  expect_str_contains(capture_print(CBASS(presidential_speech, alg.type="cbassviz")),
                      "Algorithm:[ ]+CBASS-VIZ")

  expect_str_contains(capture_print(CBASS(presidential_speech, alg.type="cbassvizl1")),
                      "Algorithm:[ ]+CBASS-VIZ \\[L1\\]")

  expect_str_contains(capture_print(CBASS(presidential_speech, alg.type="cbass", t=1.5)),
                      "Algorithm:[ ]+CBASS \\(t = 1.5\\)")

  expect_str_contains(capture_print(CBASS(presidential_speech, alg.type="cbassl1", t=1.5)),
                      "Algorithm:[ ]+CBASS \\(t = 1.5\\) \\[L1\\]")
})
