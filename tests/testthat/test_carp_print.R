context("test print.CARP")

test_that("print.CARP works for default settings", {
    carp_fit <- CARP(presidential_speech)
    carp_print <- capture_print(carp_fit)

    expect_str_contains(carp_print, "Algorithm:[ ]+CARP")
    expect_str_contains(carp_print, "Number of Observations:[ ]+44")
    expect_str_contains(carp_print, "Number of Variables:[ ]+75")
    expect_str_contains(carp_print, "Columnwise centering:[ ]+TRUE")
    expect_str_contains(carp_print, "Columnwise scaling:[ ]+FALSE")
    expect_str_contains(carp_print, "Source: Radial Basis Function Kernel Weights")
    expect_str_contains(carp_print, "Distance Metric: Euclidean")
    expect_str_contains(carp_print, stringr::fixed("Scale parameter (phi): 0.01 [Data-Driven]"))
    expect_str_contains(carp_print, stringr::fixed("Sparsified: 4 Nearest Neighbors [Data-Driven]"))
})

test_that("print.CARP works for other algorithms", {
  expect_str_contains(capture_print(CARP(presidential_speech, back_track = TRUE)),
                      "Algorithm:[ ]+CARP-VIZ")

  expect_str_contains(capture_print(CARP(presidential_speech, back_track = TRUE, norm = 1)),
                      "Algorithm:[ ]+CARP-VIZ \\[Back-Tracking Fusion Search\\] \\[L1\\]")

  expect_str_contains(capture_print(CARP(presidential_speech, back_track = FALSE, t = 1.5)),
                      "Algorithm:[ ]+CARP \\(t = 1.5\\)")

  expect_str_contains(capture_print(CARP(presidential_speech, back_track = FALSE, t = 1.5, norm = 1)),
                      "Algorithm:[ ]+CARP \\(t = 1.5\\) \\[L1\\]")
})
