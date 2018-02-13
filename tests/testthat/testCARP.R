library(clustRviz)

context("CARP Function")
data("presidential_speech")
Xdat <- presidential_speech$X[1:10,1:4]

test_that("CARP Runs",{
  expect_error(
    carp.fit <- CARP(
        X=Xdat,
        obs.labels=presidential_speech$labels[1:10]),
    NA
  )
})
