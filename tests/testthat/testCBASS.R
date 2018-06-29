library(clustRviz)

context("CBASS Function")
data("presidential_speech")
Xdat <- presidential_speech[1:10,1:4]

test_that("CBASS Runs",{
  expect_error(
    cbass.fit <- CBASS(X=Xdat),
    NA
  )
})
cbass.fit <- CBASS(X=Xdat)

test_that("print.CBASS Runs",{
  expect_error(
    print(cbass.fit),
    NA
  )
})

test_that("plot.CBASS Runs",{
  expect_error(
    plot(cbass.fit),
    NA
  )
})

test_that("plot.CBASS obs.dendrogram Runs",{
  expect_error(
    plot(cbass.fit,type='obs.dendrogram'),
    NA
  )
})

test_that("plot.CBASS var.dendrogram Runs",{
  expect_error(
    plot(cbass.fit,type='var.dendrogram'),
    NA
  )
})

test_that("plot.CBASS heatmap Runs",{
  expect_error(
    plot(cbass.fit,type='heatmap'),
    NA
  )
})

test_that("plot.CBASS interactive Runs",{
  expect_error(
    plot(cbass.fit,type='interactive'),
    NA
  )
})

test_that("Clustering CBASS  Runs",{
  expect_error(
    Clustering(cbass.fit,k.obs=1),
    NA
  )
})
