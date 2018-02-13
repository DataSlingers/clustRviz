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

carp.fit <- CARP(X=Xdat,
                 obs.labels=presidential_speech$labels[1:10])
test_that("print.CARP Runs",{
  expect_error(
    print(carp.fit),
    NA
  )
})

test_that("Clustering.CARP w/ k Runs",{
  expect_error(
    Clustering(carp.fit,k=4),
    NA
  )
})

test_that("Clustering.CARP w/ percent Runs",{
  expect_error(
    Clustering(carp.fit,percent=.5),
    NA
  )
})

test_that("Clustering.CARP w/o args Runs",{
  expect_error(
    Clustering(carp.fit),
    NA
  )
})
