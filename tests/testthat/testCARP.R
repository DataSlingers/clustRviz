library(clustRviz)

context("CARP Function")
data("presidential_speech")
Xdat <- presidential_speech[1:10,1:4]




context("CARP fitting")
test_that("CARP Runs",{
  expect_error(
    carp.fit <- CARP(
        X=Xdat),
    NA
  )
})
test_that("CARP Runs w/ carpviz",{
  expect_error(
    carp.fit <- CARP(
        X=Xdat,
        control = list(alg.type='carpviz')),
    NA
  )
})
test_that("CARP Runs w/ carpvizl1",{
  expect_error(
    carp.fit <- CARP(
        X=Xdat,
        control = list(alg.type='carpvizl1')),
    NA
  )
})
test_that("CARP Runs w/ carp 1.05",{
  expect_error(
    carp.fit <- CARP(
        X=Xdat,
        control = list(alg.type='carp',t=1.05)),
    NA
  )
})
test_that("CARP Runs w/ carp 1.1",{
  expect_error(
    carp.fit <- CARP(
        X=Xdat,
        control = list(alg.type='carp',t=1.1)),
    NA
  )
})
test_that("CARP Runs w/ carpl1 1.05",{
  expect_error(
    carp.fit <- CARP(
        X=Xdat,
        control = list(alg.type='carpl1',t=1.05)),
    NA
  )
})
test_that("CARP Runs w/ carpl1 1.1",{
  expect_error(
    carp.fit <- CARP(
        X=Xdat,
        control = list(alg.type='carpl1',t=1.1)),
    NA
  )
})






context("CARP printing")
carp.fit <- CARP(X=Xdat)
carp.fitl1 <- CARP(X=Xdat,control = list(alg.type='carpvizl1'))
test_that("print.CARP Runs",{
  expect_error(
    print(carp.fit),
    NA
  )
})



context("CARP Clustering")
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





context("CARP Saving")
test_that("saveviz Runs w/o args",{
  expect_error(
    saveviz(carp.fit,file.name = 'saveviz1.png'),
    NA
  )
})
test_that("saveviz Runs w dynamic path",{
  expect_error(
    saveviz(carp.fit,file.name = 'saveviz2.png',plot.type = 'path',image.type = 'dynamic'),
    NA
  )
})
test_that("saveviz Runs w dynamic dendrogram",{
  expect_error(
    saveviz(carp.fit,file.name = 'saveviz3.png',plot.type = 'dendrogram',image.type = 'dynamic'),
    NA
  )
})
test_that("saveviz Runs w static path",{
  expect_error(
    saveviz(carp.fit,file.name = 'saveviz4.png',plot.type = 'path',image.type = 'static'),
    NA
  )
})
test_that("saveviz Runs w static path",{
  expect_error(
    saveviz(carp.fitl1,file.name = 'saveviz45.png',plot.type = 'path',image.type = 'static'),
    NA
  )
})
test_that("saveviz Runs w static dendrogram",{
  expect_error(
    saveviz(carp.fit,file.name = 'saveviz5.png',plot.type = 'dendrogram',image.type = 'static'),
    NA
  )
})
