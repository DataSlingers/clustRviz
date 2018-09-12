library(clustRviz)

context("CARP Function")
data("presidential_speech")
Xdat <- presidential_speech[1:10,1:4]
nobs <- nrow(Xdat)




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


context("CARP fitting w/ uniform weights")
test_that("CARP Runs w/ carpviz",{
  expect_error(
    carp.fit <- CARP(
        X=Xdat,
        control = list(alg.type='carpviz'),
        weights = rep(1,times=choose(nobs,2))),
    NA
  )
})
test_that("CARP Runs w/ carpvizl1",{
  expect_error(
    carp.fit <- CARP(
        X=Xdat,
        control = list(alg.type='carpvizl1'),
        weights = rep(1,times=choose(nobs,2))),
    NA
  )
})
test_that("CARP Runs w/ carp 1.05",{
  expect_error(
    carp.fit <- CARP(
        X=Xdat,
        control = list(alg.type='carp',t=1.05),
        weights = rep(1,times=choose(nobs,2))),
    NA
  )
})
test_that("CARP Runs w/ carp 1.1",{
  expect_error(
    carp.fit <- CARP(
        X=Xdat,
        control = list(alg.type='carp',t=1.1),
        weights = rep(1,times=choose(nobs,2))),
    NA
  )
})
test_that("CARP Runs w/ carpl1 1.05",{
  expect_error(
    carp.fit <- CARP(
        X=Xdat,
        control = list(alg.type='carpl1',t=1.05),
        weights = rep(1,times=choose(nobs,2))),
    NA
  )
})
test_that("CARP Runs w/ carpl1 1.1",{
  expect_error(
    carp.fit <- CARP(
        X=Xdat,
        control = list(alg.type='carpl1',t=1.1),
        weights = rep(1,times=choose(nobs,2))),
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



context("CARP clustering")
test_that("clustering.CARP w/ k Runs",{
  expect_error(
    clustering(carp.fit,k=4),
    NA
  )
})
test_that("clustering.CARP w/ percent Runs",{
  expect_error(
    clustering(carp.fit,percent=.5),
    NA
  )
})
test_that("clustering.CARP w/o args Runs",{
  expect_error(
    clustering(carp.fit),
    NA
  )
})





# context("CARP Saving")
# test_that("saveviz Runs w dynamic path",{
#   expect_error(
#     saveviz(carp.fit,file.name = 'carp_path_dynamic.png',plot.type = 'path',image.type = 'dynamic'),
#     NA
#   )
# })
# test_that("saveviz Runs w dynamic dendrogram",{
#   expect_error(
#     saveviz(carp.fit,file.name = 'carp_dend_dynamic.png',plot.type = 'dendrogram',image.type = 'dynamic'),
#     NA
#   )
# })
#
# test_that("saveviz Runs w static path, percent",{
#   expect_error(
#     saveviz(carp.fit,file.name = 'carp_path_static.png',plot.type = 'path',image.type = 'static',percent = .1),
#     NA
#   )
# })
# test_that("saveviz Runs w static path, k",{
#   expect_error(
#     saveviz(carp.fit,file.name = 'carp_path_static.png',plot.type = 'path',image.type = 'static',k=4),
#     NA
#   )
# })
# test_that("saveviz Runs w static dendrogram, percent",{
#   expect_error(
#     saveviz(carp.fit,file.name = 'carp_dend_static.png',plot.type = 'dendrogram',image.type = 'static',percent = .1),
#     NA
#   )
# })
# test_that("saveviz Runs w static dendrogram, k",{
#   expect_error(
#     saveviz(carp.fit,file.name = 'carp_dend_static.png',plot.type = 'dendrogram',image.type = 'static',k=4),
#     NA
#   )
# })
