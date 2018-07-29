library(clustRviz)

context("CBASS Function")
data("presidential_speech")
Xdat <- presidential_speech[1:10,1:4]




context("CBASS fitting")
test_that("CBASS Runs",{
  expect_error(
    cbass.fit <- CBASS(X=Xdat),
    NA
  )
})
test_that("CBASS Runs w/ cbassviz",{
  expect_error(
    cbass.fit <- CBASS(X=Xdat,control = list(alg.type='cbassviz')),
    NA
  )
})
test_that("CBASS Runs w/ cbassvizl1",{
  expect_error(
    cbass.fit <- CBASS(X=Xdat,control = list(alg.type='cbassvizl1')),
    NA
  )
})
test_that("CBASS Runs w/ cbass",{
  expect_error(
    cbass.fit <- CBASS(X=Xdat,control = list(alg.type='cbass')),
    NA
  )
})
test_that("CBASS Runs w/ cbassl1",{
  expect_error(
    cbass.fit <- CBASS(X=Xdat,control = list(alg.type='cbassl1')),
    NA
  )
})






context('CBASS Printing')
cbass.fit <- CBASS(X=Xdat)
test_that("print.CBASS Runs",{
  expect_error(
    print(cbass.fit),
    NA
  )
})




context('CBASS Plotting')
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






context('CBASS Clustering')
test_that("Clustering CBASS  Runs k.obs=1",{
  expect_error(
    Clustering(cbass.fit,k.obs=1),
    NA
  )
})
test_that("Clustering CBASS  Runs k.obs=2",{
  expect_error(
    Clustering(cbass.fit,k.obs=2),
    NA
  )
})
test_that("Clustering CBASS  Runs k.var=2",{
  expect_error(
    Clustering(cbass.fit,k.var=2),
    NA
  )
})
test_that("Clustering CBASS  Runs percent=.3",{
  expect_error(
    Clustering(cbass.fit,percent=.3),
    NA
  )
})




# context('CBASS Saving')
# test_that("saveviz runs with static obs.dend, k.obs",{
#   expect_error(
#     saveviz(cbass.fit,
#             file.name = 'cbass_obsdend_static_kobs.png',
#             plot.type = 'obs.dendrogram',
#             image.type = 'static',
#             k.obs=3
#             ),
#     NA
#   )
# })
# test_that("saveviz runs with static obs.dend, k.var",{
#   expect_error(
#     saveviz(cbass.fit,
#             file.name = 'cbass_obsdend_static_kvar.png',
#             plot.type = 'obs.dendrogram',
#             image.type = 'static',
#             k.var=3
#             ),
#     NA
#   )
# })
#
# test_that("saveviz runs with static obs.dend, percent",{
#   expect_error(
#     saveviz(cbass.fit,
#             file.name = 'cbass_obsdend_static_percent.png',
#             plot.type = 'obs.dendrogram',
#             image.type = 'static',
#             percent=.3
#             ),
#     NA
#   )
# })
#
#
# test_that("saveviz runs with static var.dend, k.obs",{
#   expect_error(
#     saveviz(cbass.fit,
#             file.name = 'cbass_vardend_static_kobs.png',
#             plot.type = 'var.dendrogram',
#             image.type = 'static',
#             k.obs=3
#             ),
#     NA
#   )
# })
# test_that("saveviz runs with static var.dend, k.var",{
#   expect_error(
#     saveviz(cbass.fit,
#             file.name = 'cbass_vardend_static_kvar.png',
#             plot.type = 'var.dendrogram',
#             image.type = 'static',
#             k.var=3
#             ),
#     NA
#   )
# })
#
# test_that("saveviz runs with static var.dend, percent",{
#   expect_error(
#     saveviz(cbass.fit,
#             file.name = 'cbass_vardend_static_percent.png',
#             plot.type = 'var.dendrogram',
#             image.type = 'static',
#             percent=.3
#             ),
#     NA
#   )
# })
#
# test_that("saveviz runs with static heatmap, k.obs",{
#   expect_error(
#     saveviz(cbass.fit,
#             file.name = 'cbass_heatmap_static_kobs.png',
#             plot.type = 'heatmap',
#             image.type = 'static',
#             k.obs = 4
#             ),
#     NA
#   )
# })
# test_that("saveviz runs with static heatmap, k.var",{
#   expect_error(
#     saveviz(cbass.fit,
#             file.name = 'cbass_heatmap_static_kvar.png',
#             plot.type = 'heatmap',
#             image.type = 'static',
#             k.var = 10
#             ),
#     NA
#   )
# })
# test_that("saveviz runs with static heatmap, percent",{
#   expect_error(
#     saveviz(cbass.fit,
#             file.name = 'cbass_heatmap_static_percent.png',
#             plot.type = 'heatmap',
#             image.type = 'static',
#             percent=.2
#             ),
#     NA
#   )
# })
#
# test_that("saveviz runs with dynamic obsdend",{
#   expect_error(
#     saveviz(cbass.fit,
#             file.name = 'cbass_obsdend_dynamic.png',
#             plot.type = 'obs.dendrogram',
#             image.type = 'dynamic'
#             ),
#     NA
#   )
# })
# test_that("saveviz runs with dynamic vardend",{
#   expect_error(
#     saveviz(cbass.fit,
#             file.name = 'cbass_vardend_dynamic.png',
#             plot.type = 'var.dendrogram',
#             image.type = 'dynamic'
#             ),
#     NA
#   )
# })
