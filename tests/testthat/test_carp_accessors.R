context("Test accessor functions for CARP objects")

test_that("Input checking works", {
  carp_fit <- CARP(presidential_speech)

  ## Get cluster labels
  ## Bad k
  expect_error(get_cluster_labels(carp_fit, k = 2.5))
  expect_error(get_cluster_labels(carp_fit, k = 0))
  expect_error(get_cluster_labels(carp_fit, k = -3))
  expect_error(get_cluster_labels(carp_fit, k = NROW(presidential_speech) + 4))
  expect_error(get_cluster_labels(carp_fit, k = NA))
  expect_error(get_cluster_labels(carp_fit, k = c(2, 5)))

  ## Bad percent
  expect_error(get_cluster_labels(carp_fit, percent = 1.5))
  expect_error(get_cluster_labels(carp_fit, percent = -0.5))
  expect_error(get_cluster_labels(carp_fit, percent = NA))
  expect_error(get_cluster_labels(carp_fit, percent = c(0.25, 0.75)))

  ## Unknown named argument
  expect_error(get_cluster_labels(carp_fit, 0.5))

  ## Unknown unnamed argument
  expect_error(get_cluster_labels(carp_fit, arg = 5))

  ## Get cluster centroids
  ## Bad k
  expect_error(get_cluster_centroids(carp_fit, k = 2.5))
  expect_error(get_cluster_centroids(carp_fit, k = 0))
  expect_error(get_cluster_centroids(carp_fit, k = -3))
  expect_error(get_cluster_centroids(carp_fit, k = NROW(presidential_speech) + 4))
  expect_error(get_cluster_centroids(carp_fit, k = NA))
  expect_error(get_cluster_centroids(carp_fit, k = c(2, 5)))

  ## Bad percent
  expect_error(get_cluster_centroids(carp_fit, percent = 1.5))
  expect_error(get_cluster_centroids(carp_fit, percent = -0.5))
  expect_error(get_cluster_centroids(carp_fit, percent = NA))
  expect_error(get_cluster_centroids(carp_fit, percent = c(0.25, 0.75)))

  ## Unknown named argument
  expect_error(get_cluster_centroids(carp_fit, 0.5))

  ## Unknown unnamed argument
  expect_error(get_cluster_centroids(carp_fit, arg = 5))

  ## Get clustered data
  ## Bad k
  expect_error(get_clustered_data(carp_fit, k = 2.5))
  expect_error(get_clustered_data(carp_fit, k = 0))
  expect_error(get_clustered_data(carp_fit, k = -3))
  expect_error(get_clustered_data(carp_fit, k = NROW(presidential_speech) + 4))
  expect_error(get_clustered_data(carp_fit, k = NA))
  expect_error(get_clustered_data(carp_fit, k = c(2, 5)))

  ## Bad percent
  expect_error(get_clustered_data(carp_fit, percent = 1.5))
  expect_error(get_clustered_data(carp_fit, percent = -0.5))
  expect_error(get_clustered_data(carp_fit, percent = NA))
  expect_error(get_clustered_data(carp_fit, percent = c(0.25, 0.75)))

  ## Unknown named argument
  expect_error(get_clustered_data(carp_fit, 0.5))

  ## Unknown unnamed argument
  expect_error(get_clustered_data(carp_fit, arg = 5))
})

test_that("get_cluster_labels.CARP works", {
  carp_fit <- CARP(presidential_speech)

  labels <- get_cluster_labels(carp_fit, k = 1)
  names(labels) <- NULL

  expect_equal(labels,
               factor(rep("cluster_1", NROW(presidential_speech))))

  ## Known k
  expect_equal(levels(get_cluster_labels(carp_fit, k = 3)),
               c("cluster_1", "cluster_2", "cluster_3"))

  ## Should have correct names
  labels <- get_cluster_labels(carp_fit, k = 3)
  expect_equal(rownames(presidential_speech), names(labels))

  ## Distinct clusters at beginning of path
  expect_equal(NROW(presidential_speech),
               num_unique(get_cluster_labels(carp_fit, k = NROW(presidential_speech))))

  expect_equal(NROW(presidential_speech),
               num_unique(get_cluster_labels(carp_fit, percent = 0)))

  ## Mono-cluster at end of path
  expect_equal(1, num_unique(get_cluster_labels(carp_fit, k = 1)))
  expect_equal(1, num_unique(get_cluster_labels(carp_fit, percent = 1)))
})

test_that("get_cluster_centroids.CARP works", {

  carp_fits <- list(
    CARP(presidential_speech, X.center = FALSE, X.scale = FALSE),
    CARP(presidential_speech, X.center = TRUE,  X.scale = FALSE),
    CARP(presidential_speech, X.center = FALSE, X.scale = TRUE),
    CARP(presidential_speech, X.center = TRUE,  X.scale = TRUE)
  )

  ## Cluster centroids are the centers from raw data, regardless of pre-processing flags
  for(carp_fit in carp_fits){
    labels <- as.integer(get_cluster_labels(carp_fit, k = 2))
    centroids <- get_cluster_centroids(carp_fit, k = 2)

    expect_equal(colnames(centroids), colnames(presidential_speech))

    for(k in c(1, 2)){
      expect_equal(centroids[k, ], colMeans(presidential_speech[labels == k, , drop=FALSE]))
    }
  }

})

test_that("get_clustered_data.CARP works", {

  carp_fits <- list(
    CARP(presidential_speech, X.center = FALSE, X.scale = FALSE),
    CARP(presidential_speech, X.center = TRUE,  X.scale = FALSE),
    CARP(presidential_speech, X.center = FALSE, X.scale = TRUE),
    CARP(presidential_speech, X.center = TRUE,  X.scale = TRUE)
  )


  for(carp_fit in carp_fits){
    ## Clustered data matrix should have rank k for k clusters
    for(k in 1:10){
      expect_equal(Matrix:::rankMatrix(get_clustered_data(carp_fit, k = k)),
                   k,
                   check.attributes = FALSE)
    }

    labels <- as.integer(get_cluster_labels(carp_fit, k = 2))
    centroids <- get_cluster_centroids(carp_fit, k = 2)
    clustered_data_matrix <- get_clustered_data(carp_fit, k = 2)
    expect_equal(colnames(clustered_data_matrix), colnames(presidential_speech))
    expect_equal(rownames(clustered_data_matrix), rownames(presidential_speech))

    ## Clustered data are the centers from raw data, regardless of pre-processing flags
    for(n in NROW(clustered_data_matrix)){
        expect_equal(clustered_data_matrix[n, ],
                     centroids[labels[n], ])
    }
  }
})
