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

  ## Error with multiple arguments
  expect_error(get_clustered_data(carp_fit, percent = 0.5, k = 3))
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

  ## Correct number of clusters
  for(pct in seq(0, 1, length.out = 21)){
    expect_equal(length(get_cluster_labels(carp_fit, percent = pct)), NROW(presidential_speech))
  }
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
    raw_centroids <- get_cluster_centroids(carp_fit, k = 2, refit = FALSE)

    expect_equal(colnames(centroids), colnames(presidential_speech))
    expect_equal(colnames(raw_centroids), colnames(presidential_speech))

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

test_that("CARP cluster centroids are correctly calculated with refit = FALSE", {
  ## All of these results should work regardless of pre-processing options
  carp_fits <- list(
    CARP(presidential_speech, X.center = FALSE, X.scale = FALSE),
    CARP(presidential_speech, X.center = TRUE,  X.scale = FALSE),
    CARP(presidential_speech, X.center = FALSE, X.scale = TRUE),
    CARP(presidential_speech, X.center = TRUE,  X.scale = TRUE)
  )

  presidential_speech_full_cluster <- presidential_speech
  for(i in seq_len(NROW(presidential_speech))){
    presidential_speech_full_cluster[i, ] <- colMeans(presidential_speech)
  }

  for(carp_fit in carp_fits){
    expect_equal(get_clustered_data(carp_fit, percent = 0, refit = FALSE), presidential_speech)
    expect_equal(get_clustered_data(carp_fit, percent = 1, refit = FALSE), presidential_speech_full_cluster)

    for(pct in seq(0.1, 1, by = 0.1)){
      ## Variance should decrease as we increase regularization
      expect_gte(sd(get_clustered_data(carp_fit, percent = pct - 0.1, refit = FALSE)),
                 sd(get_clustered_data(carp_fit, percent = pct, refit = FALSE)))
    }
  }

  carp_fits_no_scale <- carp_fits[c(1, 2)]
  carp_fits_scale    <- carp_fits[c(3, 4)]
  ## Results should be the same regardless of centering
  ## Scaling features can make a difference in cluster assignment
  for(pct in seq(0, 1, length.out = 11)){
    expect_true(list_all_equal(lapply(carp_fits_no_scale, get_clustered_data, percent = pct, refit = FALSE)))
  }

  ## Column means should always be the same
  for(pct in seq(0, 1, length.out = 11)){
    for(carp_fit in carp_fits){
      expect_equal(colMeans(presidential_speech),
                   colMeans(get_clustered_data(carp_fit, percent = pct, refit = FALSE)))
    }
  }

  ## Column-wise standard deviations are decreasing
  colSds <- function(x) apply(x, 2, sd)
  for(pct in seq(0.1, 1, length.out = 10)){
    for(carp_fit in carp_fits){
      expect_gt(sum(colSds(get_clustered_data(carp_fit, percent = pct - 0.1, refit = FALSE))),
                sum(colSds(get_clustered_data(carp_fit, percent = pct, refit = FALSE))))
    }
  }
})
