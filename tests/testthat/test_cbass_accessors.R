context("Test accessor functions for CBASS objects")

test_that("Input checking works", {
  cbass_fit <- CBASS(presidential_speech)

  ## Exactly one argument required
  expect_error(clustering(cbass_fit))
  expect_error(clustering(cbass_fit, k.row = 5, k.col = 5))
  expect_error(clustering(cbass_fit, k.row = 5, percent = 0.5))
  expect_error(clustering(cbass_fit, k.col = 5, percent = 0.5))
  expect_error(clustering(cbass_fit, k.row = 5, k.col = 5, percent = 0.5))

  ## Get cluster labels

  ## Exactly one argument required
  expect_error(get_cluster_labels(cbass_fit))
  expect_error(get_cluster_labels(cbass_fit, k.row = 5, k.col = 5))
  expect_error(get_cluster_labels(cbass_fit, k.row = 5, percent = 0.5))
  expect_error(get_cluster_labels(cbass_fit, k.col = 5, percent = 0.5))
  expect_error(get_cluster_labels(cbass_fit, k.row = 5, k.col = 5, percent = 0.5))

  ## Bad k.row
  expect_error(get_cluster_labels(cbass_fit, k.row = 2.5))
  expect_error(get_cluster_labels(cbass_fit, k.row = 0))
  expect_error(get_cluster_labels(cbass_fit, k.row = -3))
  expect_error(get_cluster_labels(cbass_fit, k.row = NROW(presidential_speech) + 4))
  expect_error(get_cluster_labels(cbass_fit, k.row = NA))
  expect_error(get_cluster_labels(cbass_fit, k.row = c(2, 5)))

  ## Bad k.col
  expect_error(get_cluster_labels(cbass_fit, k.col = 2.5))
  expect_error(get_cluster_labels(cbass_fit, k.col = 0))
  expect_error(get_cluster_labels(cbass_fit, k.col = -3))
  expect_error(get_cluster_labels(cbass_fit, k.col = NCOL(presidential_speech) + 4))
  expect_error(get_cluster_labels(cbass_fit, k.col = NA))
  expect_error(get_cluster_labels(cbass_fit, k.col = c(2, 5)))

  ## Bad percent
  expect_error(get_cluster_labels(cbass_fit, percent = 1.5))
  expect_error(get_cluster_labels(cbass_fit, percent = -0.5))
  expect_error(get_cluster_labels(cbass_fit, percent = NA))
  expect_error(get_cluster_labels(cbass_fit, percent = c(0.25, 0.75)))

  ## Unknown named argument
  expect_error(get_cluster_labels(cbass_fit, 0.5))

  ## Unknown unnamed argument
  expect_error(get_cluster_labels(cbass_fit, arg = 5))

  ## Get clustered data

  ## Exactly one argument required
  expect_error(get_clustered_data(cbass_fit))
  expect_error(get_clustered_data(cbass_fit, k.row = 5, k.col = 5))
  expect_error(get_clustered_data(cbass_fit, k.row = 5, percent = 0.5))
  expect_error(get_clustered_data(cbass_fit, k.col = 5, percent = 0.5))
  expect_error(get_clustered_data(cbass_fit, k.row = 5, k.col = 5, percent = 0.5))

  ## Bad k.row
  expect_error(get_clustered_data(cbass_fit, k.row = 2.5))
  expect_error(get_clustered_data(cbass_fit, k.row = 0))
  expect_error(get_clustered_data(cbass_fit, k.row = -3))
  expect_error(get_clustered_data(cbass_fit, k.row = NROW(presidential_speech) + 4))
  expect_error(get_clustered_data(cbass_fit, k.row = NA))
  expect_error(get_clustered_data(cbass_fit, k.row = c(2, 5)))

  ## Bad k.col
  expect_error(get_clustered_data(cbass_fit, k.col = 2.5))
  expect_error(get_clustered_data(cbass_fit, k.col = 0))
  expect_error(get_clustered_data(cbass_fit, k.col = -3))
  expect_error(get_clustered_data(cbass_fit, k.col = NCOL(presidential_speech) + 4))
  expect_error(get_clustered_data(cbass_fit, k.col = NA))
  expect_error(get_clustered_data(cbass_fit, k.col = c(2, 5)))

  ## Bad percent
  expect_error(get_clustered_data(cbass_fit, percent = 1.5))
  expect_error(get_clustered_data(cbass_fit, percent = -0.5))
  expect_error(get_clustered_data(cbass_fit, percent = NA))
  expect_error(get_clustered_data(cbass_fit, percent = c(0.25, 0.75)))

  ## Unknown named argument
  expect_error(get_clustered_data(cbass_fit, 0.5))

  ## Unknown unnamed argument
  expect_error(get_clustered_data(cbass_fit, arg = 5))
})

test_that("get_cluster_labels.CBASS works on row labels", {
  cbass_fit <- CBASS(presidential_speech)

  labels <- get_cluster_labels(cbass_fit, k.row = 1, type = "row")
  names(labels) <- NULL

  expect_equal(labels,
               factor(rep("cluster_1", NROW(presidential_speech))))

  ## Known k
  expect_equal(levels(get_cluster_labels(cbass_fit, k.row = 3, type = "row")),
               c("cluster_1", "cluster_2", "cluster_3"))

  for (k in 1:10) {
    expect_equal(k, nlevels(get_cluster_labels(cbass_fit, k.row = k, type = "row")))
  }

  ## Should have correct names
  labels <- get_cluster_labels(cbass_fit, k.row = 3, type = "row")
  expect_equal(rownames(presidential_speech), names(labels))

  ## Distinct clusters at beginning of path
  expect_equal(NROW(presidential_speech),
               num_unique(get_cluster_labels(cbass_fit,
                                             k.row = NROW(presidential_speech),
                                             type = "row")))

  expect_equal(NROW(presidential_speech),
               num_unique(get_cluster_labels(cbass_fit, percent = 0, type = "row")))

  ## Mono-cluster at end of path
  expect_equal(1, num_unique(get_cluster_labels(cbass_fit, k.row = 1, type = "row")))
  expect_equal(1, num_unique(get_cluster_labels(cbass_fit, percent = 1, type = "row")))

  ## Correct number of labels returned
  for(pct in seq(0, 1, length.out = 21)){
    expect_equal(length(get_cluster_labels(cbass_fit, percent = pct, type = "row")), NROW(presidential_speech))
  }
})

test_that("get_cluster_labels.CBASS works on column labels", {
  cbass_fit <- CBASS(presidential_speech)

  labels <- get_cluster_labels(cbass_fit, k.col = 1, type = "col")
  names(labels) <- NULL

  expect_equal(labels,
               factor(rep("cluster_1", NCOL(presidential_speech))))

  ## Known k
  expect_equal(levels(get_cluster_labels(cbass_fit, k.col = 3, type = "col")),
               c("cluster_1", "cluster_2", "cluster_3"))

  for (k in 1:10) {
    expect_equal(k, nlevels(get_cluster_labels(cbass_fit, k.col = k, type = "col")))
  }

  ## Should have correct names
  labels <- get_cluster_labels(cbass_fit, k.col = 3, type = "col")
  expect_equal(colnames(presidential_speech), names(labels))

  ## Distinct clusters at beginning of path
  expect_equal(NCOL(presidential_speech),
               num_unique(get_cluster_labels(cbass_fit,
                                             k.col = NCOL(presidential_speech),
                                             type = "col")))

  expect_equal(NCOL(presidential_speech),
               num_unique(get_cluster_labels(cbass_fit, percent = 0, type = "col")))

  ## Mono-cluster at end of path
  expect_equal(1, num_unique(get_cluster_labels(cbass_fit, k.col = 1, type = "col")))
  expect_equal(1, num_unique(get_cluster_labels(cbass_fit, percent = 1, type = "col")))

  ## Correct number of labels returned
  for(pct in seq(0, 1, length.out = 21)){
    expect_equal(length(get_cluster_labels(cbass_fit, percent = pct, type = "col")), NCOL(presidential_speech))
  }
})

test_that("get_cluster_centroids.CBASS works", {
  cbass_fits <- list(
    CBASS(presidential_speech, X.center.global = FALSE),
    CBASS(presidential_speech, X.center.global = TRUE)
  )

  for(cbass_fit in cbass_fits){
    ## Cluster centroids matrix has k1-times-k2 distinct elements
    ## for k1 row clusters and k2 col clusters
    ##
    ## These elements are typically unique, but weird things can happen
    ## (particularly near the end of the path) so we don't test
    for(k in 1:10){
      row_labels <- num_unique(get_cluster_labels(cbass_fit, k.row = k, type = "row"))
      col_labels <- num_unique(get_cluster_labels(cbass_fit, k.row = k, type = "col"))
      expect_equal(length(get_cluster_centroids(cbass_fit, k.row = k)), row_labels * col_labels)

      row_labels <- num_unique(get_cluster_labels(cbass_fit, k.col = k, type = "row"))
      col_labels <- num_unique(get_cluster_labels(cbass_fit, k.col = k, type = "col"))
      expect_equal(length(get_cluster_centroids(cbass_fit, k.col = k)), row_labels * col_labels)

      expect_is(get_cluster_centroids(cbass_fit, k.row = k), "matrix")
      expect_is(get_cluster_centroids(cbass_fit, k.col = k), "matrix")
    }

    ## At full fusion, we should get a single element that is the grand mean
    expect_equal(mean(presidential_speech),
                 as.vector(get_cluster_centroids(cbass_fit, percent = 1)))

    ## At no fusion, we get the original data (we aren't enforcing ordering)
    expect_equal(sort(as.vector(presidential_speech)),
                 sort(as.vector(get_cluster_centroids(cbass_fit, percent = 0))))

    ## No row or column names
    expect_null(rownames(get_cluster_centroids(cbass_fit, percent = 0.5)))
    expect_null(colnames(get_cluster_centroids(cbass_fit, percent = 0.5)))
  }
})

test_that("get_clustered_data.CBASS works", {
  cbass_fits <- list(
    CBASS(presidential_speech, X.center.global = FALSE),
    CBASS(presidential_speech, X.center.global = TRUE)
  )

  for(cbass_fit in cbass_fits){
    ## Clustered data matrix has at most k1-times-k2 distinct elements
    ## for k1 row clusters and k2 col clusters
    ##
    ## It can have fewer if a pair of clusters on the "same side" (row or col) are
    ## nested in a single cluster on the other side
    for(k in 1:10){
      row_labels <- num_unique(get_cluster_labels(cbass_fit, k.row = k, type = "row"))
      col_labels <- num_unique(get_cluster_labels(cbass_fit, k.row = k, type = "col"))
      expect_lte(num_unique(get_clustered_data(cbass_fit, k.row = k)), row_labels * col_labels)

      row_labels <- num_unique(get_cluster_labels(cbass_fit, k.col = k, type = "row"))
      col_labels <- num_unique(get_cluster_labels(cbass_fit, k.col = k, type = "col"))
      expect_lte(num_unique(get_clustered_data(cbass_fit, k.col = k)), row_labels * col_labels)

      # Correct size
      expect_equal(dim(presidential_speech), dim(get_clustered_data(cbass_fit, k.row = k)))
      expect_equal(dim(presidential_speech), dim(get_clustered_data(cbass_fit, k.col = k)))

      # Correct names
      expect_equal(dimnames(presidential_speech), dimnames(get_clustered_data(cbass_fit, k.row = k)))
      expect_equal(dimnames(presidential_speech), dimnames(get_clustered_data(cbass_fit, k.col = k)))
    }

    # Expect mono-cluster at end of path
    presidential_speech_full_cluster <- presidential_speech * 0 + mean(presidential_speech)

    ## The variables are fully clustered before the observations, so this doesn't actully work
    # expect_equal(presidential_speech_full_cluster, get_clustered_data(cbass_fit, k.col = 1))
    expect_equal(presidential_speech_full_cluster, get_clustered_data(cbass_fit, k.row = 1))
    expect_equal(presidential_speech_full_cluster, get_clustered_data(cbass_fit, percent = 1))

    ## Clustered data are the centers from raw data, regardless of pre-processing flags
    ## Fix k.row
    row_labels <- as.integer(get_cluster_labels(cbass_fit, k.row = 2, type = "row"))
    col_labels <- as.integer(get_cluster_labels(cbass_fit, k.row = 2, type = "col"))
    clustered_data_matrix <- get_clustered_data(cbass_fit, k.row = 2)
    expect_equal(colnames(clustered_data_matrix), colnames(presidential_speech))
    expect_equal(rownames(clustered_data_matrix), rownames(presidential_speech))

    for(r in unique(row_labels)){
      for(c in unique(col_labels)){
        expect_equal(num_unique(clustered_data_matrix[row_labels == r, col_labels == c]), 1)
        expect_equal(mean(clustered_data_matrix[row_labels == r, col_labels == c]),
                     mean(presidential_speech[row_labels == r, col_labels == c]))
      }
    }

    ## Fix k.col
    row_labels <- as.integer(get_cluster_labels(cbass_fit, k.col = 2, type = "row"))
    col_labels <- as.integer(get_cluster_labels(cbass_fit, k.col = 2, type = "col"))
    clustered_data_matrix <- get_clustered_data(cbass_fit, k.col = 2)
    expect_equal(colnames(clustered_data_matrix), colnames(presidential_speech))
    expect_equal(rownames(clustered_data_matrix), rownames(presidential_speech))

    ## Clustered data are the centers from raw data, regardless of pre-processing flags
    for(r in unique(row_labels)){
      for(c in unique(col_labels)){
        expect_equal(num_unique(clustered_data_matrix[row_labels == r, col_labels == c]), 1)
        expect_equal(mean(clustered_data_matrix[row_labels == r, col_labels == c]),
                     mean(presidential_speech[row_labels == r, col_labels == c]))
      }
    }

    ## Fix percent
    row_labels <- as.integer(get_cluster_labels(cbass_fit, percent = 0.9, type = "row"))
    col_labels <- as.integer(get_cluster_labels(cbass_fit, percent = 0.9, type = "col"))
    clustered_data_matrix <- get_clustered_data(cbass_fit, percent = 0.9)
    expect_equal(colnames(clustered_data_matrix), colnames(presidential_speech))
    expect_equal(rownames(clustered_data_matrix), rownames(presidential_speech))

    ## Clustered data are the centers from raw data, regardless of pre-processing flags
    for(r in unique(row_labels)){
      for(c in unique(col_labels)){
        expect_equal(num_unique(clustered_data_matrix[row_labels == r, col_labels == c]), 1)
        expect_equal(mean(clustered_data_matrix[row_labels == r, col_labels == c]),
                     mean(presidential_speech[row_labels == r, col_labels == c]))
      }
    }
  }
})

test_that("CBASS cluster centroids are correctly calculated with refit = FALSE", {
  ## All of these results should work regardless of pre-processing options
  cbass_fits <- list(
    CBASS(presidential_speech, X.center.global = FALSE),
    CBASS(presidential_speech, X.center.global = TRUE)
  )

  presidential_speech_full_cluster <- presidential_speech * 0 + mean(presidential_speech)

  for(cbass_fit in cbass_fits){
    expect_equal(get_clustered_data(cbass_fit, percent = 0, refit = FALSE), presidential_speech)
    expect_equal(get_clustered_data(cbass_fit, percent = 1, refit = FALSE), presidential_speech_full_cluster)

    for(pct in seq(0.1, 1, by = 0.1)){
      ## Variance should decrease as we increase regularization
      expect_gte(sd(get_clustered_data(cbass_fit, percent = pct - 0.1, refit = FALSE)),
                 sd(get_clustered_data(cbass_fit, percent = pct, refit = FALSE)))
    }
  }

  ## Results should be the same regardless of centering
  for(pct in seq(0, 1, length.out = 11)){
    expect_equal(get_clustered_data(cbass_fits[[1]], refit = FALSE, percent = pct),
                 get_clustered_data(cbass_fits[[2]], refit = FALSE, percent = pct))
  }

  ## grand means should always be the same
  for(pct in seq(0, 1, length.out = 11)){
    for(cbass_fit in cbass_fits){
      expect_equal(mean(presidential_speech),
                   mean(get_clustered_data(cbass_fit, percent = pct, refit = FALSE)))
    }
  }

})
