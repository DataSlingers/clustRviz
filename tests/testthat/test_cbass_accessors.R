context("Test accessor functions for CBASS objects")

test_that("Input checking works", {
  cbass_fit <- CBASS(presidential_speech)

  ## Get cluster labels
  ## Bad k.obs
  expect_error(get_cluster_labels(cbass_fit, k.obs = 2.5))
  expect_error(get_cluster_labels(cbass_fit, k.obs = 0))
  expect_error(get_cluster_labels(cbass_fit, k.obs = -3))
  expect_error(get_cluster_labels(cbass_fit, k.obs = NROW(presidential_speech) + 4))
  expect_error(get_cluster_labels(cbass_fit, k.obs = NA))
  expect_error(get_cluster_labels(cbass_fit, k.obs = c(2, 5)))

  ## Bad k.var
  expect_error(get_cluster_labels(cbass_fit, k.var = 2.5))
  expect_error(get_cluster_labels(cbass_fit, k.var = 0))
  expect_error(get_cluster_labels(cbass_fit, k.var = -3))
  expect_error(get_cluster_labels(cbass_fit, k.var = NCOL(presidential_speech) + 4))
  expect_error(get_cluster_labels(cbass_fit, k.var = NA))
  expect_error(get_cluster_labels(cbass_fit, k.var = c(2, 5)))

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
  ## Bad k.obs
  expect_error(get_clustered_data(cbass_fit, k.obs = 2.5))
  expect_error(get_clustered_data(cbass_fit, k.obs = 0))
  expect_error(get_clustered_data(cbass_fit, k.obs = -3))
  expect_error(get_clustered_data(cbass_fit, k.obs = NROW(presidential_speech) + 4))
  expect_error(get_clustered_data(cbass_fit, k.obs = NA))
  expect_error(get_clustered_data(cbass_fit, k.obs = c(2, 5)))

  ## Bad k.var
  expect_error(get_clustered_data(cbass_fit, k.var = 2.5))
  expect_error(get_clustered_data(cbass_fit, k.var = 0))
  expect_error(get_clustered_data(cbass_fit, k.var = -3))
  expect_error(get_clustered_data(cbass_fit, k.var = NCOL(presidential_speech) + 4))
  expect_error(get_clustered_data(cbass_fit, k.var = NA))
  expect_error(get_clustered_data(cbass_fit, k.var = c(2, 5)))

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

test_that("get_cluster_labels.CBASS works on observation labels", {
  cbass_fit <- CBASS(presidential_speech)

  labels <- get_cluster_labels(cbass_fit, k.obs = 1, type = "obs")
  names(labels) <- NULL

  expect_equal(labels,
               factor(rep("cluster_1", NROW(presidential_speech))))

  ## Known k
  expect_equal(levels(get_cluster_labels(cbass_fit, k.obs = 3, type = "obs")),
               c("cluster_1", "cluster_2", "cluster_3"))

  ## Should have correct names
  labels <- get_cluster_labels(cbass_fit, k.obs = 3, type = "obs")
  expect_equal(rownames(presidential_speech), names(labels))

  ## Distinct clusters at beginning of path
  expect_equal(NROW(presidential_speech),
               num_unique(get_cluster_labels(cbass_fit,
                                             k.obs = NROW(presidential_speech),
                                             type = "obs")))

  expect_equal(NROW(presidential_speech),
               num_unique(get_cluster_labels(cbass_fit, percent = 0, type = "obs")))

  ## Mono-cluster at end of path
  expect_equal(1, num_unique(get_cluster_labels(cbass_fit, k.obs = 1, type = "obs")))
  expect_equal(1, num_unique(get_cluster_labels(cbass_fit, percent = 1, type = "obs")))
})

test_that("get_cluster_labels.CBASS works on observation labels", {
  cbass_fit <- CBASS(presidential_speech)

  labels <- get_cluster_labels(cbass_fit, k.obs = 1, type = "obs")
  names(labels) <- NULL

  expect_equal(labels,
               factor(rep("cluster_1", NROW(presidential_speech))))

  ## Known k
  expect_equal(levels(get_cluster_labels(cbass_fit, k.obs = 3, type = "obs")),
               c("cluster_1", "cluster_2", "cluster_3"))

  ## Should have correct names
  labels <- get_cluster_labels(cbass_fit, k.obs = 3, type = "obs")
  expect_equal(rownames(presidential_speech), names(labels))

  ## Distinct clusters at beginning of path
  expect_equal(NROW(presidential_speech),
               num_unique(get_cluster_labels(cbass_fit,
                                             k.obs = NROW(presidential_speech),
                                             type = "obs")))

  expect_equal(NROW(presidential_speech),
               num_unique(get_cluster_labels(cbass_fit, percent = 0, type = "obs")))

  ## Mono-cluster at end of path
  expect_equal(1, num_unique(get_cluster_labels(cbass_fit, k.obs = 1, type = "obs")))
  expect_equal(1, num_unique(get_cluster_labels(cbass_fit, percent = 1, type = "obs")))
})

test_that("get_cluster_labels.CBASS works on variable labels", {
  cbass_fit <- CBASS(presidential_speech)

  labels <- get_cluster_labels(cbass_fit, k.var = 1, type = "var")
  names(labels) <- NULL

  expect_equal(labels,
               factor(rep("cluster_1", NCOL(presidential_speech))))

  ## Known k
  expect_equal(levels(get_cluster_labels(cbass_fit, k.var = 3, type = "var")),
               c("cluster_1", "cluster_2", "cluster_3"))

  ## Should have correct names
  labels <- get_cluster_labels(cbass_fit, k.var = 3, type = "var")
  expect_equal(colnames(presidential_speech), names(labels))

  ## Distinct clusters at beginning of path
  expect_equal(NCOL(presidential_speech),
               num_unique(get_cluster_labels(cbass_fit,
                                             k.var = NCOL(presidential_speech),
                                             type = "var")))

  expect_equal(NCOL(presidential_speech),
               num_unique(get_cluster_labels(cbass_fit, percent = 0, type = "var")))

  ## Mono-cluster at end of path
  expect_equal(1, num_unique(get_cluster_labels(cbass_fit, k.var = 1, type = "var")))
  expect_equal(1, num_unique(get_cluster_labels(cbass_fit, percent = 1, type = "var")))
})

test_that("get_clustered_data.CBASS works", {
  cbass_fits <- list(
    CBASS(presidential_speech, X.center.global = FALSE),
    CBASS(presidential_speech, X.center.global = TRUE)
  )

  for(cbass_fit in cbass_fits){
    ## Clustered data matrix has at most k1-times-k2 distinct elements
    ## for k1 obs clusters and k2 var clusters
    ##
    ## It can have fewer if a pair of clusters on the "same side" (obs or var) are
    ## nested in a single cluster on the other side
    for(k in 1:10){
      obs_labels <- num_unique(get_cluster_labels(cbass_fit, k.obs = k, type = "obs"))
      var_labels <- num_unique(get_cluster_labels(cbass_fit, k.obs = k, type = "var"))
      expect_lte(num_unique(get_clustered_data(cbass_fit, k.obs = k)), obs_labels * var_labels)

      obs_labels <- num_unique(get_cluster_labels(cbass_fit, k.var = k, type = "obs"))
      var_labels <- num_unique(get_cluster_labels(cbass_fit, k.var = k, type = "var"))
      expect_lte(num_unique(get_clustered_data(cbass_fit, k.var = k)), obs_labels * var_labels)
    }

    ## Clustered data are the centers from raw data, regardless of pre-processing flags
    ## Fix k.obs
    obs_labels <- as.integer(get_cluster_labels(cbass_fit, k.obs = 2, type = "obs"))
    var_labels <- as.integer(get_cluster_labels(cbass_fit, k.obs = 2, type = "var"))
    clustered_data_matrix <- get_clustered_data(cbass_fit, k.obs = 2)
    expect_equal(colnames(clustered_data_matrix), colnames(presidential_speech))
    expect_equal(rownames(clustered_data_matrix), rownames(presidential_speech))

    for(o in unique(obs_labels)){
      for(v in unique(var_labels)){
        expect_equal(num_unique(clustered_data_matrix[obs_labels == o, var_labels == v]), 1)
        expect_equal(mean(clustered_data_matrix[obs_labels == o, var_labels == v]),
                     mean(presidential_speech[obs_labels == o, var_labels == v]))
      }
    }

    ## Fix k.var
    obs_labels <- as.integer(get_cluster_labels(cbass_fit, k.var = 2, type = "obs"))
    var_labels <- as.integer(get_cluster_labels(cbass_fit, k.var = 2, type = "var"))
    clustered_data_matrix <- get_clustered_data(cbass_fit, k.var = 2)
    expect_equal(colnames(clustered_data_matrix), colnames(presidential_speech))
    expect_equal(rownames(clustered_data_matrix), rownames(presidential_speech))

    ## Clustered data are the centers from raw data, regardless of pre-processing flags
    for(o in unique(obs_labels)){
      for(v in unique(var_labels)){
        expect_equal(num_unique(clustered_data_matrix[obs_labels == o, var_labels == v]), 1)
        expect_equal(mean(clustered_data_matrix[obs_labels == o, var_labels == v]),
                     mean(presidential_speech[obs_labels == o, var_labels == v]))
      }
    }

    ## Fix percent
    obs_labels <- as.integer(get_cluster_labels(cbass_fit, percent = 0.9, type = "obs"))
    var_labels <- as.integer(get_cluster_labels(cbass_fit, percent = 0.9, type = "var"))
    clustered_data_matrix <- get_clustered_data(cbass_fit, percent = 0.9)
    expect_equal(colnames(clustered_data_matrix), colnames(presidential_speech))
    expect_equal(rownames(clustered_data_matrix), rownames(presidential_speech))

    ## Clustered data are the centers from raw data, regardless of pre-processing flags
    for(o in unique(obs_labels)){
      for(v in unique(var_labels)){
        expect_equal(num_unique(clustered_data_matrix[obs_labels == o, var_labels == v]), 1)
        expect_equal(mean(clustered_data_matrix[obs_labels == o, var_labels == v]),
                     mean(presidential_speech[obs_labels == o, var_labels == v]))
      }
    }
  }
})
