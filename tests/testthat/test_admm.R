context("Full ADMM Converges")

test_that("Full ADMM converges for CARP", {
  skip_on_cran()
  skip_if_not_installed("cvxclustr")
  data(mammals, package = "cvxclustr")

  ## Example modified from help pages of cvxclustr
  X <- as.matrix(mammals[,-1])
  carp_fit <- CARP(X, exact = TRUE, X.center = FALSE, X.scale = FALSE)

  ## Calculate matching `cvxclustr` solution
  Xt <- t(X)

  ## Match CARP() selected weights
  w <- clustRviz:::weight_mat_to_vec(carp_fit$weights)
  gamma <- unique(carp_fit$cluster_membership$Gamma)

  ## Perform clustering
  suppressWarnings(cvxclust_fit <- cvxclustr::cvxclust(Xt, w, gamma, tol = 1e-7))

  ## cvxclustr seems to use a pretty loose stopping tolerance, so this is a loose check...
  for(i in seq_along(gamma)){
    expect_equal(carp_fit$U[,,i], t(cvxclust_fit$U[[i]]),
                 check.attributes = FALSE, tolerance = 1e-4)
  }
})

test_that("Full Back-Tracking ADMM converges for CARP", {
  skip_on_cran()
  skip_if_not_installed("cvxclustr")
  data(mammals, package = "cvxclustr")

  ## Example modified from help pages of cvxclustr
  X <- as.matrix(mammals[,-1])
  carp_fit <- CARP(X, exact = TRUE, back_track = TRUE, X.center = FALSE, X.scale = FALSE)

  ## Calculate matching `cvxclustr` solution
  Xt <- t(X)

  ## Match CARP() selected weights
  w <- clustRviz:::weight_mat_to_vec(carp_fit$weights)
  gamma <- unique(carp_fit$cluster_membership$Gamma)

  ## Perform clustering
  suppressWarnings(cvxclust_fit <- cvxclustr::cvxclust(Xt, w, gamma, tol = 1e-7))

  ## cvxclustr seems to use a pretty loose stopping tolerance, so this is a loose check...
  for(i in seq_along(gamma)){
    expect_equal(carp_fit$U[,,i], t(cvxclust_fit$U[[i]]),
                 check.attributes = FALSE, tolerance = 1e-4)
  }
})

test_that("Full ADMM converges for CBASS", {
  skip_on_cran()
  skip_if_not_installed("cvxclustr")
  skip_if_not_installed("cvxbiclustr")
  skip_if_not_installed("Matrix")
  data(mammals, package = "cvxclustr")
  library(Matrix)

  clustRviz_options(max_iter = 5e6, keep_debug_info = TRUE)
  on.exit(clustRviz_reset_options())

  ## Example modified from help pages of cvxclustr
  X <- as.matrix(mammals[,-1])
  cbass_fit <- CBASS(X, exact = TRUE, X.center.global = FALSE, t = 1.05)

  cbass_gamma <- cbass_fit$debug$path$gamma_path
  cbass_U     <- array(cbass_fit$debug$path$u_path, c(NROW(X), NCOL(X), length(cbass_gamma)))

  ## Calculate matching `cvxbiclustr` solution

  ## Match CBASS() selected edges and weights
  D_row <- cbass_fit$row_fusions$D
  D_col <- t(cbass_fit$col_fusions$D)

  w_row <- cbass_fit$row_fusions$weights; w_row <- clustRviz:::weight_mat_to_vec(w_row); w_row <- w_row[w_row != 0]
  w_col <- cbass_fit$col_fusions$weights; w_col <- clustRviz:::weight_mat_to_vec(w_col); w_col <- w_col[w_col != 0]

  w_row <- w_row / sum(w_row); w_row <- w_row / sqrt(NROW(X))
  w_col <- w_col / sum(w_col); w_col <- w_col / sqrt(NCOL(X))

  ## Perform clustering
  capture.output(cobra_fit <- cvxbiclustr::cobra(X,
                                                 E_row = Matrix(D_row, sparse = TRUE),
                                                 E_col = Matrix(D_col, sparse = TRUE),
                                                 w_row = w_row,
                                                 w_col = w_col,
                                                 gamma = cbass_gamma))

  obj <- function(U, lambda, l1 = FALSE){
    DU <- D_row %*% U; DTU <- D_col %*% t(U)
    if(l1){
      sum((X - U)^2)/2 + lambda * sum(w_row * abs(DU)) + lambda * sum(w_col * abs(t(DTU)))
    } else {
      sum((X - U)^2)/2 + lambda * sum(w_row * sqrt(rowSums(DU^2))) + lambda * sum(w_col * sqrt(rowSums(DTU^2)))
    }
  }

  ## Check that the objectives are of (roughly) equal quality...
  for(ix in seq_along(cbass_gamma)){
    lambda <- cbass_gamma[ix]

    my_U <- cbass_U[,,ix]
    ec_U <- cobra_fit$U[[ix]]

    expect_true(obj(my_U, lambda) <= 1.01 * obj(ec_U, lambda))
  }
})
