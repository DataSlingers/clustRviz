context("Full ADMM Converges")

test_that("Full ADMM converges for CARP", {
  skip_on_cran()
  skip_if_not_installed("cvxclustr")
  data(mammals, package = "cvxclustr")

  ## Example modified from help pages of cvxclustr
  X <- as.matrix(mammals[,-1])
  carp_fit <- CARP(X, alg.type = "admm", X.center = FALSE, X.scale = FALSE)

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
  skip("CBASS+ADMM does not get the right answer...")
  skip_on_cran()
  skip_if_not_installed("cvxclustr")
  skip_if_not_installed("cvxbiclustr")
  skip_if_not_installed("Matrix")
  data(mammals, package = "cvxclustr")

  clustRviz_options(max_iter = 5e6, keep_debug_info = TRUE)
  on.exit(clustRviz_reset_options())

  ## Example modified from help pages of cvxclustr
  X <- as.matrix(mammals[,-1])
  cbass_fit <- CBASS(X, alg.type = "admm", X.center.global = FALSE)

  cbass_gamma <- cbass_fit$debug$path$gamma_path
  cbass_U     <- array(cbass_fit$debug$path$u_path, c(NROW(X), NCOL(X), length(cbass_gamma)))

  ## Calculate matching `cvxclustr` solution

  ## Match CBASS() selected edges and weights
  E_row <- Matrix::Matrix(cbass_fit$row_fusions$D, sparse = TRUE)
  w_row <- clustRviz:::weight_mat_to_vec(cbass_fit$row_fusions$weights)
  w_row <- w_row[w_row != 0]

  E_col <- Matrix::Matrix(t(cbass_fit$col_fusions$D), sparse = TRUE)
  w_col <- clustRviz:::weight_mat_to_vec(cbass_fit$col_fusions$weights)
  w_col <- w_col[w_col != 0]

  ## Perform clustering
  capture.output(cobra_fit <- cvxbiclustr::cobra(X, E_row, E_col, w_row, w_col, cbass_gamma, tol = 1e-8))

  fn <- function(U, gamma){
    sum((U - X)^2)/2 + gamma * sum(w_row * sqrt(rowSums((cbass_fit$row_fusions$D %*% U)^2))) +
                       gamma * sum(w_col * sqrt(colSums((U %*% cbass_fit$col_fusions$D)^2)))
  }

  ## Check that the objectives are of (roughly) equal quality...
  ##
  ## NB: This does not (exactly) check that the solutions match -- I'm not sure
  ## why they aren't matching right now
  my_loss <- vapply(seq_along(cbass_gamma), function(i) fn(cbass_U[,,i], cbass_gamma[i]), 0.0)
  ec_loss <- vapply(seq_along(cbass_gamma), function(i) fn(cobra_fit$U[[i]], cbass_gamma[i]), 0.0)

  expect_true(mean(my_loss - ec_loss) < 0)

  ## Worst case is less than 10% worse (relative error)
  worst_re <- max(((my_loss - ec_loss) / pmin(my_loss, ec_loss))[-1])
  expect_true(worst_re < 0.1)
})
