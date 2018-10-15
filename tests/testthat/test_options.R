context("Options Handling")

test_that("Test option error handling", {
  ## ADMM relaxation parameter
  expect_error(clustRviz_options(rho = 0))
  expect_error(clustRviz_options(rho = -1))
  expect_error(clustRviz_options(rho = "a"))
  expect_error(clustRviz_options(rho = NA))
  expect_error(clustRviz_options(rho = c(1, 2)))

  ## Initial (burn-in phase) regularization parameter
  expect_error(clustRviz_options(epsilon = 0))
  expect_error(clustRviz_options(epsilon = -1))
  expect_error(clustRviz_options(epsilon = "a"))
  expect_error(clustRviz_options(epsilon = NA))
  expect_error(clustRviz_options(epsilon = c(1, 2)))

  ## Back-tracking parameters
  expect_error(clustRviz_options(viz_initial_step = 1))
  expect_error(clustRviz_options(viz_initial_step = 0.5))
  expect_error(clustRviz_options(viz_initial_step = 0))
  expect_error(clustRviz_options(viz_initial_step = -1))
  expect_error(clustRviz_options(viz_initial_step = -1.5))
  expect_error(clustRviz_options(viz_initial_step = "a"))
  expect_error(clustRviz_options(viz_initial_step = NA))
  expect_error(clustRviz_options(viz_initial_step = c(1.5, 2)))

  expect_error(clustRviz_options(viz_small_step = 1))
  expect_error(clustRviz_options(viz_small_step = 0.5))
  expect_error(clustRviz_options(viz_small_step = 0))
  expect_error(clustRviz_options(viz_small_step = -1))
  expect_error(clustRviz_options(viz_small_step = -1.5))
  expect_error(clustRviz_options(viz_small_step = "a"))
  expect_error(clustRviz_options(viz_small_step = NA))
  expect_error(clustRviz_options(viz_small_step = c(1.5, 2)))

  expect_error(clustRviz_options(viz_max_inner_iter = 0))
  expect_error(clustRviz_options(viz_max_inner_iter = -5))
  expect_error(clustRviz_options(viz_max_inner_iter = 35.5))
  expect_error(clustRviz_options(viz_max_inner_iter = "a"))
  expect_error(clustRviz_options(viz_max_inner_iter = NA))
  expect_error(clustRviz_options(viz_max_inner_iter = c(500, 600)))

  # Stopping and storage parameters
  expect_error(clustRviz_options(max_iter = 0))
  expect_error(clustRviz_options(max_iter = -5))
  expect_error(clustRviz_options(max_iter = 35.5))
  expect_error(clustRviz_options(max_iter = "a"))
  expect_error(clustRviz_options(max_iter = NA))
  expect_error(clustRviz_options(max_iter = c(500, 600)))

  expect_error(clustRviz_options(burn_in = 0))
  expect_error(clustRviz_options(burn_in = -5))
  expect_error(clustRviz_options(burn_in = 35.5))
  expect_error(clustRviz_options(burn_in = "a"))
  expect_error(clustRviz_options(burn_in = NA))
  expect_error(clustRviz_options(burn_in = c(500, 600)))

  expect_error(clustRviz_options(keep = 0))
  expect_error(clustRviz_options(keep = -5))
  expect_error(clustRviz_options(keep = 35.5))
  expect_error(clustRviz_options(keep = "a"))
  expect_error(clustRviz_options(keep = NA))
  expect_error(clustRviz_options(keep = c(500, 600)))

  expect_error(clustRviz_options(keep_debug_info = 0))
  expect_error(clustRviz_options(keep_debug_info = 3))
  expect_error(clustRviz_options(keep_debug_info = -5))
  expect_error(clustRviz_options(keep_debug_info = 35.5))
  expect_error(clustRviz_options(keep_debug_info = "a"))
  expect_error(clustRviz_options(keep_debug_info = NA))
  expect_error(clustRviz_options(keep_debug_info = c(500, 600)))
})

test_that("clustRviz_reset_options works", {
  base_opts <- clustRviz_options()

  clustRviz_options(rho = 2, keep = 5, max_iter = 500, burn_in = 20)
  expect_false(isTRUE(all.equal(base_opts, clustRviz_options())))

  clustRviz_reset_options()
  expect_equal(base_opts, clustRviz_options())
})

test_that("Setting clustRviz options works", {
  clustRviz_options(rho = 2)
  expect_equal(clustRviz:::.clustRvizOptionsEnv[["rho"]], 2)

  clustRviz_options(viz_initial_step = 2)
  expect_equal(clustRviz:::.clustRvizOptionsEnv[["viz_initial_step"]], 2)

  clustRviz_options(viz_max_inner_iter = 20)
  expect_equal(clustRviz:::.clustRvizOptionsEnv[["viz_max_inner_iter"]], 20)

  clustRviz_options(viz_small_step = 1.2)
  expect_equal(clustRviz:::.clustRvizOptionsEnv[["viz_small_step"]], 1.2)

  clustRviz_reset_options()
})

test_that("Options warnings work", {
  expect_warning(clustRviz_options(viz_initial_step = 1.2, viz_small_step = 1.2))
  expect_warning(clustRviz_options(viz_initial_step = 1.2, viz_small_step = 2))

  expect_warning(clustRviz_options(max_iter = 100, burn_in = 100))
  expect_warning(clustRviz_options(max_iter = 100, burn_in = 150))

  clustRviz_reset_options()
})
