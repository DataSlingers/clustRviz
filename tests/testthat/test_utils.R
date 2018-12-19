context("Test Utility Functions")

capitalize_string <- function(x){
  x <- gsub("_", " ", x)
  vapply(strsplit(x, " "),
         function(x) paste(paste0(toupper(substring(x, 1, 1)), substring(x, 2)), collapse = " "),
         character(1))
}

test_that("Validators work", {
  is_logical_scalar  <- clustRviz:::is_logical_scalar
  is_numeric_scalar  <- clustRviz:::is_numeric_scalar
  is_integer_scalar  <- clustRviz:::is_integer_scalar
  is_percent_scalar  <- clustRviz:::is_percent_scalar
  is_positive_scalar <- clustRviz:::is_positive_scalar
  is_positive_integer_scalar <- clustRviz:::is_positive_integer_scalar
  is_character_scalar <- clustRviz:::is_character_scalar
  is_nonempty_character_scalar <- clustRviz:::is_nonempty_character_scalar

  expect_true(is_logical_scalar(TRUE))
  expect_true(is_logical_scalar(FALSE))
  expect_false(is_logical_scalar(NA))
  expect_false(is_logical_scalar(0))
  expect_false(is_logical_scalar("a"))
  expect_false(is_logical_scalar(c(TRUE, TRUE)))

  expect_true(is_numeric_scalar(3))
  expect_true(is_numeric_scalar(3.5))
  expect_true(is_numeric_scalar(0))
  expect_true(is_numeric_scalar(-4))
  expect_false(is_numeric_scalar(NA))
  expect_false(is_numeric_scalar(c(2, 5)))
  expect_false(is_numeric_scalar("a"))

  expect_true(is_integer_scalar(3))
  expect_false(is_integer_scalar(3.5))
  expect_true(is_integer_scalar(0))
  expect_true(is_integer_scalar(-4))
  expect_false(is_integer_scalar(NA))
  expect_false(is_integer_scalar(c(2, 5)))
  expect_false(is_integer_scalar("a"))

  expect_true(is_percent_scalar(0.3))
  expect_true(is_percent_scalar(1))
  expect_true(is_percent_scalar(0))
  expect_false(is_percent_scalar(1.5))
  expect_false(is_percent_scalar(-1.5))
  expect_false(is_percent_scalar(NA))
  expect_false(is_percent_scalar(c(0.2, 0.5)))
  expect_false(is_percent_scalar("a"))

  expect_true(is_positive_scalar(0.3))
  expect_true(is_positive_scalar(1))
  expect_false(is_positive_scalar(0))
  expect_true(is_positive_scalar(1.5))
  expect_false(is_positive_scalar(-1.5))
  expect_false(is_positive_scalar(NA))
  expect_false(is_positive_scalar(c(0.2, 0.5)))
  expect_false(is_positive_scalar("a"))

  expect_true(is_positive_integer_scalar(3))
  expect_false(is_positive_integer_scalar(3.5))
  expect_false(is_positive_integer_scalar(0))
  expect_false(is_positive_integer_scalar(-4))
  expect_false(is_positive_integer_scalar(NA))
  expect_false(is_positive_integer_scalar(c(2, 5)))
  expect_false(is_positive_integer_scalar("a"))

  expect_false(is_character_scalar(3))
  expect_false(is_character_scalar(3.5))
  expect_false(is_character_scalar(0))
  expect_false(is_character_scalar(-4))
  expect_false(is_character_scalar(NA))
  expect_false(is_character_scalar(c(2, 5)))
  expect_true(is_character_scalar(""))
  expect_true(is_character_scalar("a"))

  expect_false(is_nonempty_character_scalar(3))
  expect_false(is_nonempty_character_scalar(3.5))
  expect_false(is_nonempty_character_scalar(0))
  expect_false(is_nonempty_character_scalar(-4))
  expect_false(is_nonempty_character_scalar(NA))
  expect_false(is_nonempty_character_scalar(c(2, 5)))
  expect_false(is_nonempty_character_scalar(""))
  expect_true(is_nonempty_character_scalar("a"))

  is_square <- clustRviz:::is_square
  expect_true(is_square(matrix(1, 5, 5)))
  expect_false(is_square(matrix(1, 5, 3)))
})

test_that("Capitalization works", {
  capitalize_string <- clustRviz:::capitalize_string

  expect_equal("A", capitalize_string("a"))
  expect_equal("Abc", capitalize_string("abc"))
  expect_equal("ABc", capitalize_string("ABc"))
  expect_equal("A Fantastic Cow", capitalize_string("a fantastic cow"))
})

test_that("Unscaling matrix works", {
  set.seed(5)
  n = 100; p = 400;

  X <- matrix(rnorm(n * p, sd = 1:25, mean = 1:50), ncol=p)
  X_std <- scale(X, center=TRUE, scale=TRUE)

  expect_equal(clustRviz:::unscale_matrix(X_std), X, check.attributes = FALSE)
})

test_that("Converting plot dimensions works", {
  set.seed(254)
  convert_units <- clustRviz:::convert_units
  x <- rexp(25, 3)

  for (from_unit in c("in", "cm", "mm", "px")) {
    for (to_unit in c("in", "cm", "mm", "px")) {
      round_trip <- convert_units(convert_units(x, from = from_unit, to = to_unit),
                                  from = to_unit,
                                  to = from_unit)
      expect_equal(x, round_trip)
    }
  }

  expect_equal(10,     convert_units(1, from = "cm", to = "mm"))
  expect_equal(0.1,    convert_units(1, from = "mm", to = "cm"))
  expect_equal(2.54,   convert_units(1, from = "in", to = "cm"))
  expect_equal(1/2.54, convert_units(1, from = "cm", to = "in"))
})

test_that("ensure_gif works", {
  ensure_gif <- clustRviz:::ensure_gif

  expect_equal("plot.gif", ensure_gif("plot.gif"))
  expect_equal("~/plot.gif", ensure_gif("~/plot.gif"))
  expect_equal("/my/long/path/plot.gif", ensure_gif("/my/long/path/plot.gif"))

  expect_warning(ensure_gif("plot.jpg"))
  expect_no_warning(ensure_gif("plot.gif"))

  expect_equal("plot.gif", suppressWarnings(ensure_gif("plot.jpg")))
  expect_equal("~/plot.gif", suppressWarnings(ensure_gif("~/plot.jpg")))
  expect_equal("/my/long/path/plot.gif", suppressWarnings(ensure_gif("/my/long/path/plot.jpg")))
})

test_that("ensure_html works", {
  ensure_html <- clustRviz:::ensure_html

  expect_equal("plot.html", ensure_html("plot.html"))
  expect_equal("~/plot.html", ensure_html("~/plot.html"))
  expect_equal("/my/long/path/plot.html", ensure_html("/my/long/path/plot.html"))

  expect_warning(ensure_html("plot.jpg"))
  expect_no_warning(ensure_html("plot.html"))

  expect_equal("plot.html", suppressWarnings(ensure_html("plot.jpg")))
  expect_equal("~/plot.html", suppressWarnings(ensure_html("~/plot.jpg")))
  expect_equal("/my/long/path/plot.html", suppressWarnings(ensure_html("/my/long/path/plot.jpg")))
})

test_that("connectedness check works", {
  is_connected_adj_mat <- clustRviz:::is_connected_adj_mat

  eye <- function(n) diag(1, nrow = n, ncol = n)

  expect_true(is_connected_adj_mat(eye(1)))
  expect_false(is_connected_adj_mat(eye(5)))

  A <- eye(3); A[1,2] <- A[2,3] <- A[2,1] <- A[3,2] <- 1
  expect_true(is_connected_adj_mat(A))

  A <- eye(3); A[1,2] <- A[2,1] <- 1
  expect_false(is_connected_adj_mat(A))
})

test_that("U smoothing for CARP works", {
  set.seed(200)

  N <- 50
  P <- 30

  U <- array(rnorm(N * P), c(N, P, 1))

  # Fake cluster assignments
  K <- 5
  membership <- sample(K, N, replace = TRUE)
  cluster_info <- list(membership = membership,
                       csize      = table(membership),
                       no         = length(unique(membership)))

  U_smoothed <- smooth_u_clustering(U, list(cluster_info))

  for(k in 1:K){
    u_row_mean <- colMeans(U[membership == k,,1])
    for(n in 1:N){
      if(membership[n] == k){
        expect_equal(U_smoothed[n,,1], u_row_mean)
      }
    }
  }
})
