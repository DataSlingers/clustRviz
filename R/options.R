## clustRviz options

clustRviz_default_options <- list(rho                = 1.0,
                                  max_iter           = as.integer(5e6),
                                  burn_in            = 50L,
                                  viz_initial_step   = 1.1,
                                  viz_small_step     = 1.01,
                                  viz_max_inner_iter = 15L,
                                  keep               = 10L,
                                  epsilon            = 0.000001,
                                  keep_debug_info    = FALSE)

.clustRvizOptionsEnv <- list2env(clustRviz_default_options)

#' \code{ClustRViz} Options
#'
#' Advanced control of algorithmic options for \code{\link{CARP}} and \code{\link{CBASS}}.
#' The \code{clustRviz_reset_options} function returns options to "factory-fresh"
#' settings.
#'
#' @param ... Options (to be passed by name). See below for available options.
#'
#' @details The following options can be set by name:\itemize{
#'   \item \code{epsilon} The initial step size (fixed during the "burn-in" period)
#'   \item \code{max_iter} An integer: the maximum number of iterations to perform.
#'   \item \code{burn_in} An integer: the number of initial iterations at a fixed
#'                       (small) value of \eqn{\gamma}
#'   \item \code{viz_initial_step} The initial (large) step size used in back-tracking
#'                                 (\code{CARP-VIZ} and \code{CBASS-VIZ}) algorithms.
#'   \item \code{viz_small_step} The secondary (small) step size used in back-tracking
#'                               (\code{CARP-VIZ} and \code{CBASS-VIZ}) algorithms.
#'   \item \code{viz_max_inner_iter} The maximum number of iterations to perform
#'                                   in the inner loop of back-tracking (\code{CARP-VIZ}
#'                                   and \code{CBASS-VIZ}) algorithms.
#'   \item \code{keep} \code{\link{CARP}} and \code{\link{CBASS}} keep every
#'                     \code{keep}-th iteration even if no fusions are detected.
#'                     Increasing this parameter may improve performance, at
#'                     the expense of returning a finer grid.
#'   \item \code{rho} For advanced users only (not advisable to change): the penalty
#'                    parameter used for the augmented Lagrangian.
#'   \item \code{keep_debug_info}: Should additional debug info (currently only the V-path)
#'                                 be kept?
#' }
#' @rdname options
#' @export
clustRviz_options <- function(...){
  dots <- list(...)

  if (length(dots) == 0){
    return(as.list(.clustRvizOptionsEnv))
  }

  if (is.list(dots[[1]])) {
    dots <- do.call(c, dots)
  }

  if ( (is.null(names(dots))) || (any(names(dots) == "")) ){
    crv_error("All arguments to ", sQuote("clustRviz_options"), " must be named.")
  }

  known_names <- names(.clustRvizOptionsEnv)

  if ( any(names(dots) %not.in% known_names) ){
    unknown_names <- which(names(dots) %not.in% known_names)
    crv_error("Unknown argument ", sQuote(names(dots)[unknown_names[1]]), " passed to ", sQuote("clustRviz_options."))
  }

  old_opts <- as.list(.clustRvizOptionsEnv)

  for(ix in seq_along(dots)){
    nm  <- names(dots)[ix]
    opt <- dots[[ix]]

    ## Validate
    if (nm %in% c("rho", "epsilon")) {
      if (!is_positive_scalar(opt)) {
        crv_error(sQuote(nm), " must be a positive scalar.")
      }
    } else if (nm %in% c("viz_initial_step", "viz_small_step")) {
      if ( (!is_positive_scalar(opt)) || (opt <= 1) ){
        crv_error(sQuote(nm), " must be greater than one.")
      }
    } else if (nm %in% c("burn_in", "max_iter", "viz_burn_in", "viz_max_inner_iter", "keep")) {
      if (!is_positive_integer_scalar(opt) ){
        crv_error(sQuote(nm), " must be a positive integer.")
      }
    } else if (nm %in% "keep_debug_info") {
      if (!is_logical_scalar(opt)) {
        crv_error(sQuote(nm), " must be a logical scalar.")
      }
    }

    ## Assign
    assign(nm, opt, .clustRvizOptionsEnv)
  }

  ## Sanity checks
  if (.clustRvizOptionsEnv[["burn_in"]] + 100 >= .clustRvizOptionsEnv[["max_iter"]]) {
    crv_warning(sQuote("burn_in"), " should typically be at least 100 less than ", sQuote("max_iter."))
  }

  if (.clustRvizOptionsEnv[["viz_small_step"]] >= .clustRvizOptionsEnv[["viz_initial_step"]]) {
    crv_warning(sQuote("viz_small_step"), " should be less than ", sQuote("viz_initial_step."))
  }

  ##
  invisible(old_opts)
}

#' @rdname options
#' @export
clustRviz_reset_options <- function(){
  do.call(clustRviz_options, clustRviz_default_options)
}
