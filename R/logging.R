# Logging infrastructure for clustRviz

## This must be kept consistent with src/clustRviz_logging.h::ClustRVizLoggerLevel
LEVELS <- c(ERROR    = 40,
            WARNING  = 30,
            MESSAGE  = 20,
            INFO     = 10,
            DEBUG    = 00)


#' \code{ClustRViz} Package Logging Functionality
#'
#' Control the global logging level for the \code{clustRviz} package.
#'
#' @importFrom stats setNames
#' @export
#' @param level The desired new log level. Available levels are \itemize{
#'     \item \code{ERROR} - corresponding to \code{base::stop};
#'     \item \code{WARNING} - corresponding to \code{base::warning};
#'     \item \code{MESSAGE} - corresponding to \code{base::message};
#'     \item \code{INFO}; and
#'     \item \code{DEBUG.}
#' } If omitted, the log level is not changed (and the current level is still
#' returned invisibly.) See below for details about the different levels.
#' @return The previous log level (invisibly).
#' @rdname clustRviz_logging
#' @aliases clustRviz_logging clustRviz_logger_level
#' @details The \code{clustRviz} package has a multi-level logging system, with a single
#'     global log level; (which applies to both \code{R} and \code{C++} level
#'     functionality.) the levels are, in decreasing order, \code{ERROR},
#'     \code{WARNING}, \code{MESSAGE} (default), \code{INFO}, \code{DEBUG}.
#'
#'     To change the amount of output from the \code{clustRviz} package, the
#'     \code{clustRviz_logger_level} function can be used to adjust the global
#'     log level. The \code{INFO} and \code{DEBUG} levels can be quite verbose
#'     and may significantly slow down the package.
#' @examples
#' # Switch to INFO level and fit somewhat loudly
#' clustRviz_logger_level("INFO")
#' CARP(presidential_speech)
#'
#' # Return to default behavior
#' clustRviz_logger_level("MESSAGE")
#' # Will give a few messages about phases of fitting
#' CARP(presidential_speech)
#'
#' # If we switch to ERRORs only, no progress message
#' clustRviz_logger_level("ERROR")
#' CARP(presidential_speech)
#'
#' # Return to default
#' clustRviz_logger_level("MESSAGE")
#'
#' # Normally this prints messages as the paths are calculated
#' CARP(presidential_speech)
#' # But it can be suppressed using standard R functions
#' suppressMessages(CARP(presidential_speech))
clustRviz_logger_level <- function(level=c("ERROR",
                                           "WARNING",
                                           "MESSAGE",
                                           "INFO",
                                           "DEBUG")){

    LEVELS_REV <- setNames(names(LEVELS), LEVELS)

    old_level <- LEVELS_REV[as.character(clustRviz_get_logger_level_cpp())]
    names(old_level) <- NULL

    if(!missing(level)){
        level <- match.arg(level);
        clustRviz_set_logger_level_cpp(LEVELS[level])
        return(invisible(old_level))
    }

    old_level
}

crv_error <- function(..., call=TRUE){
    msg <- paste(list(...), collapse="")

    ## Try to add R level calling info
    if(identical(call, TRUE)){
        tryCatch({
            msg <- paste0(msg, " (Called from ", as.character(as.list(sys.call(-1))[[1]]), ")")
        }, error=function(e){})
    } else if(is.character(call)){
        msg <- paste0(msg, " (Called from ", call, ")")
    }

    clustRviz_log_cpp(LEVELS["ERROR"], msg)
}

crv_warning <- function(..., call=TRUE){
    msg <- paste(list(...), collapse="")

    ## Try to add R level calling info
    if(identical(call, TRUE)){
        tryCatch({
            msg <- paste0(msg, " (Called from ", as.character(as.list(sys.call(-1))[[1]]), ")")
        }, error=function(e){})
    } else if(is.character(call)){
        msg <- paste0(msg, " (Called from ", call, ")")
    }

    clustRviz_log_cpp(LEVELS["WARNING"], msg)
}

crv_message <- function(...){
    msg <- paste(list(...), collapse="")
    clustRviz_log_cpp(LEVELS["MESSAGE"], msg)
}

crv_info <- function(...){
    msg <- paste(list(...), collapse="")
    clustRviz_log_cpp(LEVELS["INFO"], msg)
}

crv_debug <- function(...){
    msg <- paste(list(...), collapse="")
    clustRviz_log_cpp(LEVELS["DEBUG"], msg)
}
