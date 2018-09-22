# Helper functions for saving

#' Load a new graphics device for plotting static base graphics visualizations
#' @noRd
#' @param file.name the file name of the graphics device
#' @param width the width (in \code{units}) of the graphics device
#' @param height the height (in \code{units}) of the graphics device
#' @param units the units
#'
#' This code is loosely based on the implementation of ggplot2::ggsave
#' @importFrom tools file_ext
#' @importFrom grDevices postscript jpeg png tiff bmp svg pictex pdf
crv_new_dev_static <- function(file.name, width, height, units = c("in", "cm", "mm", "px")){

  units <- match.arg(units)
  ## Convert to inches for consistency
  if(units == "cm"){
    width  <- width / 2.54
    height <- height / 2.54
  } else if (units == "mm") {
    width  <- width / 25.4
    height <- height / 25.4
  } else if (units == "px") {
    width  <- width / 300 # Assume default 300 dpi
    height <- height / 300
  }

  ext <- file_ext(file.name)

  dev_func <- switch(ext,
                     ps   = function(filename, ...) postscript(file = filename, ..., onefile = FALSE, horizontal = FALSE, paper = "special"),
                     pdf  = function(filename, ...) pdf(file = filename, ..., version = 1.4),
                     jpeg = function(filename, ...) jpeg(filename, ..., units = "in", res = 300) ,
                     jpg  = function(filename, ...) jpeg(filename, ..., units = "in", res = 300),
                     png  = function(filename, ...) png(filename, ..., units = "in", res = 300),
                     tiff = function(filename, ...) tiff(filename, ..., units = "in", res = 300),
                     bmp  = function(filename, ...) bmp(filename, ..., units = "in", res = 300),
                     svg  = function(filename, ...) svg(filename, ...),
                     tex  = function(filename, ...) pictex(file = filename, ...),
                     crv_error("File extension ", sQuote(ext), " not currently supported by ", sQuote("saveviz."))
                     )

  dev_func(file.name, width = width, height = height)
}
