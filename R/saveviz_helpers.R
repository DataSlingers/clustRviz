# Helper functions for saving

#' @noRd
#' Convert from one set of plotting units to another
convert_units <- function(units, from = c("in", "cm", "mm", "px"),
                                 to   = c("in", "cm", "mm", "px"),
                          dpi = 300){
  from <- match.arg(from)
  to   <- match.arg(to)

  switch(from,
         `in` = switch(to,
                       `in` = units,
                        cm  = units * 2.54,
                        mm  = units * 25.4,
                        px  = units * dpi),
         cm   = switch(to,
                       `in` = units / 2.54,
                       cm   = units,
                       mm   = units * 10,
                       px   = units * dpi / 2.54),
         mm   = switch(to,
                       `in` = units / 25.4,
                       cm   = units / 10,
                       mm   = units,
                       px   = units * dpi / 25.4),
         px   = switch(to,
                       `in` = units / dpi,
                       cm   = units / dpi * 2.54,
                       mm   = units / dpi * 25.4,
                       px   = units))
}

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

  width  <- convert_units(width,  from = units, to = "in")
  height <- convert_units(height, from = units, to = "in")

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

#' Convert a file name to ".gif" and warn if necessary
#' @noRd
#' @importFrom tools file_ext file_path_sans_ext
ensure_gif <- function(file_name){

  current_extension <- file_ext(file_name)

  if (current_extension != "gif") {
    new_file_name <- paste(file_path_sans_ext(file_name), ".gif", sep = "")
    crv_warning("File extension for dynamic visualization is not ", sQuote(".gif"),
                " -- changing from ", sQuote(file_name), " to ", sQuote(new_file_name))
    new_file_name
  } else {
    file_name
  }
}

ensure_html <- function(file_name){

  current_extension <- file_ext(file_name)

  if (current_extension != "html") {
    new_file_name <- paste(file_path_sans_ext(file_name), ".html", sep = "")
    crv_warning("File extension for dynamic visualization is not ", sQuote(".html"),
                " -- changing from ", sQuote(file_name), " to ", sQuote(new_file_name))
    new_file_name
  } else {
    file_name
  }
}
