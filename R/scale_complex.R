#' @export
#' @author R Core (\code{base} package)
#
# This function is basically just scale.default, but with is.numeric replaced by
# (is.numeric | is.complex) so that it scales complex matrices correctly and returns
# the appropriate scale factors as attributes.
scale.complex <- function(x, center = TRUE, scale = TRUE){
  is.nc <- function(x) is.numeric(x) | is.complex(x)

  x <- as.matrix(x)
  nc <- ncol(x)
  if (is.logical(center)) {
    if (center) {
      center <- colMeans(x, na.rm=TRUE)
      x <- sweep(x, 2L, center, check.margin=FALSE)
    }
  }
  else {
    if(!is.nc(center)) center <- as.complex(center)
    if (length(center) == nc)
      x <- sweep(x, 2L, center, check.margin=FALSE)
    else
      stop("length of 'center' must equal the number of columns of 'x'")
  }
  if (is.logical(scale)) {
    if (scale) {
      f <- function(v) {
        v <- v[!is.na(v)]
        # Use squared modulus in SD calc here rather than just squaring
        sqrt(sum(abs(v)^2) / max(1, length(v) - 1L))
      }
      scale <- apply(x, 2L, f)
      x <- sweep(x, 2L, scale, "/", check.margin=FALSE)
    }
  }
  else {
    if(!is.nc(scale)) scale <- as.complex(scale)
    if (length(scale) == nc)
      x <- sweep(x, 2L, scale, "/", check.margin=FALSE)
    else
      stop("length of 'scale' must equal the number of columns of 'x'")
  }
  if(is.nc(center)) attr(x, "scaled:center") <- center
  if(is.nc(scale )) attr(x, "scaled:scale" ) <- scale
  x
}
