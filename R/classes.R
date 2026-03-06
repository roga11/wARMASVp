# S3 class constructors and validators for wARMASVp
# Classes are created by their respective user-facing function svp()
# This file contains any shared validation logic.

#' @keywords internal
validate_svp <- function(x) {
  stopifnot(is.list(x))
  stopifnot("phi" %in% names(x))
  stopifnot("sigv" %in% names(x))
  stopifnot("sigy" %in% names(x))
  stopifnot("mu" %in% names(x))
  invisible(x)
}
