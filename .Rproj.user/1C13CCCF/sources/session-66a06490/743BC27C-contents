.onLoad <- function(libname, pkgname) {
  # gsignal must be on the search path because the C++ code accesses

  # gsignal::poly via Rcpp::Environment("package:gsignal")
  if (!("package:gsignal" %in% search())) {
    suppressPackageStartupMessages(
      requireNamespace("gsignal", quietly = TRUE)
    )
    # Attach gsignal to the search path so C++ can find it
    attachNamespace("gsignal")
  }
}

.onUnload <- function(libpath) {
  library.dynam.unload("wARMASVp", libpath)
}
