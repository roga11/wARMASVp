# --------------------------------------------------------------------------- #
# Precomputed GMKF mixture tables and their run-time interpolation
#
# The tables in inst/extdata/ store, for each error family and each
# K in {7, 10, 12, 16}, one Gaussian mixture per nu node approximating the
# density of the centred measurement noise eps = log(z^2) - E[log(z^2)]
# (Gaussian: a single mixture, no shape parameter; GED: 200 nodes on
# [0.1, 20]; Student-t: 500 nodes on [2.01, 500] -- the estimators' nu
# clamps, so a fitted model can never fall outside a table). The build
# scripts and a replication guide ship in inst/gmkf-tables/.
#
# Interpolation between nodes runs cubic splines through the stored
# component paths, in the coordinate the grid was laid out in (log(nu) for
# GED, 1/nu for Student-t) and in the parameterisation the tables were
# fitted in: log(w_k / w_1) for k = 2..K, m_k, log(v_k). Any point
# interpolated there maps back to positive weights summing to one and
# positive variances, so the result is always a proper mixture. The
# component columns are stored in matched order across neighbouring nodes,
# which is what makes splining down a column meaningful.
#
# Ported from GMKFTable/Table_GMKF_interpolator.R
# (md5 a5b870e571be060e44bcc3c21c169dc3 at port time, 2026-08-05).
# --------------------------------------------------------------------------- #

# Session cache: raw tables and built interpolants, loaded/built once each.
.gmkf_cache <- new.env(parent = emptyenv())

#' Access a shipped GMKF mixture table
#'
#' Returns the precomputed Gaussian-mixture table that
#' \code{\link{filter_svp}} (and everything downstream of it) uses when
#' \code{mixture = "optimal"}. Each table stores, for every mixture size
#' \code{K} in \{7, 10, 12, 16\}, the fitted component weights, means and
#' variances -- one mixture for the Gaussian family, and one per \eqn{\nu}
#' node for GED (200 nodes on [0.1, 20]) and Student-t (500 nodes on
#' [2.01, 500]). The mixtures approximate the density of the centred
#' measurement noise \eqn{\log(z^2) - E[\log(z^2)]} and were fitted by
#' minimising Kullback-Leibler divergence against the exact density
#' (quadrature, not sampling).
#'
#' The returned object carries a \code{"provenance"} attribute recording the
#' md5 of the build script and helpers that produced it; the scripts
#' themselves ship in \code{system.file("gmkf-tables", package = "wARMASVp")}
#' together with a README describing how to rebuild and verify the tables.
#'
#' @param family One of \code{"gaussian"}, \code{"ged"}, \code{"studentt"}.
#'
#' @return A named list with one entry per mixture size (names \code{"7"},
#'   \code{"10"}, \code{"12"}, \code{"16"}). For \code{"gaussian"} each entry
#'   holds vectors \code{w}, \code{m}, \code{v} (and build diagnostics); for
#'   \code{"ged"} and \code{"studentt"} each entry holds the node grid
#'   \code{nu} and matrices \code{w}, \code{m}, \code{v} with one row per
#'   node.
#'
#' @examples
#' tab <- gmkf_table("gaussian")
#' names(tab)
#' attr(tab, "provenance")
#'
#' @export
gmkf_table <- function(family = c("gaussian", "ged", "studentt")) {
  family <- match.arg(family)
  key <- paste0("table_", family)
  if (is.null(.gmkf_cache[[key]])) {
    path <- system.file("extdata", paste0(family, "_mixtures.rds"),
                        package = "wARMASVp")
    if (!nzchar(path)) {
      stop("GMKF table for '", family, "' not found in the installed package.")
    }
    .gmkf_cache[[key]] <- readRDS(path)
  }
  .gmkf_cache[[key]]
}

# Build the spline interpolant for one table (internal; see file header).
# `nodes` restricts the build to a subset of node indices (diagnostics only).
#' @importFrom stats splinefun
#' @noRd
gmkf_interpolant <- function(tab, family = c("ged", "studentt", "gaussian"),
                             method = "fmm", nodes = NULL) {
  family <- match.arg(family)
  if (family == "gaussian")
    return(structure(list(family = family, entries = tab),
                     class = "gmkf_interpolant"))
  coord <- switch(family, ged = log, studentt = function(nu) 1 / nu)
  entries <- lapply(tab, function(tb) {
    K    <- tb$K
    keep <- if (is.null(nodes)) seq_along(tb$nu) else sort(unique(as.integer(nodes)))
    t    <- coord(tb$nu[keep])
    o    <- order(t); t <- t[o]; keep <- keep[o]
    W <- tb$w[keep, , drop = FALSE]
    M <- tb$m[keep, , drop = FALSE]
    V <- tb$v[keep, , drop = FALSE]
    Y <- cbind(log(W[, -1, drop = FALSE] / W[, 1]), M, log(V))
    list(K = K, t_range = range(t),
         sp = lapply(seq_len(ncol(Y)),
                     function(j) splinefun(t, Y[, j], method = method)))
  })
  structure(list(family = family, coord = coord, entries = entries,
                 nu_range = range(tab[[1]]$nu[if (is.null(nodes))
                   seq_along(tab[[1]]$nu) else as.integer(nodes)])),
            class = "gmkf_interpolant")
}

# Evaluate an interpolant at nu for mixture size K (internal).
# Returns list(w, m, v) -- a proper mixture.
#' @noRd
gmkf_mixture <- function(itp, nu = NULL, K) {
  e <- itp$entries[[as.character(K)]]
  if (is.null(e))
    stop("K = ", K, " is not in this table (available: ",
         paste(names(itp$entries), collapse = ", "), ")")
  if (itp$family == "gaussian")
    return(list(w = e$w, m = e$m, v = e$v))
  if (is.null(nu)) stop("nu is required for the ", itp$family, " table")
  t <- itp$coord(nu)
  lo <- e$t_range[1]; hi <- e$t_range[2]; tol <- 1e-12 * (hi - lo)
  if (t < lo - tol || t > hi + tol)
    stop(sprintf("nu = %g is outside the table's range [%g, %g]",
                 nu, itp$nu_range[1], itp$nu_range[2]))
  t <- min(max(t, lo), hi)
  y <- vapply(e$sp, function(f) f(t), numeric(1))
  Kn  <- e$K
  eta <- c(0, y[seq_len(Kn - 1)])
  ew  <- exp(eta - max(eta))
  list(w = ew / sum(ew),
       m = y[Kn:(2 * Kn - 1)],
       v = exp(y[(2 * Kn):(3 * Kn - 1)]))
}

# Shared validation of the (mixture, K, distribution) combination. Single
# source of the error messages; called by filter_svp() for method = "mixture"
# and by svp_IC() before its tryCatch (only when a Kalman criterion is
# requested with filter_method = "mixture"), so an invalid combination errors
# loudly instead of surfacing as silent NA information criteria.
# `dist_name` uses the .get_filter_params vocabulary:
# "gaussian" / "student_t" / "ged".
.validate_mixture_choice <- function(mixture, K, dist_name) {
  if (mixture == "optimal") {
    if (!(is.numeric(K) && length(K) == 1 && K == as.integer(K) &&
          as.integer(K) %in% c(7L, 10L, 12L, 16L))) {
      stop("GMKF tables exist for K in {7, 10, 12, 16}; use one of these ",
           "or mixture = \"em\" (any K).")
    }
  } else if (mixture == "ksc") {
    if (dist_name != "gaussian" || !(is.numeric(K) && length(K) == 1 && K == 7)) {
      stop("mixture = \"ksc\" is the published KSC (1998) table, available ",
           "only for Gaussian errors with K = 7; use mixture = \"optimal\" ",
           "or \"em\".")
    }
  }
  invisible(TRUE)
}

# Mixture from the precomputed tables, in the fit_ksc_mixture() return shape
# (minus KL_div). `distribution` uses the package codes
# ("gaussian" / "student_t" / "ged"); the table-family spelling ("studentt")
# stays internal to this function, and the returned `distribution` echoes the
# package code so filter output looks the same under "optimal" and "em".
.gmkf_table_mixture <- function(distribution, nu = NULL, K = 7) {
  .validate_mixture_choice("optimal", K, distribution)
  family <- switch(distribution,
                   gaussian  = "gaussian",
                   student_t = "studentt",
                   ged       = "ged",
                   stop("unknown distribution '", distribution, "'"))
  key <- paste0("interp_", family)
  if (is.null(.gmkf_cache[[key]])) {
    .gmkf_cache[[key]] <- gmkf_interpolant(gmkf_table(family), family)
  }
  g <- gmkf_mixture(.gmkf_cache[[key]],
                    nu = if (family == "gaussian") NULL else nu, K = K)
  list(weights = g$w,
       means   = g$m,
       vars    = g$v,
       nu      = if (family == "gaussian") NULL else nu,
       distribution = distribution)
}
