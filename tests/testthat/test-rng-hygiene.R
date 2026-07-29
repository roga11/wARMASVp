# RNG-hygiene regression tests (2026-07-29). Package computations must never
# reset the caller's ambient .Random.seed as a side effect, while their own
# outputs must stay reproducible via their internal `seed` argument. Guards
# against the fit_ksc_mixture / bootstrap-particle-filter RNG-clobber bug
# (SVpTest LAB_NOTES 2026-07-28 "THE BUG"; mirrors SVpTest test P10).

make_fit <- function(errorType, nu = NULL) {
  set.seed(99)
  y <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.5,
               errorType = errorType, nu = nu, burnin = 200)$y
  suppressWarnings(svp(y, p = 1, errorType = errorType, del = 1e-300))
}

fits <- list(
  Gaussian    = make_fit("Gaussian"),
  `Student-t` = make_fit("Student-t", nu = 5),
  GED         = make_fit("GED", nu = 1.5)
)
methods <- c("corrected", "mixture", "particle")

test_that("filter evaluations do not reset the caller's ambient RNG stream", {
  draw_after <- function(fit, method, pre) {
    set.seed(pre)
    invisible(suppressWarnings(
      filter_svp(fit, method = method, K = 7, M = 200, seed = 42)$loglik))
    runif(3)
  }
  for (et in names(fits)) {
    for (m in methods) {
      x <- draw_after(fits[[et]], m, 111)
      y <- draw_after(fits[[et]], m, 222)
      expect_false(identical(x, y),
        info = sprintf("%s / %s: ambient RNG was reset", et, m))
    }
  }
})

test_that("filter log-likelihood is independent of the ambient RNG state", {
  # The internal `seed` must keep the returned value reproducible regardless of
  # the caller's stream (this is what makes the save/restore safe).
  for (et in names(fits)) {
    for (m in methods) {
      set.seed(111)
      l1 <- suppressWarnings(filter_svp(fits[[et]], method = m, K = 7,
                                        M = 200, seed = 42)$loglik)
      set.seed(222)
      l2 <- suppressWarnings(filter_svp(fits[[et]], method = m, K = 7,
                                        M = 200, seed = 42)$loglik)
      expect_equal(l1, l2,
        info = sprintf("%s / %s: loglik depends on ambient seed", et, m))
    }
  }
})

test_that("svp_IC / forecast_svp also leave the ambient RNG stream intact", {
  draw_after <- function(expr_fun, pre) {
    set.seed(pre); invisible(suppressWarnings(expr_fun())); runif(3)
  }
  fit <- fits[["GED"]]
  ic_call <- function() svp_IC(fit, criteria = "BIC_Kalman",
                               filter_method = "mixture", seed = 42)
  fc_call <- function() forecast_svp(fit, H = 2, filter_method = "mixture",
                                     seed = 42)
  expect_false(identical(draw_after(ic_call, 1), draw_after(ic_call, 2)))
  expect_false(identical(draw_after(fc_call, 1), draw_after(fc_call, 2)))
})

test_that(".with_preserved_seed restores a pre-existing seed exactly", {
  set.seed(123)
  s_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  wARMASVp:::.with_preserved_seed({ set.seed(42); rnorm(10) })
  s_after <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  expect_identical(s_before, s_after)
})

test_that("ambient RNG is preserved under RNGkind('L'Ecuyer-CMRG')", {
  # L'Ecuyer-CMRG is the RNG kind used by the cluster / mclapply SV simulations,
  # and its .Random.seed is a longer integer vector than the default Mersenne
  # Twister. Verify save/restore handles it byte-for-byte and leaves RNGkind
  # untouched -- both via the helper and through the public filter API (which is
  # what those simulations call between simulated draws).
  old_kind <- RNGkind()
  on.exit(suppressWarnings(do.call(RNGkind, as.list(old_kind))), add = TRUE)
  RNGkind("L'Ecuyer-CMRG")

  # helper directly
  set.seed(7)
  s_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  wARMASVp:::.with_preserved_seed({ set.seed(42); rnorm(100) })
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE), s_before)
  expect_identical(RNGkind()[1], "L'Ecuyer-CMRG")

  # through the public filter API (mixture + particle)
  set.seed(7)
  s_api <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  invisible(suppressWarnings(
    filter_svp(fits[["GED"]], method = "mixture", K = 7, seed = 42)$loglik))
  invisible(suppressWarnings(
    filter_svp(fits[["GED"]], method = "particle", M = 200, seed = 42)$loglik))
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE), s_api)
  expect_identical(RNGkind()[1], "L'Ecuyer-CMRG")
})

test_that("fresh session: internal seeding leaves no .Random.seed behind", {
  # Remove-if-absent branch of .with_preserved_seed(): if the caller had never
  # drawn a random number (.Random.seed absent), the internal set.seed() must
  # not leave a deterministic RNG state behind. Tested at the helper /
  # fit_ksc_mixture level. NOTE: this property cannot be asserted through
  # filter_svp(), because Rcpp's automatic RNGScope re-materialises
  # .Random.seed (with an entropy-based value, not a deterministic reset) on
  # every C++ kernel return -- standard Rcpp behaviour, and confirmed to occur
  # even for the RNG-free "corrected" filter. What matters for the bug is that
  # an *existing* stream is never reset to a fixed state (tested above).
  ge <- globalenv()

  # the helper directly
  if (exists(".Random.seed", envir = ge, inherits = FALSE)) rm(".Random.seed", envir = ge)
  invisible(wARMASVp:::.with_preserved_seed({ set.seed(42); rnorm(5) }))
  expect_false(exists(".Random.seed", envir = ge, inherits = FALSE))

  # fit_ksc_mixture (GED reaches the internal sampling path)
  if (exists(".Random.seed", envir = ge, inherits = FALSE)) rm(".Random.seed", envir = ge)
  invisible(wARMASVp:::fit_ksc_mixture("ged", nu = 1.5, K = 7, seed = 42))
  expect_false(exists(".Random.seed", envir = ge, inherits = FALSE))
})
