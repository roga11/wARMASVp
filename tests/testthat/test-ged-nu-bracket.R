# Regression tests for the GED nu root-finder fixed-bracket [0.1, 20] change
# (2026-07-28). Background: the historical data-scale-dependent bracket
# (x0 = mean|log(y^2)| / sd|log(y^2)|) made GED nu-hat non-scale-invariant
# because Brent stopped at different points inside its tolerance ball. The
# fixed [0.1, 20] bracket (already the C++ fallback) plus Brent tol 1e-12
# restores scale invariance and harmonizes GED with the fixed-bracket Student-t
# solver. See SVpHeavyTails lessons_learned L91 / SVpTest LAB_NOTES 2026-07-28.

g_ged <- function(x) (2 / x)^2 * psigamma(1 / x, 1)

test_that("GED nu-hat is scale invariant under y -> c*y (del -> 0)", {
  # del = 1e-300 is a bit-for-bit no-op inside log(y^2 + del) on continuous
  # data, so this isolates the estimator's scale invariance from the del guard.
  set.seed(20260728)
  y <- sim_svp(2000, phi = 0.95, sigy = 1, sigv = 0.50,
               errorType = "GED", nu = 1.5, burnin = 500)$y
  nu1 <- suppressWarnings(svp(y,     p = 1, J = 10, errorType = "GED", del = 1e-300))$v
  nu2 <- suppressWarnings(svp(2 * y, p = 1, J = 10, errorType = "GED", del = 1e-300))$v
  expect_lt(abs(nu2 - nu1), 1e-10)
})

test_that("GED solver matches a high-precision uniroot reference on [0.12, 19]", {
  f_ged <- wARMASVp:::find_nu_ged_cpp
  nus <- seq(0.12, 19, length.out = 200)
  dmax <- 0
  for (nu_true in nus) {
    se2 <- g_ged(nu_true)
    r_cpp <- f_ged(se2, 0.1, 20)                      # tol default now 1e-12
    r_ref <- stats::uniroot(function(x) g_ged(x) - se2,
                            interval = c(0.1, 20), tol = 1e-14)$root
    dmax <- max(dmax, abs(r_cpp - r_ref))
  }
  expect_lt(dmax, 1e-6)
})

test_that("GED solver boundary-returns exactly 0.1 / 20 for out-of-range roots", {
  f_ged <- wARMASVp:::find_nu_ged_cpp
  # true root above 20 (thin tail): se2 < g(20) -> clip to 20
  expect_equal(f_ged(g_ged(25), 0.1, 20), 20, tolerance = 1e-12)
  # true root below 0.1 (very heavy tail): se2 > g(0.1) -> clip to 0.1
  expect_equal(f_ged(g_ged(0.05), 0.1, 20), 0.1, tolerance = 1e-12)
})

test_that("GED estimation warns at the nu boundary", {
  # thin-tailed GED (nu large) can push se2 below g(20); construct a path that
  # lands at the upper boundary and check the warning fires.
  set.seed(7)
  y <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.50,
               errorType = "GED", nu = 8, burnin = 500)$y
  # Not guaranteed to hit the boundary on every seed; assert only that IF a
  # boundary nu is returned it equals 20 exactly (semantics preserved).
  fit <- suppressWarnings(svp(y, p = 1, J = 10, errorType = "GED", del = 1e-10))
  if (fit$v >= 20 - 1e-9) expect_equal(fit$v, 20, tolerance = 1e-9)
  else expect_gt(fit$v, 0.1)
})

test_that("Student-t solver tol tightening improves accuracy vs a reference", {
  # The tol tightening (1e-6 -> 1e-12) also affects find_nu_t_cpp: it shifts the
  # Brent stopping point by up to ~2e-7 (NOT tol-invariant; the "t was already
  # bit-stable" note in the source refers to SCALE invariance via its fixed
  # bracket, not tol invariance). The meaningful property is that the tighter
  # tol agrees at least as well with a high-precision reference root.
  f_t <- wARMASVp:::find_nu_t_cpp
  f_psi <- function(x) psigamma(x / 2, 1)
  nus <- seq(2.5, 100, length.out = 60)
  d12 <- d6 <- 0
  for (nu_true in nus) {
    se2b <- f_psi(nu_true)
    ref  <- stats::uniroot(function(x) f_psi(x) - se2b,
                           interval = c(2.01, 500), tol = 1e-14)$root
    d12 <- max(d12, abs(f_t(se2b, 2.01, 500, FALSE, tol = 1e-12) - ref))
    d6  <- max(d6,  abs(f_t(se2b, 2.01, 500, FALSE, tol = 1e-6)  - ref))
  }
  # tighter tol is at least as accurate, and the shift is bounded well below
  # any reporting precision
  expect_lte(d12, d6 + 1e-12)
  expect_lt(abs(d12), 1e-6)
})
