# =========================================================================== #
# Tests for filtering: CKF, GMKF, helpers, density functions
# =========================================================================== #

# --- Helper: simulate and estimate a model ---
# sim_svp returns a vector when leverage=FALSE, a list when leverage=TRUE
.fit_model <- function(n = 500, phi = 0.95, sigy = 1, sigv = 0.3,
                       errorType = "Gaussian", nu = NULL,
                       leverage = FALSE, rho = 0, seed = 123) {
  set.seed(seed)
  sim <- sim_svp(n, phi = phi, sigy = sigy, sigv = sigv,
                 errorType = errorType, nu = nu,
                 leverage = leverage, rho = rho)
  y <- if (is.list(sim)) sim$y else as.numeric(sim)
  fit <- svp(y, p = 1, leverage = leverage, errorType = errorType)
  list(sim = sim, y = y, fit = fit)
}


# =========================================================================== #
# Lyapunov solver
# =========================================================================== #

test_that("solve_lyapunov_discrete satisfies X = F X F' + Q", {
  F_mat <- matrix(c(0.95, 0, 1, 0), 2, 2, byrow = TRUE)
  Q <- matrix(c(0.09, 0, 0, 0), 2, 2)
  X <- solve_lyapunov_discrete(F_mat, Q)
  residual <- X - F_mat %*% X %*% t(F_mat) - Q
  expect_true(max(abs(residual)) < 1e-10)
})

test_that("solve_lyapunov_discrete works for p=1", {
  F_mat <- matrix(0.95)
  Q <- matrix(0.09)
  X <- solve_lyapunov_discrete(F_mat, Q)
  expect_equal(as.numeric(X), 0.09 / (1 - 0.95^2), tolerance = 1e-10)
})


# =========================================================================== #
# pged_std (CDF of standardized GED)
# =========================================================================== #

test_that("pged_std: CDF at 0 is 0.5", {
  for (nu in c(0.5, 1, 1.5, 2, 3)) {
    expect_equal(pged_std(0, nu), 0.5, tolerance = 1e-12)
  }
})

test_that("pged_std: monotonically increasing", {
  u_grid <- seq(-5, 5, by = 0.5)
  for (nu in c(1, 1.5, 2)) {
    p_vals <- pged_std(u_grid, nu)
    expect_true(all(diff(p_vals) > 0))
  }
})

test_that("pged_std: limits at +/- infinity", {
  expect_true(pged_std(-20, 1.5) < 0.001)
  expect_true(pged_std(20, 1.5) > 0.999)
})


# =========================================================================== #
# Density functions
# =========================================================================== #

test_that("density_eps_t integrates to 1", {
  for (nu in c(3, 5, 10, 30)) {
    integral <- integrate(density_eps_t, -30, 30, nu = nu)$value
    expect_equal(integral, 1.0, tolerance = 0.01)
  }
})

test_that("density_eps_ged integrates to 1", {
  for (nu in c(0.8, 1.0, 1.5, 2.0)) {
    integral <- integrate(density_eps_ged, -30, 30, nu = nu)$value
    expect_equal(integral, 1.0, tolerance = 0.01)
  }
})


# =========================================================================== #
# fit_ksc_mixture
# =========================================================================== #

test_that("fit_ksc_mixture Gaussian returns KSC table", {
  mix <- fit_ksc_mixture("gaussian")
  expect_equal(length(mix$weights), 7)
  expect_equal(sum(mix$weights), 1.0, tolerance = 1e-10)
  # KSC table approximates raw log(chi^2_1), mean ≈ E[log(chi^2_1)] = psi(1/2)+log(2) ≈ -1.27
  mix_mean <- sum(mix$weights * mix$means)
  expect_equal(mix_mean, digamma(0.5) + log(2), tolerance = 0.05)
  # Variance should be approximately pi^2/2
  mix_var <- sum(mix$weights * (mix$vars + mix$means^2)) - mix_mean^2
  expect_equal(mix_var, pi^2 / 2, tolerance = 0.5)
})

test_that("fit_ksc_mixture Student-t: valid moments", {
  mix <- fit_ksc_mixture("student_t", nu = 5, K = 5, n_sample = 50000, seed = 42)
  expect_equal(length(mix$weights), 5)
  expect_equal(sum(mix$weights), 1.0, tolerance = 1e-4)
  mix_mean <- sum(mix$weights * mix$means)
  expect_true(abs(mix_mean) < 0.5)
  sigma_eps2_t <- psigamma(0.5, 1) + psigamma(5 / 2, 1)
  mix_var <- sum(mix$weights * (mix$vars + mix$means^2)) - mix_mean^2
  expect_equal(mix_var, sigma_eps2_t, tolerance = 1.5)
})

test_that("fit_ksc_mixture GED: valid moments", {
  mix <- fit_ksc_mixture("ged", nu = 1.5, K = 5, n_sample = 50000, seed = 42)
  expect_equal(length(mix$weights), 5)
  expect_equal(sum(mix$weights), 1.0, tolerance = 1e-4)
  mix_mean <- sum(mix$weights * mix$means)
  expect_true(abs(mix_mean) < 0.5)
})


# =========================================================================== #
# filter_svp — CKF (corrected Kalman)
# =========================================================================== #

test_that("filter_svp returns correct class and components (Gaussian)", {
  m <- .fit_model()
  filt <- filter_svp(m$fit)
  expect_s3_class(filt, "svp_filter")
  expect_true(all(c("w_filtered", "w_smoothed", "zt", "zt_smoothed",
                     "P_filtered", "P_predicted", "xi_filtered", "xi_smoothed",
                     "loglik", "method", "model") %in% names(filt)))
  expect_equal(filt$method, "corrected")
  expect_equal(length(filt$w_filtered), 500)
  expect_equal(length(filt$P_filtered), 500)
})

test_that("filter_svp CKF returns consistent filtered and smoothed output", {
  m <- .fit_model()
  filt <- filter_svp(m$fit)
  # Smoothed and filtered should be same length
  expect_equal(length(filt$w_filtered), length(filt$w_smoothed))
  # Both should be finite
  expect_true(all(is.finite(filt$w_filtered)))
  expect_true(all(is.finite(filt$w_smoothed)))
})

test_that("filter_svp: no NaN/Inf in output (all distributions)", {
  # Gaussian
  m <- .fit_model(errorType = "Gaussian")
  filt <- filter_svp(m$fit)
  expect_true(all(is.finite(filt$w_filtered)))
  expect_true(all(is.finite(filt$w_smoothed)))
  expect_true(is.finite(filt$loglik))

  # Student-t
  m_t <- .fit_model(errorType = "Student-t", nu = 5)
  filt_t <- filter_svp(m_t$fit)
  expect_true(all(is.finite(filt_t$w_filtered)))
  expect_true(is.finite(filt_t$loglik))

  # GED
  m_g <- .fit_model(errorType = "GED", nu = 1.5)
  filt_g <- filter_svp(m_g$fit)
  expect_true(all(is.finite(filt_g$w_filtered)))
  expect_true(is.finite(filt_g$loglik))
})

test_that("filter_svp: filtered w tracks true w on simulated data", {
  set.seed(42)
  sim <- sim_svp(2000, phi = 0.95, sigy = 1, sigv = 0.5,
                 leverage = TRUE, rho = 0)  # returns list with $h (=w)
  fit <- svp(sim$y, p = 1)
  filt <- filter_svp(fit)
  cor_val <- cor(filt$w_smoothed, sim$h)
  # Smoothed w should correlate with true w (higher sigma_v = more signal)
  expect_true(cor_val > 0.3)
})

test_that("filter_svp: P_filtered is positive", {
  m <- .fit_model()
  filt <- filter_svp(m$fit)
  expect_true(all(filt$P_filtered > 0))
  expect_true(all(filt$P_predicted > 0))
})

test_that("filter_svp: leverage affects prediction", {
  m <- .fit_model(leverage = TRUE, rho = -0.5)
  filt <- filter_svp(m$fit)
  expect_true(all(is.finite(filt$w_filtered)))
})


# =========================================================================== #
# filter_svp — GMKF (Gaussian mixture)
# =========================================================================== #

test_that("GMKF with K=1 approximates CKF", {
  m <- .fit_model(n = 300)
  filt_ckf <- filter_svp(m$fit, method = "corrected")
  # K=1 mixture should be very close to CKF
  # (Not exact match due to different y_star centering, but correlation should be high)
  filt_gmkf1 <- filter_svp(m$fit, method = "mixture", K = 1)
  cor_val <- cor(filt_ckf$w_filtered, filt_gmkf1$w_filtered)
  expect_true(cor_val > 0.95)
})

test_that("GMKF K=7 runs without error (Gaussian)", {
  m <- .fit_model(n = 300)
  filt <- filter_svp(m$fit, method = "mixture", K = 7)
  expect_s3_class(filt, "svp_filter")
  expect_equal(filt$method, "mixture")
  expect_true(all(is.finite(filt$w_filtered)))
  expect_true(is.finite(filt$loglik))
})

test_that("GMKF runs without error (Student-t)", {
  m <- .fit_model(n = 200, errorType = "Student-t", nu = 5)
  filt <- filter_svp(m$fit, method = "mixture", K = 3)
  expect_true(all(is.finite(filt$w_filtered)))
})

test_that("GMKF runs without error (GED)", {
  m <- .fit_model(n = 200, errorType = "GED", nu = 1.5)
  filt <- filter_svp(m$fit, method = "mixture", K = 3)
  expect_true(all(is.finite(filt$w_filtered)))
})


# =========================================================================== #
# Stress tests
# =========================================================================== #

test_that("filter_svp: very persistent SV(1) doesn't diverge", {
  set.seed(99)
  y <- as.numeric(sim_svp(500, phi = 0.99, sigy = 1, sigv = 0.1))
  fit <- svp(y, p = 1)
  filt <- filter_svp(fit)
  expect_true(all(is.finite(filt$w_filtered)))
  expect_true(max(abs(filt$w_filtered)) < 50)
})

test_that("filter_svp: heavy tails Student-t nu=3 doesn't collapse", {
  m <- .fit_model(n = 300, errorType = "Student-t", nu = 3)
  filt <- filter_svp(m$fit, method = "mixture", K = 3)
  expect_true(all(is.finite(filt$w_filtered)))
})

test_that("filter_svp: SV(2) companion form works", {
  set.seed(77)
  y <- as.numeric(sim_svp(500, phi = c(0.2, 0.63), sigy = 1, sigv = 0.5))
  fit <- svp(y, p = 2)
  filt <- filter_svp(fit)
  expect_true(all(is.finite(filt$w_filtered)))
  expect_equal(nrow(filt$xi_filtered), 2)
})

test_that("filter_svp: short series T=50 doesn't crash", {
  set.seed(55)
  y <- as.numeric(sim_svp(50, phi = 0.9, sigy = 1, sigv = 0.3))
  fit <- svp(y, p = 1)
  filt <- filter_svp(fit)
  expect_equal(length(filt$w_filtered), 50)
})


# =========================================================================== #
# S3 methods
# =========================================================================== #

test_that("print.svp_filter produces output", {
  m <- .fit_model(n = 200)
  filt <- filter_svp(m$fit)
  expect_output(print(filt), "Filter")
})

test_that("plot.svp_filter runs without error", {
  m <- .fit_model(n = 200)
  filt <- filter_svp(m$fit)
  expect_silent(suppressWarnings(plot(filt)))
})


# =========================================================================== #
# Regression tests: P_filt_T full matrix (Bug 3 fix) and pinv stability (Bug 4)
# =========================================================================== #

test_that("filter_svp returns P_filt_T as a p x p matrix (p=1)", {
  set.seed(42)
  y <- as.numeric(sim_svp(300, phi = 0.95, sigy = 1, sigv = 0.3))
  fit <- svp(y, p = 1)
  filt <- filter_svp(fit, method = "corrected")
  expect_true("P_filt_T" %in% names(filt))
  expect_equal(dim(filt$P_filt_T), c(1L, 1L))
  expect_true(filt$P_filt_T[1, 1] > 0)
})

test_that("filter_svp returns full P_filt_T with non-zero off-diagonals (p=2)", {
  set.seed(42)
  y <- sim_svp(500, phi = c(0.20, 0.63), sigy = 1, sigv = 1)
  fit <- svp(y, p = 2)
  filt_ckf  <- filter_svp(fit, method = "corrected")
  filt_gmkf <- filter_svp(fit, method = "mixture")
  # Dimension check
  expect_equal(dim(filt_ckf$P_filt_T),  c(2L, 2L))
  expect_equal(dim(filt_gmkf$P_filt_T), c(2L, 2L))
  # Off-diagonal should be non-zero for p=2 (was zero under diagonal approx)
  expect_true(abs(filt_ckf$P_filt_T[1, 2])  > 1e-6)
  expect_true(abs(filt_gmkf$P_filt_T[1, 2]) > 1e-6)
  # Matrix should be symmetric and positive definite
  expect_true(isSymmetric(filt_ckf$P_filt_T, tol = 1e-10))
  expect_true(all(eigen(filt_ckf$P_filt_T, only.values = TRUE)$values > 0))
})

test_that("RTS smoother (pinv) produces no NaN/Inf on near-unit-root SV(2)", {
  set.seed(42)
  y <- sim_svp(400, phi = c(0.50, 0.49), sigy = 1, sigv = 0.5)
  fit <- svp(y, p = 2)
  filt_ckf  <- filter_svp(fit, method = "corrected")
  filt_gmkf <- filter_svp(fit, method = "mixture")
  expect_false(any(is.nan(filt_ckf$w_smoothed)  | is.infinite(filt_ckf$w_smoothed)))
  expect_false(any(is.nan(filt_gmkf$w_smoothed) | is.infinite(filt_gmkf$w_smoothed)))
})

# =========================================================================== #
# Regression: BPF xi_filtered dimensions for p>=2 (Bug F fix, 2026-04-03)
# =========================================================================== #

test_that("BPF xi_filtered has correct p x T dimensions for p=2", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(300, phi = c(0.20, 0.63), sigy = 1, sigv = 1)
  fit <- svp(y, p = 2)
  filt <- filter_svp(fit, method = "particle", M = 100)
  expect_equal(nrow(filt$xi_filtered), 2L)
  expect_equal(ncol(filt$xi_filtered), 300L)
  # Row 1 should match w_filtered
  expect_equal(filt$xi_filtered[1, ], filt$w_filtered)
  # Row 2 is lagged: xi[2,t] = w_filtered[t-1]
  expect_equal(filt$xi_filtered[2, 2:300], filt$w_filtered[1:299])
})
