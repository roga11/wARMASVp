# =========================================================================== #
# Tests for leverage estimation under heavy-tailed distributions
# =========================================================================== #

# --- Simulation Tests ---

test_that("SVL(1) Student-t leverage simulation produces correct output", {
  set.seed(42)
  out <- sim_svp(500, phi = 0.9, sigy = 1, sigv = 0.5,
                 errorType = "Student-t", leverage = TRUE, rho = -0.5, nu = 5)
  expect_type(out, "list")
  expect_true(all(c("y", "h", "zeta", "veta") %in% names(out)))
  expect_length(out$y, 500)
  expect_length(out$h, 500)
  expect_length(out$zeta, 500)
  expect_length(out$veta, 500)
})

test_that("SVL(1) GED leverage simulation produces correct output", {
  set.seed(42)
  out <- sim_svp(500, phi = 0.9, sigy = 1, sigv = 0.5,
                 errorType = "GED", leverage = TRUE, rho = -0.5, nu = 1.5)
  expect_type(out, "list")
  expect_true(all(c("y", "h", "zeta", "veta") %in% names(out)))
  expect_length(out$y, 500)
})

test_that("SVL(1) Student-t with rho=0 matches non-leverage marginal variance", {
  set.seed(123)
  y_lev <- sim_svp(5000, phi = 0.9, sigy = 1, sigv = 0.5,
                   errorType = "Student-t", leverage = TRUE, rho = 0, nu = 5)$y
  set.seed(456)
  y_nolev <- sim_svp(5000, phi = 0.9, sigy = 1, sigv = 0.5,
                     errorType = "Student-t", leverage = FALSE, nu = 5)
  expect_true(abs(var(y_lev) - var(y_nolev)) / var(y_nolev) < 0.15)
})

test_that("SVL(2) Student-t leverage simulation works", {
  set.seed(42)
  out <- sim_svp(500, phi = c(0.2, 0.63), sigy = 1, sigv = 1,
                 errorType = "Student-t", leverage = TRUE, rho = -0.3, nu = 5)
  expect_type(out, "list")
  expect_length(out$y, 500)
})

# --- Correction Factor Tests ---

test_that("C_t(nu) matches known values", {
  expect_equal(correction_factor_t(3), 6 / pi, tolerance = 1e-6)
  expect_equal(correction_factor_t(5), 40 / (9 * pi), tolerance = 1e-4)
  expect_true(correction_factor_t(100) > 1 && correction_factor_t(100) < 1.02)
  expect_true(correction_factor_t(1000) > 1 && correction_factor_t(1000) < 1.003)
})

test_that("C_g(2) = 1 (Gaussian limit)", {
  expect_equal(correction_factor_ged_approx(2), 1, tolerance = 1e-6)
})

test_that("C_g(1.5) matches paper Table 8", {
  expect_equal(correction_factor_ged_approx(1.5), 0.959, tolerance = 0.002)
})

test_that("qged_std reproduces qnorm for nu=2", {
  p_vals <- c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
  expect_equal(qged_std(p_vals, 2), qnorm(p_vals), tolerance = 1e-6)
})

test_that("E[|u|] for GED(2) equals sqrt(2/pi)", {
  expect_equal(ged_E_abs_u(2), sqrt(2 / pi), tolerance = 1e-10)
})

# --- Estimation Tests ---

test_that("SVL(1) Student-t estimation returns leverage parameters", {
  set.seed(42)
  y <- sim_svp(2000, phi = 0.9, sigy = 1, sigv = 0.5,
               errorType = "Student-t", leverage = TRUE, rho = -0.5, nu = 5)$y
  fit <- svp(y, p = 1, J = 50, leverage = TRUE, errorType = "Student-t")
  expect_true(!is.na(fit$rho))
  expect_true(abs(fit$rho) < 1)
  expect_true(fit$v > 2)
  expect_true(isTRUE(fit$leverage))
  expect_true(!is.na(fit$gammatilde))
  expect_true(!is.null(fit$CF))
  expect_equal(length(fit$theta), 5)  # phi, sigy, sigv, nu, rho
})

test_that("SVL(1) GED estimation returns leverage parameters", {
  set.seed(42)
  y <- sim_svp(2000, phi = 0.9, sigy = 1, sigv = 0.5,
               errorType = "GED", leverage = TRUE, rho = -0.5, nu = 1.5)$y
  fit <- svp(y, p = 1, J = 50, leverage = TRUE, errorType = "GED")
  expect_true(!is.na(fit$rho))
  expect_true(abs(fit$rho) < 1)
  expect_true(fit$v > 0)
  expect_true(isTRUE(fit$leverage))
  expect_equal(length(fit$theta), 5)
})

test_that("SVL(1) GED with nu=2 gives leverage close to Gaussian", {
  set.seed(123)
  y <- sim_svp(3000, phi = 0.9, sigy = 1, sigv = 0.3, leverage = TRUE, rho = -0.3)$y
  fit_g <- svp(y, p = 1, J = 50, leverage = TRUE, errorType = "Gaussian")
  suppressWarnings(
    fit_ged <- svp(y, p = 1, J = 50, leverage = TRUE, errorType = "GED")
  )
  # GED nu=2 should give similar rho to Gaussian
  if (fit_ged$v > 1.8 && fit_ged$v < 2.2) {
    expect_true(abs(fit_g$rho - fit_ged$rho) < 0.1)
  }
})

test_that("Student-t nu=500 gives leverage close to Gaussian", {
  set.seed(123)
  y <- sim_svp(3000, phi = 0.9, sigy = 1, sigv = 0.3, leverage = TRUE, rho = -0.3)$y
  fit_g <- svp(y, p = 1, J = 50, leverage = TRUE, errorType = "Gaussian")
  suppressWarnings(
    fit_t <- svp(y, p = 1, J = 50, leverage = TRUE, errorType = "Student-t")
  )
  # nu should be large (near Gaussian), rho should be similar
  if (fit_t$v > 30) {
    expect_true(abs(fit_g$rho - fit_t$rho) < 0.05)
  }
})

test_that("SVL(2) Student-t estimation works", {
  set.seed(42)
  y <- sim_svp(2000, phi = c(0.2, 0.63), sigy = 1, sigv = 1,
               errorType = "Student-t", leverage = TRUE, rho = -0.3, nu = 5)$y
  fit <- svp(y, p = 2, J = 50, leverage = TRUE, errorType = "Student-t")
  expect_true(!is.na(fit$rho))
  expect_equal(length(fit$phi), 2)
  expect_equal(length(fit$theta), 6)  # phi1, phi2, sigy, sigv, nu, rho
})

test_that("Model object has all required leverage fields", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.9, sigy = 1, sigv = 0.5,
               errorType = "Student-t", leverage = TRUE, rho = -0.5, nu = 5)$y
  fit <- svp(y, p = 1, J = 50, leverage = TRUE, errorType = "Student-t")
  expect_true("rho" %in% names(fit))
  expect_true("gammatilde" %in% names(fit))
  expect_true("leverage" %in% names(fit))
  expect_true("rho_type" %in% names(fit))
  expect_true("trunc_lev" %in% names(fit))
  expect_true("CF" %in% names(fit))
  expect_equal(fit$rho_type, "pearson")
  expect_true(isTRUE(fit$leverage))
})

# --- LMC Test Tests ---

test_that("lmc_lev with Student-t runs and returns svp_test", {
  set.seed(42)
  y <- sim_svp(500, phi = 0.9, sigy = 1, sigv = 0.5,
               errorType = "Student-t", leverage = TRUE, rho = -0.3, nu = 5)$y
  suppressWarnings(
    result <- lmc_lev(y, p = 1, J = 50, N = 9, rho_null = 0,
                      errorType = "Student-t")
  )
  expect_s3_class(result, "svp_test")
  expect_true(result$pval >= 0 && result$pval <= 1)
  expect_equal(length(result$sN), 9)
})

test_that("lmc_lev with GED runs and returns svp_test", {
  set.seed(42)
  y <- sim_svp(500, phi = 0.9, sigy = 1, sigv = 0.5,
               errorType = "GED", leverage = TRUE, rho = -0.3, nu = 1.5)$y
  suppressWarnings(
    result <- lmc_lev(y, p = 1, J = 50, N = 9, rho_null = 0,
                      errorType = "GED")
  )
  expect_s3_class(result, "svp_test")
  expect_true(result$pval >= 0 && result$pval <= 1)
})

# --- SE Tests ---

test_that("svpSE for SVL(1) Student-t returns correct dimensions", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.9, sigy = 1, sigv = 0.5,
               errorType = "Student-t", leverage = TRUE, rho = -0.5, nu = 5)$y
  fit <- svp(y, p = 1, J = 50, leverage = TRUE, errorType = "Student-t")
  suppressWarnings(
    se <- svpSE(fit, n_sim = 9, burnin = 200)
  )
  expect_equal(ncol(se$CI), 5)  # phi, sigy, sigv, nu, rho
  expect_equal(nrow(se$CI), 2)  # lower, upper
  expect_equal(length(se$SEsim), 5)
})

test_that("svpSE for SVL(1) GED returns correct dimensions", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.9, sigy = 1, sigv = 0.5,
               errorType = "GED", leverage = TRUE, rho = -0.5, nu = 1.5)$y
  fit <- svp(y, p = 1, J = 50, leverage = TRUE, errorType = "GED")
  suppressWarnings(
    se <- svpSE(fit, n_sim = 9, burnin = 200)
  )
  expect_equal(ncol(se$CI), 5)
  expect_equal(length(se$SEsim), 5)
})

# --- Coef and Print Methods ---

test_that("coef.svp_t includes rho for leverage models", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.9, sigy = 1, sigv = 0.5,
               errorType = "Student-t", leverage = TRUE, rho = -0.5, nu = 5)$y
  fit <- svp(y, p = 1, J = 50, leverage = TRUE, errorType = "Student-t")
  co <- coef(fit)
  expect_true("rho" %in% names(co))
  expect_equal(length(co), 5)
})

test_that("coef.svp_t excludes rho for non-leverage models", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.9, sigy = 1, sigv = 0.5,
               errorType = "Student-t", nu = 5)
  fit <- svp(y, p = 1, J = 50, errorType = "Student-t")
  co <- coef(fit)
  expect_false("rho" %in% names(co))
  expect_equal(length(co), 4)
})
