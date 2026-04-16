test_that("svp Gaussian SV(1) returns correct class and reasonable estimates", {
  set.seed(42)
  y <- sim_svp(5000, phi = 0.95, sigy = 1, sigv = 0.3)$y
  fit <- svp(y, p = 1, J = 10)
  expect_s3_class(fit, "svp")
  expect_length(fit$phi, 1)
  # phi should be within 0.1 of true value

  expect_true(abs(fit$phi - 0.95) < 0.1)
  # sigy should be positive
  expect_true(fit$sigy > 0)
})

test_that("svp Gaussian SV(2) returns 2 phi values", {
  set.seed(42)
  y <- sim_svp(3000, phi = c(0.20, 0.63), sigy = 1, sigv = 0.5)$y
  fit <- svp(y, p = 2, J = 10)
  expect_s3_class(fit, "svp")
  expect_length(fit$phi, 2)
})

test_that("svp leverage returns rho estimate", {
  set.seed(42)
  sim <- sim_svp(3000, phi = 0.95, sigy = 1, sigv = 0.3,
                 leverage = TRUE, rho = -0.5)
  fit <- svp(sim$y, p = 1, leverage = TRUE)
  expect_s3_class(fit, "svp")
  expect_false(is.na(fit$rho))
  # rho should be negative (true = -0.5)
  expect_true(fit$rho < 0)
})

test_that("svp Student-t returns correct class and reasonable nu", {
  set.seed(42)
  y <- sim_svp(5000, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "Student-t", nu = 5)$y
  fit <- suppressWarnings(svp(y, p = 1, errorType = "Student-t"))
  expect_s3_class(fit, "svp_t")
  expect_true(!is.na(fit$v))
  # nu should be in reasonable range (true = 5)
  expect_true(fit$v > 2 && fit$v < 50)
  # sigy should be positive and not too far from 1
  expect_true(fit$sigy > 0)
  expect_true(abs(fit$sigy - 1) < 0.5)
})

test_that("svp GED returns correct class and reasonable nu", {
  set.seed(42)
  y <- sim_svp(5000, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "GED", nu = 1.5)$y
  fit <- suppressWarnings(svp(y, p = 1, errorType = "GED"))
  expect_s3_class(fit, "svp_ged")
  expect_true(!is.na(fit$v))
  # nu should be in reasonable range (true = 1.5)
  expect_true(fit$v > 0.3 && fit$v < 5)
})

test_that("svp Student-t SV(2) works", {
  set.seed(42)
  y <- sim_svp(3000, phi = c(0.20, 0.63), sigy = 1, sigv = 0.5,
               errorType = "Student-t", nu = 10)$y
  fit <- suppressWarnings(svp(y, p = 2, errorType = "Student-t"))
  expect_s3_class(fit, "svp_t")
  expect_length(fit$phi, 2)
})

test_that("svpSE returns CI and SE components", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.3)$y
  fit <- svp(y, p = 1)
  se <- svpSE(fit, n_sim = 19)
  expect_type(se, "list")
  expect_true("CI" %in% names(se))
  expect_true("SEsim" %in% names(se))
  # CI should have 2 rows (lower, upper) and 3 columns (phi, sigy, sigv)
  expect_equal(nrow(se$CI), 2)
  expect_equal(ncol(se$CI), 3)
})

# =========================================================================== #
# Regression test: nonstationary_ind propagated to svp_t and svp_ged (Bug 6 fix)
# =========================================================================== #

test_that("svp_t object contains nonstationary_ind field", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "Student-t", nu = 5)$y
  fit <- suppressWarnings(svp(y, errorType = "Student-t"))
  expect_true("nonstationary_ind" %in% names(fit))
  expect_type(fit$nonstationary_ind, "logical")
})

test_that("svp_ged object contains nonstationary_ind field", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "GED", nu = 1.5)$y
  fit <- suppressWarnings(svp(y, errorType = "GED"))
  expect_true("nonstationary_ind" %in% names(fit))
  expect_type(fit$nonstationary_ind, "logical")
})

test_that("coef method works for all model types", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.3)$y
  fit <- svp(y)
  expect_true(is.numeric(coef(fit)))

  yt <- sim_svp(1000, phi = 0.90, sigy = 1, sigv = 0.3,
                errorType = "Student-t", nu = 5)$y
  fit_t <- suppressWarnings(svp(yt, errorType = "Student-t"))
  if (!is.na(fit_t$v)) {
    expect_true(is.numeric(coef(fit_t)))
  }
})
