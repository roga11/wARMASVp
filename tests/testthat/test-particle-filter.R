# =========================================================================== #
# Tests for Bootstrap Particle Filter (BPF)
# =========================================================================== #

test_that("PF returns correct class and components (Gaussian)", {
  set.seed(42)
  y <- as.numeric(sim_svp(300, phi = 0.95, sigy = 1, sigv = 0.3))
  fit <- svp(y, p = 1)
  pf <- filter_svp(fit, method = "particle", M = 200, seed = 42)
  expect_s3_class(pf, "svp_filter")
  expect_equal(pf$method, "particle")
  expect_equal(length(pf$w_filtered), 300)
  expect_equal(length(pf$P_filtered), 300)
  expect_equal(length(pf$ESS), 300)
  expect_true(is.finite(pf$loglik))
})

test_that("PF: ESS doesn't collapse (Gaussian)", {
  set.seed(42)
  y <- as.numeric(sim_svp(300, phi = 0.95, sigy = 1, sigv = 0.3))
  fit <- svp(y, p = 1)
  pf <- filter_svp(fit, method = "particle", M = 500, seed = 42)
  # Mean ESS should be > M/10 = 50
  expect_true(mean(pf$ESS) > 500 / 10)
})

test_that("PF: filtered states close to CKF (Gaussian)", {
  set.seed(42)
  y <- as.numeric(sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.3))
  fit <- svp(y, p = 1)
  ckf <- filter_svp(fit, method = "corrected")
  pf <- filter_svp(fit, method = "particle", M = 1000, seed = 42)
  cor_val <- cor(pf$w_filtered, ckf$w_filtered)
  expect_true(cor_val > 0.7)
})

test_that("PF: no NaN/Inf in output (all distributions)", {
  # Gaussian
  set.seed(42)
  y <- as.numeric(sim_svp(200, phi = 0.95, sigy = 1, sigv = 0.3))
  fit <- svp(y, p = 1)
  pf <- filter_svp(fit, method = "particle", M = 200, seed = 42)
  expect_true(all(is.finite(pf$w_filtered)))
  expect_true(all(is.finite(pf$ESS)))

  # Student-t
  set.seed(42)
  y_t <- as.numeric(sim_svp(200, phi = 0.95, sigy = 1, sigv = 0.3,
                             errorType = "Student-t", nu = 5))
  fit_t <- svp(y_t, p = 1, errorType = "Student-t")
  pf_t <- filter_svp(fit_t, method = "particle", M = 200, seed = 42)
  expect_true(all(is.finite(pf_t$w_filtered)))

  # GED
  set.seed(42)
  y_g <- as.numeric(sim_svp(200, phi = 0.95, sigy = 1, sigv = 0.3,
                             errorType = "GED", nu = 1.5))
  fit_g <- suppressWarnings(svp(y_g, p = 1, errorType = "GED"))
  pf_g <- filter_svp(fit_g, method = "particle", M = 200, seed = 42)
  expect_true(all(is.finite(pf_g$w_filtered)))
})

test_that("PF: ESS doesn't collapse (Student-t)", {
  set.seed(42)
  y <- as.numeric(sim_svp(300, phi = 0.95, sigy = 1, sigv = 0.3,
                           errorType = "Student-t", nu = 5))
  fit <- svp(y, p = 1, errorType = "Student-t")
  pf <- filter_svp(fit, method = "particle", M = 500, seed = 42)
  expect_true(mean(pf$ESS) > 500 / 10)
})

test_that("PF: leverage (Gaussian)", {
  set.seed(42)
  sim <- sim_svp(300, phi = 0.95, sigy = 1, sigv = 0.3,
                 leverage = TRUE, rho = -0.5)
  fit <- svp(sim$y, p = 1, leverage = TRUE)
  pf <- filter_svp(fit, method = "particle", M = 300, seed = 42)
  expect_true(all(is.finite(pf$w_filtered)))
  expect_true(mean(pf$ESS) > 300 / 10)
})

test_that("PF: leverage (Student-t, exact z recovery via lambda sampling)", {
  set.seed(42)
  sim <- sim_svp(300, phi = 0.95, sigy = 1, sigv = 0.3,
                 errorType = "Student-t", nu = 5,
                 leverage = TRUE, rho = -0.5)
  fit <- svp(sim$y, p = 1, leverage = TRUE, errorType = "Student-t")
  pf <- filter_svp(fit, method = "particle", M = 300, seed = 42)
  expect_true(all(is.finite(pf$w_filtered)))
  expect_true(mean(pf$ESS) > 300 / 10)
})

test_that("PF: leverage (GED, exact z via copula inversion)", {
  set.seed(42)
  sim <- sim_svp(300, phi = 0.95, sigy = 1, sigv = 0.3,
                 errorType = "GED", nu = 1.5,
                 leverage = TRUE, rho = -0.5)
  fit <- svp(sim$y, p = 1, leverage = TRUE, errorType = "GED")
  pf <- filter_svp(fit, method = "particle", M = 300, seed = 42)
  expect_true(all(is.finite(pf$w_filtered)))
  expect_true(mean(pf$ESS) > 300 / 10)
})

test_that("PF: SV(2) companion form works", {
  set.seed(77)
  y <- as.numeric(sim_svp(300, phi = c(0.2, 0.63), sigy = 1, sigv = 0.5))
  fit <- svp(y, p = 2)
  pf <- filter_svp(fit, method = "particle", M = 300, seed = 42)
  expect_true(all(is.finite(pf$w_filtered)))
  expect_true(mean(pf$ESS) > 300 / 10)
})

test_that("PF: reproducibility with same seed", {
  set.seed(42)
  y <- as.numeric(sim_svp(200, phi = 0.95, sigy = 1, sigv = 0.3))
  fit <- svp(y, p = 1)
  pf1 <- filter_svp(fit, method = "particle", M = 200, seed = 123)
  pf2 <- filter_svp(fit, method = "particle", M = 200, seed = 123)
  expect_equal(pf1$w_filtered, pf2$w_filtered)
  expect_equal(pf1$loglik, pf2$loglik)
})

test_that("print.svp_filter works for particle filter", {
  set.seed(42)
  y <- as.numeric(sim_svp(200, phi = 0.95, sigy = 1, sigv = 0.3))
  fit <- svp(y, p = 1)
  pf <- filter_svp(fit, method = "particle", M = 200, seed = 42)
  expect_output(print(pf), "particle")
})
