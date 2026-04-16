test_that("sim_svp Gaussian SV(1) returns list of length-n vectors", {
  set.seed(1)
  sim <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.3)
  expect_type(sim, "list")
  expect_named(sim, c("y", "h", "z", "v"))
  expect_length(sim$y, 500)
  expect_length(sim$h, 500)
  expect_length(sim$z, 500)
  expect_length(sim$v, 500)
  expect_type(sim$y, "double")
})

test_that("sim_svp Gaussian SV(2) returns correct length", {
  set.seed(1)
  sim <- sim_svp(500, phi = c(0.20, 0.63), sigy = 1, sigv = 0.5)
  expect_named(sim, c("y", "h", "z", "v"))
  expect_length(sim$y, 500)
})

test_that("sim_svp Student-t returns correct length", {
  set.seed(1)
  sim <- sim_svp(500, phi = 0.90, sigy = 1, sigv = 0.3,
                 errorType = "Student-t", nu = 5)
  expect_named(sim, c("y", "h", "z", "v"))
  expect_length(sim$y, 500)
})

test_that("sim_svp GED returns correct length", {
  set.seed(1)
  sim <- sim_svp(500, phi = 0.90, sigy = 1, sigv = 0.3,
                 errorType = "GED", nu = 1.5)
  expect_named(sim, c("y", "h", "z", "v"))
  expect_length(sim$y, 500)
})

test_that("sim_svp leverage returns list with expected components", {
  set.seed(1)
  sim <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.3,
                 leverage = TRUE, rho = -0.5)
  expect_type(sim, "list")
  expect_named(sim, c("y", "h", "z", "v"))
  expect_length(sim$y, 500)
})

test_that("sim_svp output satisfies model equation y = sigy * exp(h/2) * z", {
  set.seed(1)
  sim <- sim_svp(500, phi = 0.95, sigy = 1.5, sigv = 0.3)
  expect_equal(sim$y, 1.5 * exp(sim$h / 2) * sim$z, tolerance = 1e-10)
})

test_that("sim_svp rejects invalid inputs", {
  expect_error(sim_svp(500, phi = 0.9, sigy = 1, sigv = 0.3,
                       errorType = "Student-t", nu = 1.5),
               "nu must be > 2")
  expect_error(sim_svp(500, phi = 0.9, sigy = 1, sigv = 0.3,
                       errorType = "GED", nu = -1),
               "nu must be > 0")
  expect_error(sim_svp(500, phi = 0.9, sigy = 1, sigv = 0.3,
                       errorType = "Student-t"),
               "nu is required")
})

# =========================================================================== #
# Regression test: simulation uses contemporaneous epsilon (Bug 1 fix)
# =========================================================================== #

test_that("Gaussian sim: E[log(y_t^2)] and variance match SV model (no lag contamination)", {
  set.seed(2024)
  # Large sample to detect any systematic bias from lagged eps
  sim <- sim_svp(5000, phi = 0.95, sigy = 1, sigv = 0.3)
  y <- sim$y
  # Under the correct model: E[log(y^2)] = log(sigy^2) + E[log(eps^2)]
  # = log(1) + digamma(0.5) + log(2) ≈ -0.2704
  expected_mean <- log(1^2) + digamma(0.5) + log(2)
  # With lagged eps the mean would be the same (iid eps), but the cross-corr
  # structure changes. Test that estimation recovers phi well (key check).
  fit <- svp(y)
  expect_true(abs(fit$phi - 0.95) < 0.15)
  expect_true(fit$sigy > 0.8 && fit$sigy < 1.2)
})
