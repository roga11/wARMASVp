test_that("lmc_ar runs and returns svp_test object", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.3)$y
  test <- lmc_ar(y, p_null = 1, p_alt = 2, N = 19)
  expect_s3_class(test, "svp_test")
  expect_true("pval" %in% names(test))
  expect_true(test$pval >= 0 && test$pval <= 1)
})

test_that("mmc_ar runs and returns valid p-value", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.3)$y
  test <- mmc_ar(y, p_null = 1, p_alt = 2, N = 9,
                 method = "pso", maxit = 5)
  expect_true(test$value >= 0 && test$value <= 1)
})

test_that("lmc_lev runs and returns svp_test object", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.3)$y
  test <- lmc_lev(y, p = 1, N = 19)
  expect_s3_class(test, "svp_test")
  expect_true(test$pval >= 0 && test$pval <= 1)
})

test_that("lmc_t runs and returns svp_test object", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "Student-t", nu = 5)$y
  test <- suppressWarnings(lmc_t(y, nu_null = 5, N = 19))
  expect_s3_class(test, "svp_test")
  expect_true(test$pval >= 0 && test$pval <= 1)
})

test_that("lmc_ged runs and returns svp_test object", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "GED", nu = 1.5)$y
  test <- suppressWarnings(lmc_ged(y, nu_null = 1.5, N = 19))
  expect_s3_class(test, "svp_test")
  expect_true(test$pval >= 0 && test$pval <= 1)
})

test_that("print.svp_test works", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.3)$y
  test <- lmc_ar(y, p_null = 1, p_alt = 2, N = 19)
  expect_output(print(test))
})


# =========================================================================== #
# Directional testing
# =========================================================================== #

test_that(".signed_root computes correctly", {
  expect_equal(wARMASVp:::.signed_root(4, 5, 30), sign(5 - 30) * sqrt(4))  # -2
  expect_equal(wARMASVp:::.signed_root(4, 35, 30), sign(35 - 30) * sqrt(4))  # 2
  expect_equal(wARMASVp:::.signed_root(0, 5, 30), 0)
})

test_that(".pvalue_directional computes correctly", {
  S_sim <- c(-3, -1, 0, 1, 2)
  # H1: theta < theta_0 → count S_sim <= S_obs
  pval_less <- wARMASVp:::.pvalue_directional(-2, S_sim, "less")
  expect_equal(pval_less, (1 + sum(S_sim <= -2)) / (5 + 1))  # (1+1)/6 = 1/3

  # H1: theta > theta_0 → count S_sim >= S_obs
  pval_greater <- wARMASVp:::.pvalue_directional(2, S_sim, "greater")
  expect_equal(pval_greater, (1 + sum(S_sim >= 2)) / (5 + 1))  # (1+1)/6 = 1/3
})

test_that("lmc_t with direction='two-sided' gives same structure as before", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "Student-t", nu = 5)$y
  test <- suppressWarnings(lmc_t(y, nu_null = 5, N = 19, direction = "two-sided"))
  expect_s3_class(test, "svp_test")
  expect_equal(test$direction, "two-sided")
  expect_null(test$S_T)
  expect_true(test$pval >= 0 && test$pval <= 1)
})

test_that("lmc_t with direction='less': p-value in [0,1] and S_T present", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "Student-t", nu = 5)$y
  test <- suppressWarnings(lmc_t(y, nu_null = 5, N = 19, direction = "less"))
  expect_s3_class(test, "svp_test")
  expect_equal(test$direction, "less")
  expect_true(!is.null(test$S_T))
  expect_true(test$pval >= 0 && test$pval <= 1)
})

test_that("lmc_ged with direction='less': p-value in [0,1]", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "GED", nu = 1.5)$y
  test <- suppressWarnings(lmc_ged(y, nu_null = 1.5, N = 19, direction = "less"))
  expect_s3_class(test, "svp_test")
  expect_equal(test$direction, "less")
  expect_true(test$pval >= 0 && test$pval <= 1)
})

test_that("print.svp_test shows direction for directional tests", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "Student-t", nu = 5)$y
  test <- suppressWarnings(lmc_t(y, nu_null = 5, N = 19, direction = "less"))
  expect_output(print(test), "less")
  expect_output(print(test), "Signed root")
})


# =========================================================================== #
# SV(p>1) tests — tail distribution and leverage
# =========================================================================== #

test_that("lmc_t with p=2 runs and returns svp_test object", {
  set.seed(42)
  y <- sim_svp(1000, phi = c(0.20, 0.63), sigy = 1, sigv = 1,
               errorType = "Student-t", nu = 3)$y
  test <- suppressWarnings(lmc_t(y, p = 2, nu_null = 3, N = 19))
  expect_s3_class(test, "svp_test")
  expect_true(test$pval >= 0 && test$pval <= 1)
})

test_that("mmc_t with p=2 runs and returns valid p-value", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = c(0.20, 0.63), sigy = 1, sigv = 1,
               errorType = "Student-t", nu = 3)$y
  test <- suppressWarnings(mmc_t(y, p = 2, nu_null = 3, N = 9,
                                  method = "pso", maxit = 5))
  expect_true(test$value >= 0 && test$value <= 1)
})

test_that("lmc_ged with p=2 runs and returns svp_test object", {
  set.seed(42)
  y <- sim_svp(1000, phi = c(0.20, 0.63), sigy = 1, sigv = 1,
               errorType = "GED", nu = 1.5)$y
  test <- suppressWarnings(lmc_ged(y, p = 2, nu_null = 1.5, N = 19))
  expect_s3_class(test, "svp_test")
  expect_true(test$pval >= 0 && test$pval <= 1)
})

test_that("mmc_ged with p=2 runs and returns valid p-value", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = c(0.20, 0.63), sigy = 1, sigv = 1,
               errorType = "GED", nu = 1.5)$y
  test <- suppressWarnings(mmc_ged(y, p = 2, nu_null = 1.5, N = 9,
                                    method = "pso", maxit = 5))
  expect_true(test$value >= 0 && test$value <= 1)
})

test_that("lmc_lev with p=2 runs and returns svp_test object", {
  set.seed(42)
  y <- sim_svp(1000, phi = c(0.20, 0.63), sigy = 1, sigv = 1)$y
  test <- suppressWarnings(lmc_lev(y, p = 2, N = 19))
  expect_s3_class(test, "svp_test")
  expect_true(test$pval >= 0 && test$pval <= 1)
})

test_that("mmc_lev with p=2 runs and returns svp_test object", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = c(0.20, 0.63), sigy = 1, sigv = 1)$y
  test <- suppressWarnings(mmc_lev(y, p = 2, N = 9,
                                    method = "pso", maxit = 5))
  expect_true(test$value >= 0 && test$value <= 1)
})

# =========================================================================== #
# Bartlett/WAmat forwarding tests (Bugs A-D fix, 2026-04-03)
# =========================================================================== #

test_that("mmc_t with Amat='Weighted' (WAmat) runs without error", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "Student-t", nu = 5)$y
  test <- suppressWarnings(mmc_t(y, nu_null = 5, N = 9,
                                  Amat = "Weighted",
                                  method = "pso", maxit = 5))
  expect_true(test$value >= 0 && test$value <= 1)
})

test_that("mmc_ged with Amat='Weighted' (WAmat) runs without error", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "GED", nu = 1.5)$y
  test <- suppressWarnings(mmc_ged(y, nu_null = 1.5, N = 9,
                                    Amat = "Weighted",
                                    method = "pso", maxit = 5))
  expect_true(test$value >= 0 && test$value <= 1)
})

test_that("mmc_lev Gaussian Amat='Weighted' runs without error", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.3)$y
  test <- suppressWarnings(mmc_lev(y, N = 9, Amat = "Weighted",
                                    method = "pso", maxit = 5))
  expect_true(test$value >= 0 && test$value <= 1)
})

test_that("mmc_lev Student-t Amat='Weighted' runs without error", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "Student-t", nu = 5)$y
  test <- suppressWarnings(mmc_lev(y, N = 9, Amat = "Weighted",
                                    errorType = "Student-t",
                                    method = "pso", maxit = 5))
  expect_true(test$value >= 0 && test$value <= 1)
})

test_that("mmc_lev GED Amat='Weighted' runs without error", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "GED", nu = 1.5)$y
  test <- suppressWarnings(mmc_lev(y, N = 9, Amat = "Weighted",
                                    errorType = "GED",
                                    method = "pso", maxit = 5))
  expect_true(test$value >= 0 && test$value <= 1)
})

test_that("mmc_t eps length validation works", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = c(0.20, 0.63), sigy = 1, sigv = 1,
               errorType = "Student-t", nu = 3)$y
  expect_error(
    mmc_t(y, p = 2, nu_null = 3, N = 9, eps = c(0.3, 0.3, 0.3)),
    "eps must have length 4"
  )
})

test_that("mmc_ged eps length validation works", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = c(0.20, 0.63), sigy = 1, sigv = 1,
               errorType = "GED", nu = 1.5)$y
  expect_error(
    mmc_ged(y, p = 2, nu_null = 1.5, N = 9, eps = c(0.3, 0.3, 0.3)),
    "eps must have length 4"
  )
})

# ---- mmc_lev heavy-tail eps tests ----

test_that("mmc_lev Student-t accepts eps of length p+3 (includes nu)", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "Student-t", nu = 5, leverage = TRUE, rho = -0.3)$y
  test <- suppressWarnings(mmc_lev(y, N = 9, errorType = "Student-t",
                                    eps = c(0.3, 0.3, 0.3, 2.0),
                                    method = "pso", maxit = 5))
  expect_true(test$value >= 0 && test$value <= 1)
  expect_length(test$par, 4)  # phi, sigy, sigv, nu
})

test_that("mmc_lev GED accepts eps of length p+3 (includes nu)", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "GED", nu = 1.5, leverage = TRUE, rho = -0.3)$y
  test <- suppressWarnings(mmc_lev(y, N = 9, errorType = "GED",
                                    eps = c(0.3, 0.3, 0.3, 0.5),
                                    method = "pso", maxit = 5))
  expect_true(test$value >= 0 && test$value <= 1)
  expect_length(test$par, 4)  # phi, sigy, sigv, nu
})

test_that("mmc_lev Student-t eps length p+2 still works (backward compat)", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "Student-t", nu = 5, leverage = TRUE, rho = -0.3)$y
  test <- suppressWarnings(mmc_lev(y, N = 9, errorType = "Student-t",
                                    eps = c(0.3, 0.3, 0.3),
                                    method = "pso", maxit = 5))
  expect_true(test$value >= 0 && test$value <= 1)
})

test_that("mmc_lev heavy-tail eps length validation rejects wrong lengths", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "Student-t", nu = 5, leverage = TRUE, rho = -0.3)$y
  # length 1 — wrong
  expect_error(
    mmc_lev(y, N = 9, errorType = "Student-t", eps = 0.3,
            method = "pso", maxit = 5),
    "eps must have length"
  )
  # length p+4 = 5 — wrong
  expect_error(
    mmc_lev(y, N = 9, errorType = "Student-t", eps = rep(0.3, 5),
            method = "pso", maxit = 5),
    "eps must have length"
  )
})

test_that("mmc_lev Student-t p=2 eps of length p+3=5 works", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = c(0.20, 0.63), sigy = 1, sigv = 1.0,
               errorType = "Student-t", nu = 5, leverage = TRUE, rho = -0.3)$y
  test <- suppressWarnings(mmc_lev(y, p = 2, N = 9, errorType = "Student-t",
                                    eps = c(0.3, 0.3, 0.3, 0.3, 2.0),
                                    method = "pso", maxit = 5))
  expect_true(test$value >= 0 && test$value <= 1)
  expect_length(test$par, 5)  # phi1, phi2, sigy, sigv, nu
})


# =========================================================================== #
# Unified Amat interface tests (2026-04-15)
# =========================================================================== #

test_that("lmc_ar with Amat='Weighted' runs (HAC via Amat)", {
  set.seed(42)
  y <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.2)$y
  test <- lmc_ar(y, p_null = 1, p_alt = 2, N = 19, Amat = "Weighted")
  expect_s3_class(test, "svp_test")
  expect_true(test$pval >= 0 && test$pval <= 1)
  expect_true(grepl("Bartlett", test$test_type))
})

test_that("mmc_ar with Amat='Weighted' runs (HAC via Amat)", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.2)$y
  test <- mmc_ar(y, p_null = 1, p_alt = 2, N = 9, Amat = "Weighted",
                 method = "pso", maxit = 5)
  expect_true(test$value >= 0 && test$value <= 1)
})

test_that("lmc_ar rejects user-supplied weighting matrix", {
  set.seed(42)
  y <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.2)$y
  expect_error(
    lmc_ar(y, p_null = 1, p_alt = 2, N = 19, Amat = diag(4)),
    "not supported"
  )
})

test_that("lmc_lev with Amat='Weighted' runs (HAC via Amat)", {
  set.seed(42)
  y <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.2,
               leverage = TRUE, rho = -0.3)$y
  test <- suppressWarnings(lmc_lev(y, p = 1, N = 19, Amat = "Weighted"))
  expect_s3_class(test, "svp_test")
  expect_true(test$pval >= 0 && test$pval <= 1)
})

test_that("mmc_lev with Amat='Weighted' runs (HAC via Amat)", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.2,
               leverage = TRUE, rho = -0.3)$y
  test <- suppressWarnings(mmc_lev(y, p = 1, N = 9, Amat = "Weighted",
                                    method = "pso", maxit = 5))
  expect_true(test$value >= 0 && test$value <= 1)
})

test_that("Amat takes precedence over Bartlett", {
  set.seed(42)
  y <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.2)$y
  # Amat = "Weighted" should override Bartlett = FALSE
  test <- lmc_ar(y, p_null = 1, p_alt = 2, N = 19,
                 Bartlett = FALSE, Amat = "Weighted")
  expect_true(grepl("Bartlett", test$test_type))
})

test_that("lmc_t default is now identity weighting (Bartlett=FALSE)", {
  set.seed(42)
  y <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.2,
               errorType = "Student-t", nu = 5)$y
  # Default call — should use identity (no HAC), same as Bartlett=FALSE
  test_default <- suppressWarnings(lmc_t(y, nu_null = 5, N = 19))
  test_explicit <- suppressWarnings(lmc_t(y, nu_null = 5, N = 19,
                                           Bartlett = FALSE, Amat = NULL))
  expect_equal(test_default$s0, test_explicit$s0)
})

test_that("lmc_t Amat='Weighted' gives HAC weighting", {
  set.seed(42)
  y <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.2,
               errorType = "Student-t", nu = 5)$y
  # Amat="Weighted" should produce same result as Bartlett=TRUE
  test_amat <- suppressWarnings(lmc_t(y, nu_null = 5, N = 19,
                                       Amat = "Weighted"))
  test_bart <- suppressWarnings(lmc_t(y, nu_null = 5, N = 19,
                                       Bartlett = TRUE))
  expect_equal(test_amat$s0, test_bart$s0)
})


# =========================================================================== #
# Heavy-tail AR-order tests (v0.2.0)
# =========================================================================== #

test_that("lmc_ar with errorType='Student-t' runs", {
  set.seed(42)
  y <- sim_svp(400, phi = 0.95, sigy = 1, sigv = 0.5,
               errorType = "Student-t", nu = 5)$y
  test <- suppressWarnings(
    lmc_ar(y, p_null = 1, p_alt = 2, J = 10, N = 19,
           Bartlett = TRUE, errorType = "Student-t")
  )
  expect_s3_class(test, "svp_test")
  expect_equal(test$errorType, "Student-t")
  expect_true(test$pval >= 0 && test$pval <= 1)
  expect_true(grepl("Student-t", test$test_type))
  expect_true(test$s0 >= 0)  # capping enforces non-negative
})

test_that("lmc_ar with errorType='GED' runs", {
  set.seed(42)
  y <- sim_svp(400, phi = 0.95, sigy = 1, sigv = 0.5,
               errorType = "GED", nu = 1.5)$y
  test <- suppressWarnings(
    lmc_ar(y, p_null = 1, p_alt = 2, J = 10, N = 19,
           Bartlett = TRUE, errorType = "GED")
  )
  expect_s3_class(test, "svp_test")
  expect_equal(test$errorType, "GED")
  expect_true(test$pval >= 0 && test$pval <= 1)
  expect_true(grepl("GED", test$test_type))
  expect_true(test$s0 >= 0)
})

test_that("mmc_ar with errorType='Student-t' runs", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(400, phi = 0.95, sigy = 1, sigv = 0.5,
               errorType = "Student-t", nu = 5)$y
  test <- suppressWarnings(
    mmc_ar(y, p_null = 1, p_alt = 2, J = 10, N = 9,
           Bartlett = TRUE, errorType = "Student-t",
           method = "pso", maxit = 5)
  )
  expect_equal(test$errorType, "Student-t")
  expect_true(test$value >= 0 && test$value <= 1)
  expect_true(test$s0 >= 0)
})

test_that("mmc_ar with errorType='GED' runs", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(400, phi = 0.95, sigy = 1, sigv = 0.5,
               errorType = "GED", nu = 1.5)$y
  test <- suppressWarnings(
    mmc_ar(y, p_null = 1, p_alt = 2, J = 10, N = 9,
           Bartlett = TRUE, errorType = "GED",
           method = "pso", maxit = 5)
  )
  expect_equal(test$errorType, "GED")
  expect_true(test$value >= 0 && test$value <= 1)
  expect_true(test$s0 >= 0)
})

test_that("lmc_ar caps negative observed test statistic at 1e-10", {
  set.seed(99)
  y <- sim_svp(300, phi = 0.7, sigy = 1, sigv = 0.4)$y
  test <- suppressWarnings(
    lmc_ar(y, p_null = 1, p_alt = 2, J = 10, N = 19, Bartlett = TRUE)
  )
  expect_true(test$s0 >= 1e-10)  # always non-negative after capping
})

test_that("lmc_ar errorType validation rejects invalid argument", {
  set.seed(42)
  y <- sim_svp(200, phi = 0.9, sigy = 1, sigv = 0.3)$y
  expect_error(
    lmc_ar(y, p_null = 1, p_alt = 2, errorType = "Cauchy")
  )
})
