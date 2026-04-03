test_that("lmc_ar runs and returns svp_test object", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.3)
  test <- lmc_ar(y, p_null = 1, p_alt = 2, N = 19)
  expect_s3_class(test, "svp_test")
  expect_true("pval" %in% names(test))
  expect_true(test$pval >= 0 && test$pval <= 1)
})

test_that("mmc_ar runs and returns valid p-value", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.3)
  test <- mmc_ar(y, p_null = 1, p_alt = 2, N = 9,
                 method = "pso", maxit = 5)
  expect_true(test$value >= 0 && test$value <= 1)
})

test_that("lmc_lev runs and returns svp_test object", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.3)
  test <- lmc_lev(y, p = 1, N = 19)
  expect_s3_class(test, "svp_test")
  expect_true(test$pval >= 0 && test$pval <= 1)
})

test_that("lmc_t runs and returns svp_test object", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "Student-t", nu = 5)
  test <- suppressWarnings(lmc_t(y, nu_null = 5, N = 19))
  expect_s3_class(test, "svp_test")
  expect_true(test$pval >= 0 && test$pval <= 1)
})

test_that("lmc_ged runs and returns svp_test object", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "GED", nu = 1.5)
  test <- suppressWarnings(lmc_ged(y, nu_null = 1.5, N = 19))
  expect_s3_class(test, "svp_test")
  expect_true(test$pval >= 0 && test$pval <= 1)
})

test_that("print.svp_test works", {
  set.seed(42)
  y <- as.numeric(sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.3))
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
  y <- as.numeric(sim_svp(1000, phi = 0.90, sigy = 1, sigv = 0.3,
                           errorType = "Student-t", nu = 5))
  test <- suppressWarnings(lmc_t(y, nu_null = 5, N = 19, direction = "two-sided"))
  expect_s3_class(test, "svp_test")
  expect_equal(test$direction, "two-sided")
  expect_null(test$S_T)
  expect_true(test$pval >= 0 && test$pval <= 1)
})

test_that("lmc_t with direction='less': p-value in [0,1] and S_T present", {
  set.seed(42)
  y <- as.numeric(sim_svp(1000, phi = 0.90, sigy = 1, sigv = 0.3,
                           errorType = "Student-t", nu = 5))
  test <- suppressWarnings(lmc_t(y, nu_null = 5, N = 19, direction = "less"))
  expect_s3_class(test, "svp_test")
  expect_equal(test$direction, "less")
  expect_true(!is.null(test$S_T))
  expect_true(test$pval >= 0 && test$pval <= 1)
})

test_that("lmc_ged with direction='less': p-value in [0,1]", {
  set.seed(42)
  y <- as.numeric(sim_svp(1000, phi = 0.90, sigy = 1, sigv = 0.3,
                           errorType = "GED", nu = 1.5))
  test <- suppressWarnings(lmc_ged(y, nu_null = 1.5, N = 19, direction = "less"))
  expect_s3_class(test, "svp_test")
  expect_equal(test$direction, "less")
  expect_true(test$pval >= 0 && test$pval <= 1)
})

test_that("print.svp_test shows direction for directional tests", {
  set.seed(42)
  y <- as.numeric(sim_svp(1000, phi = 0.90, sigy = 1, sigv = 0.3,
                           errorType = "Student-t", nu = 5))
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
               errorType = "Student-t", nu = 3)
  test <- suppressWarnings(lmc_t(y, p = 2, nu_null = 3, N = 19))
  expect_s3_class(test, "svp_test")
  expect_true(test$pval >= 0 && test$pval <= 1)
})

test_that("mmc_t with p=2 runs and returns valid p-value", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = c(0.20, 0.63), sigy = 1, sigv = 1,
               errorType = "Student-t", nu = 3)
  test <- suppressWarnings(mmc_t(y, p = 2, nu_null = 3, N = 9,
                                  method = "pso", maxit = 5))
  expect_true(test$value >= 0 && test$value <= 1)
})

test_that("lmc_ged with p=2 runs and returns svp_test object", {
  set.seed(42)
  y <- sim_svp(1000, phi = c(0.20, 0.63), sigy = 1, sigv = 1,
               errorType = "GED", nu = 1.5)
  test <- suppressWarnings(lmc_ged(y, p = 2, nu_null = 1.5, N = 19))
  expect_s3_class(test, "svp_test")
  expect_true(test$pval >= 0 && test$pval <= 1)
})

test_that("mmc_ged with p=2 runs and returns valid p-value", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = c(0.20, 0.63), sigy = 1, sigv = 1,
               errorType = "GED", nu = 1.5)
  test <- suppressWarnings(mmc_ged(y, p = 2, nu_null = 1.5, N = 9,
                                    method = "pso", maxit = 5))
  expect_true(test$value >= 0 && test$value <= 1)
})

test_that("lmc_lev with p=2 runs and returns svp_test object", {
  set.seed(42)
  y <- sim_svp(1000, phi = c(0.20, 0.63), sigy = 1, sigv = 1)
  test <- suppressWarnings(lmc_lev(y, p = 2, N = 19))
  expect_s3_class(test, "svp_test")
  expect_true(test$pval >= 0 && test$pval <= 1)
})

test_that("mmc_lev with p=2 runs and returns svp_test object", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = c(0.20, 0.63), sigy = 1, sigv = 1)
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
               errorType = "Student-t", nu = 5)
  test <- suppressWarnings(mmc_t(y, nu_null = 5, N = 9,
                                  Amat = "Weighted",
                                  method = "pso", maxit = 5))
  expect_true(test$value >= 0 && test$value <= 1)
})

test_that("mmc_ged with Amat='Weighted' (WAmat) runs without error", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "GED", nu = 1.5)
  test <- suppressWarnings(mmc_ged(y, nu_null = 1.5, N = 9,
                                    Amat = "Weighted",
                                    method = "pso", maxit = 5))
  expect_true(test$value >= 0 && test$value <= 1)
})

test_that("mmc_lev Gaussian Bartlett=TRUE runs without error", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.3)
  test <- suppressWarnings(mmc_lev(y, N = 9, Bartlett = TRUE,
                                    method = "pso", maxit = 5))
  expect_true(test$value >= 0 && test$value <= 1)
})

test_that("mmc_lev Student-t Bartlett=TRUE runs without error", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "Student-t", nu = 5)
  test <- suppressWarnings(mmc_lev(y, N = 9, Bartlett = TRUE,
                                    errorType = "Student-t",
                                    method = "pso", maxit = 5))
  expect_true(test$value >= 0 && test$value <= 1)
})

test_that("mmc_lev GED Bartlett=TRUE runs without error", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "GED", nu = 1.5)
  test <- suppressWarnings(mmc_lev(y, N = 9, Bartlett = TRUE,
                                    errorType = "GED",
                                    method = "pso", maxit = 5))
  expect_true(test$value >= 0 && test$value <= 1)
})

test_that("mmc_t eps length validation works", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = c(0.20, 0.63), sigy = 1, sigv = 1,
               errorType = "Student-t", nu = 3)
  expect_error(
    mmc_t(y, p = 2, nu_null = 3, N = 9, eps = c(0.3, 0.3, 0.3)),
    "eps must have length 4"
  )
})

test_that("mmc_ged eps length validation works", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(500, phi = c(0.20, 0.63), sigy = 1, sigv = 1,
               errorType = "GED", nu = 1.5)
  expect_error(
    mmc_ged(y, p = 2, nu_null = 1.5, N = 9, eps = c(0.3, 0.3, 0.3)),
    "eps must have length 4"
  )
})
