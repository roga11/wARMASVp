test_that("lmc_ar runs and returns svp_test object", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.3)
  test <- lmc_ar(y, p_null = 1, p_alt = 2, N = 19)
  expect_s3_class(test, "svp_test")
  expect_true("pval" %in% names(test))
  expect_true(test$pval >= 0 && test$pval <= 1)
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
