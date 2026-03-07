# wARMASVp 0.1.0

Initial release.

## Estimation

* `svp()`: Closed-form W-ARMA-SV estimation for SV(p) models of any order.
* Gaussian, Student-t, and GED innovation distributions supported for all p.
* Leverage estimation (Gaussian only) via Pearson or Kendall correlation.
* `svpSE()`: Simulation-based standard errors and confidence intervals.

## Simulation

* `sim_svp()`: Simulate SV(p) processes with Gaussian, Student-t, or GED
  innovations, with optional leverage effects (Gaussian).

## Hypothesis Testing

* Local Monte Carlo (LMC) and Maximized Monte Carlo (MMC) tests:
  - `lmc_ar()` / `mmc_ar()`: AR order selection.
  - `lmc_lev()` / `mmc_lev()`: Leverage effects.
  - `lmc_t()` / `mmc_t()`: Student-t vs. Gaussian.
  - `lmc_ged()` / `mmc_ged()`: GED vs. Gaussian.

## Forecasting

* `kalman_filter()`: Kalman filtering and smoothing for SV(p) state-space form.
* `forecast_svp()`: h-step-ahead volatility forecasts with confidence bands.

## Bug Fixes

* Corrected Student-t mean-of-log-squared formula:
  `mu_bar(nu) = psi(1/2) - psi(nu/2) + log(nu - 2)`.
  The previous formula omitted the standardization correction for unit-variance
  Student-t innovations.
