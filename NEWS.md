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

## Convention Changes

* Switched Student-t innovations from standardized (unit variance) to
  unstandardized (raw t(nu) with Var = nu/(nu-2)), matching the SV-t
  literature (Chib, Nardari & Shephard 2002; Jacquier, Polson & Rossi 2004)
  and the SVHT reference paper (Ahsan, Dufour & Rodriguez Rondon 2025b).
  The mean-of-log-squared formula is now:
  `mu_bar(nu) = psi(1/2) - psi(nu/2) + log(nu)`.
  Simulation no longer divides raw Student-t samples by sqrt(nu/(nu-2)).
  GED innovations remain standardized (unit variance), following Nelson (1991).
