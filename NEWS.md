# wARMASVp 0.1.0.9000

* MMC: Fixed observed test statistic S0 is now kept fixed during optimization,
  per Dufour (2006, eq 4.22). Previously recomputed at each optimizer
  evaluation in leverage, Student-t, and GED tests.
* MMC: Default `eps[sigma_y] = 0` in all MMC functions (was 0.3). The null
  distribution is sigma_y-invariant, so varying it is unnecessary.
* DESCRIPTION: Added DOI for JTSA 2025 reference per CRAN reviewer feedback.

# wARMASVp 0.1.0

Initial release.

## Estimation

* `svp()`: Closed-form W-ARMA-SV estimation for SV(p) models of any order.
* Gaussian, Student-t, and GED innovation distributions supported for all p.
* Leverage estimation for all distributions: closed-form for Gaussian and
  Student-t, exact root-finding for GED.
* `svpSE()`: Simulation-based standard errors and confidence intervals.

## Simulation

* `sim_svp()`: Simulate SV(p) processes with Gaussian, Student-t, or GED
  innovations, with optional leverage effects for all distributions.

## Hypothesis Testing

* Local Monte Carlo (LMC) and Maximized Monte Carlo (MMC) tests based on
  Dufour (2006), with fixed-innovation MMC for exact finite-sample inference:
  - `lmc_ar()` / `mmc_ar()`: AR order selection.
  - `lmc_lev()` / `mmc_lev()`: Leverage effects (all distributions).
  - `lmc_t()` / `mmc_t()`: Student-t vs. Gaussian (with directional testing).
  - `lmc_ged()` / `mmc_ged()`: GED vs. Gaussian (with directional testing).
* All test procedures support general SV(p) (any order).

## Filtering

* `filter_svp()`: Kalman filtering and smoothing with three methods:
  - Corrected Kalman Filter (CKF): Gaussian approximation, fast.
  - Gaussian Mixture Kalman Filter (GMKF): KSC (1998) 7-component mixture,
    recommended.
  - Bootstrap Particle Filter (BPF): exact density weights, benchmark.

## Forecasting

* `forecast_svp()`: Multi-step ahead volatility forecasts with MSE-based
  confidence bands. Supports log-variance, variance, and volatility output
  scales.

## Convention Changes

* Switched Student-t innovations from standardized (unit variance) to
  unstandardized (raw t(nu) with Var = nu/(nu-2)), matching the SV-t
  literature (Chib, Nardari & Shephard 2002; Jacquier, Polson & Rossi 2004)
  and the SVHT reference paper (Ahsan, Dufour & Rodriguez Rondon 2025b).
  The mean-of-log-squared formula is now:
  `mu_bar(nu) = psi(1/2) - psi(nu/2) + log(nu)`.
  Simulation no longer divides raw Student-t samples by sqrt(nu/(nu-2)).
  GED innovations remain standardized (unit variance), following Nelson (1991).
