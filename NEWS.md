# wARMASVp 0.2.0

## New features

* `svp_IC()` and `svp_AR_order()`: AR-order selection for SV(p) models via
  information criteria. Four criteria are returned by default (`BIC_Kalman`,
  `AIC_Kalman`, `BIC_HR`, `AIC_HR`), spanning state-space QML and
  Hannan-Rissanen estimation families; four more (`AICc_Kalman`, `BIC_Whittle`,
  `BIC_YW`, `AIC_YW`) are available opt-in via the `criteria` argument.
  `svp_AR_order()` sweeps over `p = 1, ..., pmax`; both functions read
  `errorType` and `leverage` from the fitted model.
* `lmc_ar()` / `mmc_ar()` now accept `errorType = "Gaussian"`, `"Student-t"`,
  or `"GED"`. The tail parameter is held fixed at the null MLE during
  simulation; innovations are pre-drawn from the corresponding distribution.

## Breaking changes

* `sim_svp()` now always returns a named list `list(y, h, z, v)` of length-n
  vectors (observed returns, log-volatility path, return innovation, volatility
  innovation). The `K` (multiple-replicate) argument has been removed; wrap the
  call in a loop for replicates. Callers that previously relied on `sim_svp()`
  returning a bare vector must now extract `$y`.

## Defaults

* `filter_svp()` and `forecast_svp()` gain a `proxy` argument and now default
  to `proxy = "bayes_optimal"` (was the paper-faithful `"u"`-proxy). For
  Student-t leverage this uses the posterior mean `E[zeta | u]` rather than the
  raw `u`-proxy, which has marginal variance `nu/(nu-2) > 1`. No effect for
  Gaussian, GED, or non-leverage models.

## Bug fixes

* GMKF: corrected the Student-t leverage parameterization in the Gaussian
  mixture Kalman filter.
* Filtering / forecasting: corrected the state-innovation variance `Q` under
  leverage. The filter uses the conditional `Q = sigma_v^2 (1 - delta^2)`; the
  forecaster uses the conditional `Q` at horizon 1 and the marginal
  `Q = sigma_v^2` at horizons >= 2.
* Bootstrap particle filter: Student-t leverage recovery now samples the
  mixing variable from its posterior rather than its prior.
* GED leverage: the CKF/GMKF leverage shift now applies the copula proxy
  rather than using the raw innovation.
* MMC: the observed test statistic S0 is kept fixed during optimization, per
  Dufour (2006, eq. 4.22). Previously recomputed at each optimizer evaluation
  in the leverage, Student-t, and GED tests.
* MMC: default `eps[sigma_y] = 0` in all MMC functions (was 0.3). The
  simulated null distribution is sigma_y-invariant, so varying it is
  unnecessary.

## Performance

* The KSC mixture EM step used by the GMKF (`fit_ksc_mixture()`) is now
  implemented in C++, giving roughly a 12x speedup for Student-t and GED
  filtering.

## Documentation

* `DESCRIPTION`: added the DOI for the JTSA 2025 reference per CRAN reviewer
  feedback.
* Updated the introductory vignette with an AR-order-selection section.

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
  and the SVHT reference paper (Ahsan, Dufour & Rodriguez-Rondon 2026).
  The mean-of-log-squared formula is now:
  `mu_bar(nu) = psi(1/2) - psi(nu/2) + log(nu)`.
  Simulation no longer divides raw Student-t samples by sqrt(nu/(nu-2)).
  GED innovations remain standardized (unit variance), following Nelson (1991).
