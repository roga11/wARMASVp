# wARMASVp (development version)

## New features

* **`svp()` objects now report the measurement-noise scale.** Every fit carries
  `sigeps` (the distribution-free moment estimate of sigma_eps, the standard
  deviation of log z_t^2; `NA` when the variance estimate is not positive),
  `sigeps_nu` (the value the fitted model actually uses: pi/sqrt(2) for
  Gaussian errors, sigma_eps(nu-hat) for Student-t and GED), `mubar`
  (E[log z_t^2] under the fitted model, so that `mu = log(sigy^2) + mubar`) and,
  for Student-t and GED, `nu_boundary`. `print()` shows them.
  `winsorize_eps` is now honoured for Gaussian fits (it affects only `sigeps`).
  No estimate changes.

* **`svp(..., sigeps = "estimated")` for Gaussian errors.** Replaces the
  log-chi-square constant `pi^2/2` by the distribution-free moment estimate of
  `sigma_eps^2`, recomputes `sigv` from it, feeds it to the corrected Kalman
  filter, and identifies `sigy` by `E[log z^2] = -sigma_eps^2/2`. Intended for
  measurement-error proxies such as realized volatility, where `var(log y^2)`
  is far below `pi^2/2` and the theoretical constant makes `sigv` come from
  the `abs()` branch and the filter over-weight the noise. Returns
  `sigeps_method` and `sigeps_fallback`. Not available with leverage,
  `sigvMethod = "hybrid"`, the mixture or particle filters, or `svpSE()`.
  Default `"theoretical"`: all existing results unchanged. `svp_AR_order()`
  gains the same `sigeps` argument.

## Bug fixes

* The Gaussian `sigy` step and the Gaussian moment conditions of the leverage
  and AR-order tests now use the exact `E[log z^2] = digamma(1/2) + log(2)`
  (-1.2703628454614782) instead of the rounded -1.2704 that the Kalman filter
  never used. Gaussian `sigy` moves by a factor of about 1 - 1.9e-5; the filter
  intercept now equals `mu` exactly and `mu == log(sigy^2) + mubar` holds for
  every fit.

* Student-t and GED `sigy` at `p = 1` used `mean(log(y^2))` without the `del`
  offset, unlike every other branch; a series containing an exact zero gave
  `sigy = 0`. All branches now use `mean(log(y^2 + del))`.

* **Precomputed optimal mixture tables for the Gaussian-mixture Kalman filter
  (GMKF), now the default.** `filter_svp()`, `forecast_svp()`, `svp_IC()` and
  `svp_AR_order()` gain a `mixture` argument selecting the source of the GMKF
  measurement mixture:

  - `mixture = "optimal"` (**new default**): tables fitted offline by
    minimising Kullback-Leibler divergence against the exact density of the
    measurement noise, for every error family (Gaussian, Student-t, GED) and
    `K` in {7, 10, 12, 16}, evaluated at the fitted `nu` by spline
    interpolation. Deterministic (no seed anywhere in the path), roughly 30x
    faster per `method = "mixture"` filter call, and measured closer to the
    exact density than the sampled EM fit at every `nu` tested — with
    heavy-tail filtering now exactly as cheap as Gaussian.
  - `mixture = "ksc"`: the Kim, Shephard and Chib (1998) Table 4 mixture
    (Gaussian errors, `K = 7` only), retained as the industry standard for
    comparability with other software. KSC built that table for their MCMC
    sampler, with those assumptions baked in; for the filtering done in this
    package, `"optimal"` is the better default.
  - `mixture = "em"`: the previous behaviour — an EM fit to 10,000 random
    draws — retained for back-compatibility. Given the same `seed`, it
    reproduces the mixture of version 0.2.0.9002 exactly; the filter output
    differs from 0.2.0.9002 only through the `sigy` fixes listed under Bug
    fixes below (of order 1e-5 in the log-likelihood). Supports any
    `K >= 1`. (For Gaussian errors at `K = 7` this path has always returned
    the KSC table, so `"em"` and `"ksc"` coincide there.)

  The shipped tables live in `inst/extdata/` and can be inspected with the
  new exported accessor `gmkf_table()`; the exact scripts that built them
  ship in `inst/gmkf-tables/` with a README describing how to rebuild and
  verify them (each table carries a provenance attribute with the md5 of the
  script and helpers that produced it). The build scripts alone use
  `matrixStats` and `clue` (added to `Suggests`).

* **Behaviour changes under default settings.** Because `"optimal"` is the
  default, results from `filter_svp()`, `forecast_svp()`, `svp_IC()` and
  `svp_AR_order()` with `method`/`filter_method = "mixture"` change
  (deterministically): for Student-t/GED the table sits at the sampled EM's
  seed-median (removing seed noise that could move `BIC_Kalman` differences
  across `p` by several nats at large `T`); for Gaussian the change from the
  KSC table is a genuine fit improvement of about 1e-3 nats per observation.
  To reproduce old results exactly, set `mixture = "em"` (all families) or
  `mixture = "ksc"` (Gaussian, `K = 7`).

* **Compatibility notes.** `filter_svp(method = "mixture")` with `K` outside
  {7, 10, 12, 16} now requires `mixture = "em"` (the default errors, naming
  the supported sizes — no silent fallback; the same validation runs inside
  `svp_IC()`/`svp_AR_order()` so an unsupported combination errors loudly
  instead of returning `NA` criteria). The `seed` argument has no effect on
  the mixture under `"optimal"`/`"ksc"`. Filter output gains a
  `mixture_source` element recording the source used; the existing (and now
  documented) `mixture` element — the mixture actually used — loses its
  `KL_div` field, which was hardcoded `NA` and read by nothing. The new
  `mixture` argument sits after `proxy`, so code calling these functions with
  purely positional arguments past `proxy` must be updated (named arguments
  are unaffected).

* The warning issued when the GED shape parameter hits its upper bound was
  **factually backwards**. It stated that the tails were "indistinguishable
  from Gaussian"; in fact GED(20) is close to *uniform* and has tails roughly
  a million times **thinner** than Gaussian (kurtosis 1.82 versus 3, and
  `P(|z| > 2)` of 4.2e-08 versus 0.0455). Gaussian is GED(2), in the middle of
  the admissible range, not at its top. The message appears to have been copied
  from the Student-t boundary warning, where `nu -> Inf` genuinely *is* Gaussian
  and the wording is correct.

  The replacement says what a boundary estimate actually means: the shape is
  **not identified** there, because `nu` in `[20, Inf)` spans a range of the
  matched moment `(2/nu)^2 * trigamma(1/nu)` of width only ~0.015 — within
  sampling error at typical sample sizes — so `nu_hat = 20` should be read as
  "large and indistinguishable from infinity" rather than as an estimate of 20.
  It also notes that near-uniform, thinner-than-Gaussian innovations are
  implausible for asset returns and usually indicate that the AR order,
  leverage specification, or error family is misspecified, and suggests
  comparing against `errorType = "Student-t"`.

  The Student-t boundary warning is unchanged and remains correct.

* `lmc_ar()` and `mmc_ar()` no longer cap the test statistic at `1e-10`. Both
  the observed statistic and every simulated null draw are now retained
  uncapped, so the simulated null distribution is continuous.

  The `Bartlett = TRUE` (HAC) statistic is an LR-type test built from GMM
  criteria rather than maximized likelihoods, and the null and alternative
  criteria are formed under separately estimated models with their own implied
  long-run variances. Nothing is maximized over a nesting parameter space, so
  the non-negativity of a genuine likelihood ratio does not carry over:
  negative values are structurally possible and informative — they indicate
  the null fits the moment conditions better than the alternative — so forcing
  `S >= 0` was not justified.

  Capping was also mildly anti-conservative. It placed a point mass at
  `1e-10` (measured: 71% of draws on a Gaussian SV(1) DGP with `T = 1000`,
  `p_null = 1` vs `p_alt = 2`), and the right-tailed Monte Carlo p-value
  `(N + 1 - #{s0 >= sN})/(N + 1)` is algebraically `(1 + #{sN > s0})/(N + 1)`,
  which *excludes* ties from the numerator. When the observed statistic also
  sat at the cap, the p-value was pulled down. On the example above the same
  simulated draws give `p = 0.24` capped versus `p = 0.79` uncapped; the latter
  is correct, since a negative observed statistic is evidence *for* the null.

  Note the cap was provably a no-op whenever `s0 > 1e-10`: capping only moves
  draws upward to the cap and they remain strictly below `s0`, leaving
  `#{sN > s0}` unchanged. Results with a strictly positive observed statistic
  are therefore **unaffected**, as are all results using the default
  `Bartlett = FALSE` statistic `T * sum(phi_extra^2)`, which is non-negative by
  construction and was never capped.

  The other six testing functions (`lmc_t`/`mmc_t`, `lmc_ged`/`mmc_ged`,
  `lmc_lev`/`mmc_lev`) still cap; whether to remove it there is under review.
  Added regression tests to `test-testing.R`.

* Filtering no longer resets the caller's random-number stream as a side
  effect. `fit_ksc_mixture()` (used by `filter_svp(method = "mixture")` for
  Student-t / GED) and the bootstrap particle filter
  (`filter_svp(method = "particle")`, all error types) each call `set.seed()`
  internally to make their own draws reproducible, but previously left the
  global `.Random.seed` in that fixed state on return. Any caller that drew
  random numbers *after* a `filter_svp()` / `svp_IC()` / `svp_AR_order()` /
  `forecast_svp()` call — e.g. a Monte Carlo loop simulating data from the
  ambient stream between statistic evaluations — would get a repeated draw.
  Both sites now save and restore `.Random.seed` around the internal seeding
  (new internal helper `.with_preserved_seed()`), so these functions are
  side-effect-free with respect to the caller's RNG. **All numerical outputs
  (log-likelihoods, IC values, forecasts) are bit-identical** — only the
  side effect is removed. Added regression tests (`test-rng-hygiene.R`).
  *Reproducibility note:* a script that (incorrectly) relied on the ambient
  stream being reset after a filter call will now see different draws than
  under previous versions; seed each replication explicitly.

* GED degrees-of-freedom estimation (`svp(..., errorType = "GED")`): the
  `nu` root-finder now uses the fixed bracket `[0.1, 20]` instead of a
  data-scale-dependent bracket derived from `mean|log(y^2)| / sd|log(y^2)|`.
  Because `|log(y^2)|` is not scale invariant, the old bracket made `nu-hat`
  differ by ~1e-8 under a pure rescaling `y -> c*y` (Brent stopped at a
  different point inside its tolerance ball), even though the estimating
  equation's root is scale invariant. The fixed bracket — already the
  fallback interval inside the C++ solver, and matching the boundary warnings
  — restores exact scale invariance and harmonizes GED with the fixed-bracket
  Student-t solver. Brent's convergence tolerance in both `nu` solvers
  (Student-t and GED) was tightened from `1e-6` to `1e-12`. On real-data
  scales the effect on estimates is below display precision (per-path
  `|Delta nu-hat|` <= ~5e-7; reported means/RMSE shift <= ~2e-8); AR
  coefficients are unaffected. Boundary classification and rates are
  preserved. Added regression tests (`test-ged-nu-bracket.R`).

## Documentation

* `filter_svp()`, `forecast_svp()`, `svp_IC()`, `svp_AR_order()`: corrected
  the `proxy` argument description. `"bayes_optimal"` is the prescribed
  proxy in Ahsan, Dufour and Rodriguez-Rondon (2026) and is the default
  everywhere; `"u"` is a simpler alternative that does not correct the
  `Var(u-hat) > 1` marginal inflation under Student-t leverage. Removed
  outdated "paper-faithful (Remark 3.5)" framing of `"u"` and an
  inconsistent author order (`Rodriguez-Rondon, Dufour and Ahsan` →
  `Ahsan, Dufour and Rodriguez-Rondon`, alphabetical).


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
