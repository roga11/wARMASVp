## Submission

This is a version update (0.2.0 -> 0.3.0). It adds precomputed mixture tables
for the Gaussian-mixture Kalman filter (new `mixture` argument, new
`gmkf_table()` accessor), an estimated measurement-noise variance for Gaussian
fits (`svp(..., sigeps = "estimated")`), and a number of fixes to filtering
under leverage and to the Monte Carlo tests. See NEWS.md for the full list.

User-visible changes that can alter existing results or calls:

* The `proxy` argument of `filter_svp()`, `forecast_svp()`, `svp_IC()` and
  `svp_AR_order()` has been removed.
* `del` now defaults to `0` (was `1e-10`); series containing exact zeros
  require an explicit positive `del`.
* The mixture filter defaults to the precomputed tables (`mixture =
  "optimal"`); `mixture = "em"` gives the previous behaviour.
* The `mmc_*` functions now return objects of class `"svp_test"`.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* macOS Sonoma 14.2.1, R 4.4.0 (local)
* Windows Server 2022, R release (win-builder) -- TO BE RUN BEFORE SUBMISSION
* Windows Server 2022, R devel (win-builder) -- TO BE RUN BEFORE SUBMISSION

## Reverse dependencies

There are no reverse dependencies.

## Notes

* Slow tests (MMC optimization, particle filter) are wrapped with
  `skip_on_cran()` to keep CRAN check times reasonable.
