## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Test environments

* macOS Sonoma 14.2.1, R 4.4.0 (local)
* Windows Server 2022, R 4.5.3 release (win-builder)
* Windows Server 2022, R 4.6.0 beta (win-builder)

## Notes

* The single NOTE is the expected "New submission" message, along with
  possibly misspelled words (Ahsan, Dufour, Kalman, LMC, MMC, Rondon,
  SV, Winsorized) which are author names and standard technical terms
  in the stochastic volatility literature.

* Slow tests (MMC optimization, particle filter) are wrapped with
  `skip_on_cran()` to keep CRAN check times reasonable.
