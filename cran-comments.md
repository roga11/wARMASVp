## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Test environments

* macOS Sonoma 14.2.1, R 4.4.0 (local)

## Notes

* The single NOTE is the expected "New submission" message.

* Slow tests (MMC optimization, particle filter) are wrapped with
  `skip_on_cran()` to keep CRAN check times reasonable.
