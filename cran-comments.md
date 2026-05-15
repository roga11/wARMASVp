## Submission

This is a version update (0.1.0 -> 0.2.0). It adds information-criterion
AR-order selection (`svp_IC()`, `svp_AR_order()`), extends the AR-order tests
to Student-t and GED innovations, fixes some filtering/forecasting bugs
under leverage, and ports the mixture-filter EM step to C++ for speed. See
NEWS.md for the full list. `sim_svp()` now always returns a named list (a
documented breaking change).

## R CMD check results

0 errors | 0 warnings | 1 note

* The NOTE is the expected "New maintainer" message. The maintainer name has
  changed from "Gabriel Rodriguez Rondon" to "Gabriel Rodriguez-Rondon"
  (hyphenated form). The email address and ORCID are unchanged; this is the
  same person, with the correct spelling of the surname.

## Test environments

* macOS Sonoma 14.2.1, R 4.4.0 (local)
* Windows Server 2022, R release (win-builder)
* Windows Server 2022, R devel (win-builder)

## Reverse dependencies

There are no reverse dependencies.

## Notes

* Slow tests (MMC optimization, particle filter) are wrapped with
  `skip_on_cran()` to keep CRAN check times reasonable.
