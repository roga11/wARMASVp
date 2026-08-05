# Replication scripts for the GMKF mixture tables

The files in this folder are the exact scripts that built the mixture tables
shipped in `inst/extdata/` (`gaussian_mixtures.rds`, `ged_mixtures.rds`,
`studentt_mixtures.rds`) — the tables `filter_svp(method = "mixture",
mixture = "optimal")` reads. They are provided so the tables can be rebuilt
and verified independently.

- `Table_GMKF_helpers.R` — shared machinery (target densities, quadrature,
  EM warm starts, trust-region refinement, component tracking, KL
  evaluation). Sourcing it defines functions and nothing else.
- `Table_GMKF_Gaussian.R` — 4 mixtures (K = 7, 10, 12, 16), no shape
  parameter.
- `Table_GMKF_GED.R` — 200 nu nodes x 4 K, nu in [0.1, 20].
- `Table_GMKF_StudentT.R` — 500 nu nodes x 4 K, nu in [2.01, 500].

## How to run

1. **Copy this whole folder to a writable directory first.** Each build
   script writes `logs/` and `csvTables/` and the output `.rds` next to
   itself, unconditionally, as soon as it runs — running it in place would
   write into the installed R library.

   ```r
   dir.create("~/gmkf-rebuild")
   file.copy(list.files(system.file("gmkf-tables", package = "wARMASVp"),
                        full.names = TRUE), "~/gmkf-rebuild")
   ```

2. Run each script with `source()`, or with `Rscript` **from inside the
   folder** (`cd ~/gmkf-rebuild && Rscript Table_GMKF_GED.R`). Under
   `Rscript`, a script that cannot locate itself falls back to the working
   directory to find `Table_GMKF_helpers.R` and stops if it is not there.

3. Run `Table_GMKF_Gaussian.R` first if you want the GED and Student-t
   cross-checks against the Gaussian table (each compares its nu-equivalent
   mixture to the independently built Gaussian one). The check is skipped
   silently when `gaussian_mixtures.rds` is not present in the folder.

## Requirements and runtime

- Packages: `matrixStats` and `clue` (listed in this package's `Suggests`;
  used only by these scripts, not by the package itself).
- The GED and Student-t scripts parallelise over `N_WORKERS = 4` fork-based
  `parallel::mclapply` workers (edit the setting at the top of each script;
  forking is not available on Windows — set `N_WORKERS = 1` there).
  The Gaussian script is sequential.
- Total runtime for all three scripts: roughly 65 minutes (Gaussian ~9 min,
  GED ~21 min, Student-t ~34 min on an Apple-silicon laptop).

## Verifying the shipped tables

Each `.rds` carries a provenance attribute recording the md5 of the build
script and of the helpers file that produced it:

```r
attr(readRDS(system.file("extdata", "ged_mixtures.rds",
                         package = "wARMASVp")), "provenance")
tools::md5sum(system.file("gmkf-tables", "Table_GMKF_GED.R",
                          package = "wARMASVp"))   # must match script_md5
```

A rebuild with these scripts reproduces the shipped tables; compare your
output to `inst/extdata/` with `tools::md5sum()` (byte identity) or
element-wise in R.

## Run-time counterpart

The package evaluates the tables at any nu via cubic-spline interpolation
through the stored component paths — see `R/gmkf_tables.R` in the package
sources (`gmkf_table()` is the exported accessor for the raw tables).
