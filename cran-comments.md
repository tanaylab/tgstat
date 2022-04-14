## R CMD check results

0 errors | 0 warnings | 0 notes

## Changes

* Fixed compilation on R pre-4.2.0.
* Changed maintainer to Aviezer Lifshitz (aviezer.lifshitz@weizmann.ac.il). Email was sent to CRAN from the previous maintainer (Michael Hoichman, misha@hoichman.com). 
* Use roxygen2 and markdown for documentation.

## rchk notes

The `rchk` notes shown at https://raw.githubusercontent.com/kalibera/cran-checks/master/rchk/results/tgstat.out reflect the fact that protection in tgstat is done via dedicated functions (rprotect and runprotect). This causes `rchk` to shoot false alarams - see https://github.com/kalibera/cran-checks/blob/master/rchk/PROTECT.md#hints-for-interpreting-the-reports.
