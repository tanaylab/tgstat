## R CMD check results

0 errors | 0 warnings | 0 notes

## Changes

* Fixed a compilation warning on M1 macs.

## rchk notes

The `rchk` notes shown at https://raw.githubusercontent.com/kalibera/cran-checks/master/rchk/results/tgstat.out reflect the fact that protection in tgstat is done via dedicated functions (rprotect and runprotect). This causes `rchk` to shoot false alarams - see https://github.com/kalibera/cran-checks/blob/master/rchk/PROTECT.md#hints-for-interpreting-the-reports.
