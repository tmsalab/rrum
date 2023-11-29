# rrum 0.2.2

## Bugfixes

- Addressed two issues in R-devel related to C++ usage (deprecation of `CXX_STD` and attribute compilation woes). ([#15](https://github.com/tmsalab/rrum/issues/15))
- Changed the format in `CITATION` to use `c()` instead of `as.personList()` and
  switched entries from `citEntry()` to `bibentry()`.

## Deployment

- Updated GitHub Actions for testing, pkgdown, and coverage ([#16](https://github.com/tmsalab/rrum/pull/16)).

# rrum 0.2.1

## Documentation

- Added a `pkgdown` website that deploys to <https://tmsalab.github.io/rrum/> ([#12](https://github.com/tmsalab/rrum/pull/12), [#14](https://github.com/tmsalab/rrum/pull/14)).
- Pruned non-user facing documentation entry ([#14](https://github.com/tmsalab/rrum/pull/14)).

## Deployment

- Changed from Travis-CI to GitHub Actions for R ([#13](https://github.com/tmsalab/rrum/pull/13), [#14](https://github.com/tmsalab/rrum/pull/14))

# rrum 0.2.0

## API Breakage

- Deprecated `rRUM_Gibbs()` in favor of `rrum()`.
- Deprecated `pi_reference()` in favor of `simcdm::attribute_classes()`. 

## Changes

- Added `CITATION` file for citing both the APM paper and package.
- Imported simulation functions from `simcdm`

## Documentation

- Improved `README` examples

## Deployment

- Added Travis-CI configuration for TMSA Lab.
- Added Unit Tests for model reproducibility.
- Added code coverage checks.

# rrum 0.1.0

- Improved documentation
- Addressed RcppExport updates

# rrum 0.0.5

- Introduced new rRUM estimation routine.
- Provided a means to simulate rRUM data.

