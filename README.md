
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rrum

[![Build
Status](https://travis-ci.org/tmsalab/rrum.svg)](https://travis-ci.org/tmsalab/rrum)
[![Package-License](http://img.shields.io/badge/license-GPL%20\(%3E=2\)-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
[![CRAN Version
Badge](http://www.r-pkg.org/badges/version/rrum)](https://cran.r-project.org/package=rrum)
[![CRAN
Status](https://cranchecks.info/badges/worst/rrum)](https://cran.r-project.org/web/checks/check_results_rrum.html)
[![Downloads](http://cranlogs.r-pkg.org/badges/rrum?color=brightgreen)](http://www.r-pkg.org/pkg/rrum)
[![Grand Total of
Downloads](http://cranlogs.r-pkg.org/badges/grand-total/rrum?color=brightgreen)](http://www.r-pkg.org/pkg/rrum)
[![Coverage
status](https://codecov.io/gh/tmsalab/rrum/branch/master/graph/badge.svg)](https://codecov.io/github/tmsalab/rrum?branch=master)

The goal of rrum is to provide an implementation of Gibbs sampling
algorithm for Bayesian Estimation of **reduced Reparametrized Unifed
Model (rRUM)**, described by Culpepper and Hudson (2017) \<doi:
10.1177/0146621617707511\>.

## Installation

You can install `rrum` from CRAN using:

``` r
install.packages("rrum")
```

Or, you can be on the cutting-edge development version on GitHub using:

``` r
if(!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("tmsalab/rrum")
```

## Usage

To use `rrum`, load the package using:

``` r
library("rrum")
```

From here, the rRUM model can be estimated using:

``` r
rrum_model = rrum(<data>, <q>)
```

Additional parameters can be accessed with:

``` r
rrum_model = rrum(<data>, <q>, chain_length = 10000L,
                  as = 1, bs = 1, ag = 1, bg = 1,
                  delta0 = rep(1, 2^ncol(Q)))
```

`rRUM` item data can be simulated using:

``` r
rrum_items = simcdm::sim_rrum_items(<q>, <rstar>, <pistar>, <alpha>)
```

## Authors

Steven Andrew Culpepper, Aaron Hudson, and James Joseph Balamuta

## Citing the `rrum` package

To ensure future development of the package, please cite `rrum` package
if used during an analysis or simulation study. Citation information for
the package may be acquired by using in *R*:

``` r
citation("rrum")
```

## License

GPL (\>= 2)
