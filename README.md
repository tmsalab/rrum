
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rrum

[![Build
Status](https://travis-ci.org/tmsalab/rrum.svg)](https://travis-ci.org/tmsalab/rrum)
[![Package-License](http://img.shields.io/badge/license-GPL%20\(%3E=2\)-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
[![CRAN](http://www.r-pkg.org/badges/version/rrum)](https://cran.r-project.org/package=rrum)
[![Downloads](http://cranlogs.r-pkg.org/badges/rrum?color=brightgreen)](http://www.r-pkg.org/pkg/rrum)

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

## Authors

Steven Andrew Culpepper, Aaron Hudson, and James Joseph Balamuta

## License

GPL (\>= 2)
