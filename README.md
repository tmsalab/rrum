
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rrum

<!-- badges: start -->

[![R-CMD-check](https://github.com/tmsalab/rrum/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tmsalab/rrum/actions/workflows/R-CMD-check.yaml)
[![Package-License](https://img.shields.io/badge/license-GPL%20(%3E=2)-brightgreen.svg?style=flat)](https://www.gnu.org/licenses/gpl-2.0.html)
[![CRAN Version
Badge](https://www.r-pkg.org/badges/version/rrum)](https://cran.r-project.org/package=rrum)
[![CRAN
Status](https://badges.cranchecks.info/worst/rrum.svg)](https://cran.r-project.org/web/checks/check_results_rrum.html)
[![RStudio CRAN Mirror’s Monthly
Downloads](https://cranlogs.r-pkg.org/badges/rrum?color=brightgreen)](https://www.r-pkg.org/pkg/rrum)
[![RStudio CRAN Mirror’s Total
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/rrum?color=brightgreen)](https://www.r-pkg.org/pkg/rrum)
<!-- badges: end -->

The goal of `rrum` is to provide an implementation of Gibbs sampling
algorithm for Bayesian Estimation of **Reduced Reparameterized Unified
Model (rrum)**, described by Culpepper and Hudson (2017) \<doi:
10.1177/0146621617707511\>.

## Installation

You can install `rrum` from CRAN using:

``` r
install.packages("rrum")
```

Or, you can be on the cutting-edge development version on GitHub using:

``` r
# install.packages('remotes')
remotes::install_github("tmsalab/rrum")
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
# Set a seed for reproducibility
set.seed(888)

# Setup Parameters
N = 15   # Number of Examinees / Subjects
J = 10   # Number of Items
K = 2    # Number of Skills / Attributes

# Simulate identifiable Q matrix
Q = sim_q_matrix(J, K)

# Penalties for failing to have each of the required attributes
rstar  = .5 * Q

# The probabilities of answering each item correctly for individuals 
# who do not lack any required attribute
pistar = rep(.9, J)

# Latent Class Probabilities
pis = c(.1, .2, .3, .4)

# Generate latent attribute profile with custom probability (N subjects by K skills)
subject_alphas = sim_subject_attributes(N, K, prob = pis)

# Simulate rrum items
rrum_items = simcdm::sim_rrum_items(Q, rstar, pistar, subject_alphas)
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
