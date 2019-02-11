## Test environments

- local OS X install, R 3.5.2
- ubuntu 14.04 (on travis-ci), R 3.5.2
- win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

Possibly mis-spelled words in DESCRIPTION:
  Culpepper (24:65)

Found the following (possibly) invalid URLs:
  URL: http://www.r-pkg.org/pkg/rrum (moved to https://www.r-pkg.org:443/pkg/rrum)
    From: README.md
    Status: 404
    Message: Not Found
  URL: https://cran.r-project.org/web/checks/check_results_rrum.html
    From: README.md
    Status: 404
    Message: Not Found
    
- This is a new release therefore there are no URLs for the package
  on CRAN or at the r-pkg organization. Once accepted, these URLs will
  become active.
- The mis-spelled word is a package author's last name.

## Feedback

We've updated the package based on feedback regarding the `DESCRIPTION` file's
title and description entries. In particular, we've removed the acronym from
the package title and embedded it in the description. Furthermore, we opted to
make the acronym in lowercase.