
<!-- README.md is generated from README.Rmd. Please edit that file -->
adapref
=======

<!-- badges: start -->
<!-- badges: end -->
The goal of adapref is to allow adaptive preferential sampling in phylodynamics. It is a package to be used jointly with phylodyn and spmrf. Many of the functions of this package have been developed starting from the work done in these two packages. We acknowledge the developers for their work.

Installation
------------

1.  Install (if necessary) package dependencies and helpers: `devtools`, `ape`, `spmrf`, `rstan`, `loo` and `phylodyn`.

2.  Load `devtools` using `library(devtools)`.

3.  Install `adapref` using

    1.  `install_github("JuliaPalacios/phylodyn")`, or

    2.  `install_github("JuliaPalacios/phylodyn", build_vignettes = TRUE)` if you want some illustrative vignettes (note: using `build_vignettes = TRUE` will make the install take longer).

Vignettes
---------

1.  [SimpleBNPR](https://github.com/mdkarcher/phylodyn/blob/master/vignettes/SimpleBNPR.Rmd): A short example showing how to use BNPR and BNPR-PS on simulated data, illustraring methodology in \[2\] and \[5\].

References
----------
