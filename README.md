rtdists: Response time distributions in R
====

## Features

* Distribution functions for response time models:
  * distribution function or probability density function (PDF)
  * cumulative distribution function (CDF)
  * random number generation (RNG)

* The following choice RT models are currently supported:
  * diffusion model (using fast-dm C code by Andreas and Jochen Voss)
  * LBA with varying distributions of the drift rate variability

## Installation

* From CRAN: In the future

* Development version from Github:  
Note, for installing the development version package `devtools` is needed which may require some additional software (see [here](http://r-pkgs.had.co.nz/intro.html) section "Getting started")
```
devtools::install_github("rtdists/rtdists")
```

## Get Started:
```
require(rtdists)
rrd(10, a=1, z=0.5, v=2, t0=0.5, d=0, sz=0, sv=0, st0=0)
example(Diffusion)

rlba_norm(10, A=0.5, b=1, v=c(1.2, 1), sv=c(1,1.5), t0 = 0.5)
rlba_gamma(10, A=0.5, b=1, v=c(1.2, 1), sv=c(1, 1.5), t0 = 0.5)
rlba_frechet(10, A=0.5, b=1, v=c(1.2, 1), sv=c(1, 1.5), t0 = 0.5)
rlba_lnorm(10, A=0.5, b=1, v=c(1.2, 1), sv=c(1, 1.5), t0 = 0.5)
example(LBA)
```

