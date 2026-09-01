# rtdists

## 0.12-x

Released September 2026

- Added racing diffusion model (RDM), see `?RDM`. Thanks to Kianté Fernandez (#23).
- Fixed a rather large bug in `pdiffusion` that occurred whenever RTs were not ordered. Thanks to Kianté Fernandez (#21).
- Large speed improvement to diffusion model functions. Change is probably most noticeable in `ddiffusion` when called with large number of RTs and varying parameters. Thanks to Kianté Fernandez (#21).
- Fixed bugs in the lognormal and gamma single-accumulator LBA functions (`dlba_lnorm`, `plba_lnorm`, `dlba_gamma`, `plba_gamma`) that produced incorrect results when any element of `A` was below `1e-10` (the A_small branch). The lognormal functions used `max` instead of `min` in several intermediate calculations, `plba_lnorm` did not apply `rem_t0` in the A_small branch, and the gamma functions used full-length vectors instead of subsetted ones.

## 0.11-6

Released June 2026.

- Included a check in `rdiffusion` preventing RTs < `t0`. Thanks to Andrew Heathcote for reporting this.

## 0.11-x

Released March 2020.

- Faster calculation of diffusion density (`ddiffusion`) and CDF (`pdiffusion`) in case exactly one set of parameters is passed. In this case, some checks are skipped now leading to a 40% to 50% speed increase.
- Diffusion quantile function `qdiffusion` should be faster now when it is called with many probabilities for the same response.
- Diffusion quantile function `qdiffusion` has a new argument `max_diff` which allows to control the minimally acceptable difference between desired and obtained probabilities.
- Added experimental `method = "qdiffusion"` to `rdiffusion()` which obtains random derivates via the quantile function and `runif()`. This method is for the time being a lot slower than the native one (i.e., `method = "fastdm"`).
- Tests should work even if `options("stringsAsFactors" = FALSE)` (in preparation for R 4.0).
- 0.11-3: removed bug that diffusion functions failed if there were NAs in the parameters.
- 0.11-3: uses RCPP `STRICT_R_HEADERS`: https://github.com/rtdists/rtdists/pull/14

## 0.10-x

Released October 2019.

- Removed another argument passing bug that occured when using `args.dist` and more than two drift rates. Reported by Glen Livingston Jr.
- Removed issue in vignette due to breaking change in a tidyverse package.

## 0.9-x

Released August 2018.

- Removed bug in `dLBA` and `pLBA` that prevented correct usage of trial-wise parameters. This bug always appeared when data with more than one response was present together with trial-wise parameters. In the following example the first call should be identical to the second and third call:

    ```r
    x1 <- dLBA(rt=c(1,1), response=c(1,2), A=1,b=list(c(1,3),c(2,4)),
            t0=0.1, mean_v=c(3,3), sd_v=c(1,1),distribution="norm")
    x2a <- dLBA(rt=c(1), response=c(1), A=1,b=list(c(1),c(2)),
            t0=0.1,mean_v=c(3,3),sd_v=c(1,1),distribution="norm")
    x2b <- dLBA(rt=c(1), response=c(2), A=1,b=list(c(3),c(4)),
            t0=0.1,mean_v=c(3,3),sd_v=c(1,1),distribution="norm")
    all(x1 == c(x2a, x2b)) ## should be TRUE
    ```

## 0.8-x

Released December 2017 & updated June 2018.

- Removed bug preventing `args.list` (e.g., `posdrift`) to be passed correctly in `dLBA`, `pLBA`, `qLBA`, and `rLBA`. Thanks to Bruno Nicenboim for reporting this. See: https://github.com/rtdists/rtdists/issues/7
- Removed an non-used `devtools::` call in tests that caused a false positive CRAN warning (June 2018).
- Deactivated a few more tests in CRAN for faster checking (June 2018).

## 0.7-x

Released May 2017.

- Performance of diffusion functions and `rLBA` increased, especially for calls with parameters that differ trialwise. As a consequence single `rlba_...` functions now return a matrix and no `data.frame`.
- All C functions are now accessed via Rcpp.
- `pdiffusion` uses the C++ CDF (no more numerical integration in R).
- `sv` can produce slow errors, and `sz` fast erros (this was the wrong way around in the documentation). Thanks to Gabriel Tillman for noticing that.
- Removed a bug in the `pdiffusion` C code letting rtdists hang indefinitely (see https://github.com/rtdists/rtdists/pull/3). Thanks to Tomas Kalibera for the fix.
- Removed bug: `meanlog_v` and `sdlog_v` were not recycled for the lnorm LBA.
- Ratcliff and Rouder (1998) vignette now uses nested `data.frame`s and `purrr::map` (i.e., more proper use of the tidyverse).
- Removed bug when non-accumulator LBA parameters where passed as a named list.

## 0.6-6

Bug-fix version, released July 2016.

- Bug when passing start point with `s != 1` removed.

## 0.6-x

Released July 2016.

- Start point `z` in diffusion model is now on absolute scale and not relative to `b` in line with `A` (start point of LBA) which is also on absolute scale. (Thanks to Steve Lewandowsky for noticing this.)
- PDFs, CDFs, and quantile functions of both models now accept a `data.frame` as first argument containing both RTs/probabilities and responses. Allows more convenient way to pass data.
- Renamed `boundary` (argument in diffusion functions) to `response` to be in line with LBA functions. (Thanks to Steve Lewandowsky for suggesting this.)
- Added diffusion constant `s` as argument to all diffusion functions.
- Added `scale_p` and `scale_max` arguments to quantile functions which automatically scale the entered probability to more conveniently obtain predicted quantiles.
- LBA functions now accept `factor` as response (which is converted via `as.numeric`). This allows to pass results from `rdiffusion` directly to LBA function.
- Changed integration routine in `pdiffusion` to `pracma::integral()` which seems to be more robust. (Thanks to Anna-Lena Schubert and Steve Lewandowsky for reporting problems with the previous version.)
- Removed bug preventing `lnorm` as distribution in `rLBA`. (Thanks to Steve Lewandowsky for reporting this bug.)

## 0.5-x

Released May 2016.

- Calculation of the CDF for the diffusion model was incorrect (this bug was present in all prior versions of rtdists). `pdiffusion` now simply integrates the PDF to obtain the CDF using R's `integrate` which provides the correct result (albeit slower).
- Added `rr98` data set: Experiment 1 from Ratcliff and Rouder (1998, Psych. Science). We thank Roger Ratcliff and Jeff Rouder for providing the data.
- Added vignette showing how to analyze the data from Ratcliff and Rouder (1998) with both diffusion and LBA model.
- Quantile functions work more robust and try `uniroot` if `optimize` does not converge.

## 0.4-x

Released April 2016.

- Added `dLBA()`, `pLBA()`, `qLBA()`, and `rLBA()`. `dLBA()` is a fully vectorized versions of `n1PDF` which has response as second argument, allowing to get the density for each response and corresponding response time in one step. As for the diffusion model (see below), this allows a likelihood function which only includes one call to the density function. `pLBA()` and `qLBA()` are the correpsonding CDF and quantile functions, respectively. `rLBA()` is a fully vectorized version of the RNG functions and should be used from now on as top-level function.
- `t0` in the LBA now accepts accumulator and trialwise parameters just as `A` and `b`. `st0` now accepts trialwise parameter (not accumulator wise).
- Diffusion model function have been renamed to `ddiffusion`, `pdiffusion`, and `rdiffusion`. Added quantile function for diffusion model, `qdiffusion`.
- Diffusion model functions are now completely vectorized and accept vectors as parameters (including for boundary). As for the LBA, this allows a likelihood function which only includes one call to the density function (see examples).
- Boundary parameter for diffusion functions accept numeric/factor vectors.
- `t0` in the diffusion model now corresponds to the lower bound of the uniform distribution from which `t0` is drawn (it was the mean before). The specifications of `t0` now agree between LBA and diffusion model.
- Density for diffusion model is now always positive.
- First argument in most functions (vector of response times) renamed to `rt` (was `t` before).
- PDF and CDF LBA functions more robust (now mostly return 0 instead of NaN) for problematic parameter values (such as `A = 0`) [2015-09-17].
