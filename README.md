

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

* From CRAN: `install.packages("rtdists")`

* Development version from Github:  
Note, for installing the development version package `devtools` is needed which may require some additional software (see [here](http://r-pkgs.had.co.nz/intro.html) section "Getting started")
`devtools::install_github("rtdists/rtdists")`


## Get Started:
```
require(rtdists)

## LBA: recovers parameters
rt1 <- rLBA(500, A=0.5, b=1, t0 = 0.5, mean_v=c(2.4, 1.6), sd_v=c(1,1.2))
head(rt1)
#          rt response
# 1 0.9618439        1
# 2 1.4018867        1
# 3 1.2188635        2
# 4 1.2192848        2
# 5 1.3690870        1
# 6 0.9478666        2
prop.table(table(rt1$response))
# 
#    1    2 
# 0.67 0.33 

objective_fun <- function(par, rt, response, distribution = "norm") {
  # simple parameters
  spar <- par[!grepl("[12]$", names(par))]  
  
  # distribution parameters:
  dist_par_names <- unique(sub("[12]$", "", grep("[12]$" ,names(par), value = TRUE)))
  dist_par <- vector("list", length = length(dist_par_names))
  names(dist_par) <- dist_par_names
  for (i in dist_par_names) dist_par[[i]] <- as.list(unname(par[grep(i, names(par))]))
  dist_par$sd_v <- c(1, dist_par$sd_v) 

  # get summed log-likelihood:
  d <- do.call(dLBA, args = c(rt=list(rt), response=list(response), spar, dist_par, 
                               distribution=distribution, silent=TRUE))
  if (any(d == 0)) return(1e6)
  else return(-sum(log(d)))
}


objective_fun(c(A=0.5, b=1, t0=0.5, mean_v1=2.4, mean_v2=1.6, sd_v1=1.2), 
              rt=rt1$rt, response=rt1$response)
# [1] -80.07828

init_par <- c(runif(3, 0, 0.5), runif(3, 0.5, 2))
names(init_par) <- c("A", "b", "t0", "mean_v1", "mean_v2", "sd_v2")
nlminb(objective_fun, start = init_par, rt=rt1$rt, response=rt1$response, lower = 0)
# $par
#         A         b        t0   mean_v1   mean_v2     sd_v1 
# 0.6132217 0.9290514 0.5308384 2.1996558 1.6219721 0.8583436 
# 
# $objective
# [1] -84.46802
# 
# $convergence
# [1] 0
# 
# $iterations
# [1] 50
# 
# $evaluations
# function gradient 
#       74      332 
# 
# $message
# [1] "relative convergence (4)"
# 

## Diffusion
example(Diffusion)

```

------
Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.
