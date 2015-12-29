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
rt1 <- rlba_norm(500, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3))
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
  for (i in dist_par_names) dist_par[[i]] <- unname(par[grep(i, names(par))])
  
  # get summed log-likelihood:
  d <- do.call(diLBA, args = c(rt=list(rt), response=list(response), spar, dist_par, 
                               distribution=distribution))
  if (any(d == 0)) return(1e6)
  else return(-sum(log(d)))
}

objective_fun(c(A=0.5, b=1, t0=0.5, mean_v1=1.2, mean_v2=1.0, sd_v1=0.2, sd_v2=0.3), 
              rt=rt1$rt, response=rt1$response)
# [1] 37.67205

init_par <- runif(7)
names(init_par) <- c("A", "b", "t0", "mean_v1", "mean_v2", "sd_v1", "sd_v2")
nlminb(objective_fun, start = init_par, rt=rt1$rt, response=rt1$response, lower = 0)
# $par
#         A         b        t0   mean_v1   mean_v2     sd_v1     sd_v2 
# 0.4747780 0.9422244 0.5135954 1.1513527 0.9466506 0.2128113 0.3178318 
# 
# $objective
# [1] 37.07786
# 
# $convergence
# [1] 0
# 
# $iterations
# [1] 42
# 
# $evaluations
# function gradient 
#       61      331 
# 
# $message
# [1] "relative convergence (4)"
# 


## Diffusion
example(Diffusion)

```

