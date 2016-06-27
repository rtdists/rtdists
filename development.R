
require(devtools)
require(testthat)
load_all()

options(error = recover)
options(error = NULL)

build_vignettes()

test()
test_package("rtdists")  # no long tests of n1 functions
test(filter = "bugs")
test(filter = "input")
test(filter = "lba")
test_file("tests/testthat/test-diffusion.R")
test_file("tests/testthat/test-diffusion-math.R")
test_file("tests/testthat/test-ilba_basics.R")
test_file("tests/testthat/test-lba-bugs.R")
test_file("tests/testthat/test-lba-math.R")
test_file("tests/testthat/test-lba_input.R")
test_file("tests/testthat/test-lba_race.R")
test_file("tests/testthat/test-lba_race_input.R")
test_file("tests/testthat/test-rrd.R")

#test_file("tests/testthat/test-lba_race_random_parameters.R")

## Analyze problematic data:

# 1. n1CDF with t0:
load("inst//extdata//n1CDF_diff_example.RData")
n1CDF(r_lba1$rt[ r_lba1$response==1 ], A = A, b = b, t0 = t0, mean_v = v1, sd_v = v2)
n1CDF(r_lba1$rt[ r_lba1$response==1 ] - t0, A = A, b = b, t0 = 0, mean_v = v1, sd_v = v2)



## Complete documentation including DESCRPTION file is written using roxygen2 and wrapper roxyPackage:
require(roxyPackage)  # install.packages("roxyPackage", repo="http://R.reaktanz.de")

R.libs <- "."

roxy.package(
  pck.source.dir = ".",
  pck.version = "0.6-2",
  pck.description = data.frame(
    Package = "rtdists",
    Type = "Package",
    Title = "Response Time Distributions",
    AuthorsR = "c(
        person(given=\"Henrik\", family=\"Singmann\", email=\"singmann+rtdists@gmail.com\", role=c(\"aut\", \"cre\")),
        person(given=\"Scott\", family=\"Brown\", role=c(\"aut\")),
        person(given=\"Matthew\", family=\"Gretton\", role=c(\"aut\")),
        person(given=\"Andrew\", family=\"Heathcote\", role=c(\"aut\")),
        person(given=\"Andreas\", family=\"Voss\", role=c(\"ctb\")),
        person(given=\"Jochen\", family=\"Voss\", role=c(\"ctb\")),
        person(given=\"Andrew\", family=\"Terry\", role=c(\"ctb\"))
    )",
    Depends = "R (>= 3.0.0)",
    Suggests = "testthat, glba, knitr, rmarkdown, dplyr, tidyr, lattice, latticeExtra, binom, RWiener",
    Imports = "evd, msm, gsl, pracma, stats, utils",
    Description = "Provides response time distributions (density/PDF, distribution function/CDF, quantile function, and random generation): (a) Ratcliff diffusion model based on C code by Andreas and Jochen Voss and (b) linear ballistic accumulator (LBA) with different distributions underlying the drift rate.",
    URL = "https://github.com/rtdists/rtdists/",
    License = "GPL (>=3)",
    VignetteBuilder = "knitr",
    stringsAsFactors = FALSE),
  actions = c("roxy"),
  R.libs = R.libs, 
  repo.root = tempdir())

