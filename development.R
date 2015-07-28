
require(devtools)
load_all()

options(error = recover)
options(error = NULL)

rrd(10, a=1, z=0.5, v=2, t0=0.5, d=0, sz=0, sv=0, st0=0)
rt1 <- rrd(10, a=c(1, 1.5, 1.2), z=0.5, v=1, t0=0.5, d=0, sz=0, sv=0, st0=0)
drd(rt1[rt1$response == "upper", "rt"], a=c(1, 1.5, 1.2), z=0.5, v=1, t0=0.5, d=0, sz=0, sv=0, st0=0)
prd(sort(rt1[rt1$response == "upper", "rt"]), a=c(1, 1.5, 1.2), z=0.5, v=1, t0=0.5, d=0, sz=0, sv=0, st0=0)

rlba_norm(10, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2, 0.4, 0.5))
rlba_gamma(10, A=0.5, b=1, t0 = 0.5, shape_v=c(1, 1.6, 0.7), scale_v=c(0.2,0.3))
rlba_frechet(10, A=0.5, b=1, t0 = 0.5, shape_v=c(1.2, 1), scale_v=c(0.2,0.3,0.8))
rlba_lnorm(10, A=0.5, b=1, t0 = 0.5, meanlog_v=c(1.2, 1), sdlog_v=c(0.2, 0.3,0.4))

dlba_frechet(seq(0.6,2,0.1), A=A, b=b, t0=t0, shape_v=1.1, scale_v=1)

rlba_norm(100, b = 1.3, A = 1, vs = c(0.8, 1.2), s = 1.2, t0 = .2, st0 = 0)

(x <- rlba_frechet(10, A=0.7, b = 0.5, t0 = 0.1, shape_v = c(0.7, 1.0), scale_v = c(0.7, 1.2)))

require(testthat)
test()
test_package("rtdists")  # no long tests of n1 functions
test_file("tests/testthat/test-lba-math.R")

#x <- .Random.seed
test_file("tests/testthat/test-lba_race.R")

.Random.seed <- x

test_file("tests/testthat/test-lba_race_random_parameters.R")

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
  pck.version = "0.4-4",
  pck.description = data.frame(
    Package = "rtdists",
    Type = "Package",
    Title = "Response Time Distributions",
    AuthorsR = "c(
        person(given=\"Scott\", family=\"Brown\", role=c(\"aut\")),
        person(given=\"Matthew\", family=\"Gretton\", role=c(\"aut\")),
        person(given=\"Andrew\", family=\"Heathcote\", role=c(\"aut\")),
        person(given=\"Andreas\", family=\"Voss\", role=c(\"ctb\")),
        person(given=\"Jochen\", family=\"Voss\", role=c(\"ctb\")),
        person(given=\"Andrew\", family=\"Terry\", role=c(\"ctb\")),
        person(given=\"Henrik\", family=\"Singmann\", email=\"singmann+rtdists@gmail.com\", role=c(\"aut\", \"cre\"))
    )",
    Depends = "R (>= 3.0.0)",
    Suggests = "testthat",
    Imports = "evd, msm, gsl, stats, utils",
    Description = "Provides response time distributions (density/PDF, distribution function/CDF, and random generation): (a) Ratcliff diffusion model based on C code by Andreas and Jochen Voss and (b) linear ballistic accumulator (LBA) with different distribution underlying the drift rate.",
    URL = "https://github.com/rtdists/rtdists/",
    License = "GPL (>=3)",
    stringsAsFactors = FALSE),
  actions = c("roxy"),
  R.libs = R.libs, 
  repo.root = tempdir())

