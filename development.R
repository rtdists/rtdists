
require(devtools)
load_all()

rlba_norm(10, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2, 0.4, 0.5))
rlba_gamma(10, A=0.5, b=1, t0 = 0.5, shape_v=c(1, 1.6, 0.7), scale_v=c(0.2,0.3))
rlba_frechet(10, A=0.5, b=1, t0 = 0.5, shape_v=c(1.2, 1), scale_v=c(0.2,0.3,0.8))
rlba_lnorm(10, A=0.5, b=1, t0 = 0.5, meanlog_v=c(1.2, 1), sdlog_v=c(0.2, 0.3,0.4))

dlba_frechet(seq(0.6,2,0.1), A=A, b=b, t0=t0, shape_v=1.1, scale_v=1)

rlba_norm(100, b = 1.3, A = 1, vs = c(0.8, 1.2), s = 1.2, t0 = .2, st0 = 0)

require(testthat)
test_package("rtdists")


## Complete documentation including DESCRPTION file is written using roxygen2 and wrapper roxyPackage:
require(roxyPackage)  # install.packages("roxyPackage", repo="http://R.reaktanz.de")

R.libs <- "."

roxy.package(
  pck.source.dir = ".",
  pck.version = "0.2-1",
  pck.description = data.frame(
    Package = "rtdists",
    Type = "Package",
    Title = "Respone Time distribtuions in R",
    AuthorsR = "c(
        person(given=\"Scott\", family=\"Brown\", role=c(\"aut\")),
        person(given=\"Matthew\", family=\"Gretton\", role=c(\"aut\")),
        person(given=\"Andrew\", family=\"Heathcote\", role=c(\"aut\")),
        person(given=\"Andreas\", family=\"Voss\", role=c(\"aut\")),
        person(given=\"Jochen\", family=\"Voss\", role=c(\"aut\")),
        person(given=\"Henrik\", family=\"Singmann\", email=\"singmann+rtdists@gmail.com\", role=c(\"aut\", \"cre\"))
    )",
    Depends = "R (>= 3.0.0)",
    Suggests = "testthat",
    Imports = "evd, msm, gsl",
    Description = "Provides response time distributions: (a) diffusion model based on C code by Andreas and Jochen Voss and (b) linear ballistic accumulator (LBA) with different distribution underlying the drift rate variability.",
    URL = "https://github.com/rtdists/rtdists/",
    License = "GPL (>=3)",
    stringsAsFactors = FALSE),
  actions = c("roxy"),
  R.libs = R.libs, 
  repo.root = tempdir())

