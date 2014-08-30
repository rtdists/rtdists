
require(devtools)
require(roxygen2)
load_all()

rrd(10, c(rep(1, 4), rep(0.1, 4)))

require(testthat)
test_package("rtdists")


## Complete documentation including DESCRPTION file is written using roxygen2 and wrapper roxyPackage:
require(roxyPackage)  # install.packages("roxyPackage", repo="http://R.reaktanz.de")

R.libs <- "."

roxy.package(
  pck.source.dir = ".",
  pck.version = "0.1-0",
  pck.description = data.frame(
    Package = "rtdists",
    Type = "Package",
    Title = "Respone Time distribtuions in R",
    AuthorsR = "c(
        person(given=\"Andrew\", family=\"Heathcote\", role=c(\"aut\")),
        person(given=\"Matthew\", family=\"Gretton\", role=c(\"aut\")),
        person(given=\"Andreas\", family=\"Voss\", role=c(\"aut\")),
        person(given=\"Henrik\", family=\"Singmann\", email=\"singmann+rtdists@gmail.com\", role=c(\"aut\", \"cre\"))
    )",
    Depends = "R (>= 3.0.0)",
    Suggests = "testthat",
    #Imports = "",
    Description = "Provides response time distributions based on C code by Andreas Voss.",
    URL = "https://github.com/rtdists/rtdists/",
    License = "GPL (>=3)",
    stringsAsFactors = FALSE),
  actions = c("roxy"),
  R.libs = R.libs, 
  repo.root = tempdir())

