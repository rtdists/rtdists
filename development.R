
require(devtools)
require(testthat)
load_all()

options(error = recover)
options(error = NULL)

build()

build_vignettes()

document()

Rcpp::compileAttributes()

test()
test_package("rtdists")  # no long tests of n1 functions
test(filter = "bugs")
test(filter = "input")
test(filter = "lba")
test(filter = "rng")
test(filter = "diffusion")
test_file("tests/testthat/test-diffusion.R")
test_file("tests/testthat/test-diffusion-math.R")
test_file("tests/testthat/test-diffusion-bugs.R")
test_file("tests/testthat/test-lba_basics.R")
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


## reverse dependency checks
revdepcheck::revdep_check(num_workers = 2) ## run in a new terminal
revdepcheck::revdep_summary()

###
rhub::list_validated_emails()
rhub::validate_email("singmann@gmail.com", 
                     token = "2928445ca3c645079033bde73b754e40")
rhub::check_for_cran("../rtdists_0.11-5.tar.gz")
rhub::check_on_debian(path = "../rtdists_0.11-4.tar.gz", 
                        check_args = "--as-cran")
rhub::check_with_rdevel(path = "../rtdists_0.11-3.tar.gz", 
                        check_args = "--as-cran")
rhub::check_with_sanitizers("../rtdists_0.11-4.tar.gz")
rhub::check(platform = "linux-x86_64-rocker-gcc-san",
            path = "../rtdists_0.11-5.tar.gz", 
            check_args = "--as-cran")
