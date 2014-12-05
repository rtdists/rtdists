context("LBA-math agrees with current implementation")

runif(1)
x <- .Random.seed
set.seed(2)

test_that("PDF and CDF", {
  n <- 10
  samples_per_run <- 100
  source(system.file("extdata", "lba-math.R", package = "rtdists"))
        
  #source("inst/extdata//lba-math.r")
  for (i in seq_len(n)) {
    A <- runif(1, 0.3, 0.9)
    b <- A+runif(1, 0, 0.5)
    t0 <- runif(1, 0.1, 0.7)
    v1 <- runif(2, 0.5, 1.5)
    v2 <- runif(2, 0.1, 0.5)
    r_lba1 <- rlba_norm(samples_per_run, A=A, b=b, t0 = t0, mean_v=v1[1:2], sd_v=v2[1:2])
    
    expect_equal(
      dlba_norm(r_lba1$rt[r_lba1$response==1], A=A, b=b, t0 = t0, mean_v=v1[1], sd_v=v2[1]), 
      fptpdf(pmax(r_lba1$rt[r_lba1$response==1]-t0[1], 0), x0max=A, chi=b, driftrate=v1[1], sddrift=v2[1])
      )
    expect_equal(
      plba_norm(r_lba1$rt[r_lba1$response==1], A=A, b=b, t0 = t0, mean_v=v1[1], sd_v=v2[1]), 
      fptcdf(pmax(r_lba1$rt[r_lba1$response==1]-t0[1], 0), x0max=A, chi=b, driftrate=v1[1], sddrift=v2[1])
    )
    
  }
  
})

test_that("small A values for 'norm'", {
  n <- 10
  samples_per_run <- 100
  source(system.file("extdata", "lba-math.R", package = "rtdists"))
  
  #source("inst/extdata//lba-math.r")
  for (i in seq_len(n)) {
    A <- runif(1, 0, 1e-10)
    b <- A+runif(1, 0, 0.5)
    t0 <- runif(1, 0.1, 0.7)
    v1 <- runif(2, 0.5, 1.5)
    v2 <- runif(2, 0.1, 0.5)
    r_lba1 <- rlba_norm(samples_per_run, A=A, b=b, t0 = t0, mean_v=v1[1:2], sd_v=v2[1:2])
    
    expect_equal(
      dlba_norm(r_lba1$rt[r_lba1$response==1], A=A, b=b, t0 = t0, mean_v=v1[1], sd_v=v2[1]), 
      fptpdf(pmax(r_lba1$rt[r_lba1$response==1]-t0[1], 0), x0max=A, chi=b, driftrate=v1[1], sddrift=v2[1])
    )
    expect_equal(
      plba_norm(r_lba1$rt[r_lba1$response==1], A=A, b=b, t0 = t0, mean_v=v1[1], sd_v=v2[1]), 
      fptcdf(pmax(r_lba1$rt[r_lba1$response==1]-t0[1], 0), x0max=A, chi=b, driftrate=v1[1], sddrift=v2[1])
    )
    
  }
  
})

test_that("Random generation", {
  n <- 10
  samples_per_run <- 100
  source(system.file("extdata", "lba-math.R", package = "rtdists"))
  #source("inst/extdata//lba-math.r")
  for (i in seq_len(n)) {
    A <- runif(1, 0.3, 0.9)
    b <- A+runif(1, 0, 0.5)
    t0 <- runif(1, 0.1, 0.7)
    v1 <- runif(2, 0.5, 1.5)
    v2 <- runif(2, 0.1, 0.5)
    
    x <- .Random.seed   
    
    r_lba1 <- rlba_norm(samples_per_run, A=A, b=b, t0 = t0, mean_v=v1[1:2], sd_v=v2[1:2])
    
    .Random.seed <<- x
    
    #r_lba2 <- rlba_norm(samples_per_run, A=A, b=b, t0 = t0, mean_v=v1[1:2], sd_v=v2[1:2])
    r_lba2 <- rlba(samples_per_run, A=A, b=b, t0 = t0, vs=v1[1:2], s=v2[1:2])
    
    expect_equal(r_lba1$rt, r_lba2$rt)
    expect_equal(r_lba1$resp, r_lba2$resp)
  }
})

test_that("n1CDF", {
  n <- 10
  samples_per_run <- 100
  source(system.file("extdata", "lba-math.R", package = "rtdists"))
  #source("inst/extdata//lba-math.r")
  for (i in seq_len(n)) {
    A <- runif(1, 0.3, 0.9)
    b <- A+runif(1, 0, 0.5)
    t0 <- runif(1, 0.1, 0.7)
    v1 <- runif(2, 0.5, 1.5)
    v2 <- runif(2, 0.1, 0.5)
    r_lba1 <- rlba_norm(samples_per_run, A=A, b=b, t0 = t0, mean_v=v1[1:2], sd_v=v2[1:2], posdrift = TRUE)
    #head(r_lba1)
    #if(!isTRUE(all.equal(n1CDF(r_lba1$rt[r_lba1$response==1], A=A, b=b, t0 = t0, mean_v=v1[1:2], sd_v=v2[1]),.n1CDF(pmax(r_lba1$rt[r_lba1$response==1]-t0[1], 0), x0max=A, chi=b, drift=v1[1:2], sdI=v2[1]) ))) browser()
      #n1CDF(r_lba1$rt[r_lba1$response==1], A=A, b=b, t0 = t0, mean_v=v1[1:2], sd_v=v2[1], browser = TRUE) 
      #n1CDF(pmax(r_lba1$rt[r_lba1$response==1]-t0[1], 0), A=A, b=b, t0 = 0, mean_v=v1[1:2], sd_v=v2[1])      
      #.n1CDF(pmax(r_lba1$rt[r_lba1$response==1]-t0[1], 0), x0max=A, chi=b, drift=v1[1:2], sdI=v2[1], browser=TRUE)     
    #save(r_lba1, A, b, t0, v1, v2, file = "n1CDF_no_diff_example_5.RData")
    
    expect_equal(
      n1CDF(r_lba1$rt[r_lba1$response==1], A=A, b=b, t0 = t0, mean_v=v1[1:2], sd_v=v2[1]), 
      .n1CDF(pmax(r_lba1$rt[r_lba1$response==1]-t0[1], 0), x0max=A, chi=b, drift=v1[1:2], sdI=v2[1])
    )
    
  }
  
})

.Random.seed <<- x
