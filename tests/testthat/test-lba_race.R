
context("LBA race functions")

test_that("n1CDF corresponds to random derivates", {
  normalised_n1CDF = function(t,...) n1CDF(t,...)/n1CDF(t=Inf,...) 
  r_lba1 <- rlba_norm(1e3, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=0.2)
  r_lba2 <- rlba_norm(1e3, A=0.5, b=1, t0 = 0.5, st0 = 0.2, mean_v=c(1.2, 1), sd_v=0.2)
  expect_less_than(ks.test(r_lba1$rt[r_lba1$response==1], "n1CDF", A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1.0), sd_v=0.2)$p.value, 0.0001)
  expect_less_than(ks.test(r_lba1$rt[r_lba1$response==1], normalised_n1CDF, A=0.6, b=1, t0 = 0.5, mean_v=c(1.2, 1.0), sd_v=0.2)$p.value, 0.0001)
  expect_less_than(ks.test(r_lba2$rt[r_lba2$response==1], normalised_n1CDF, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1.0), sd_v=0.2)$p.value, 0.0001)
  expect_less_than(ks.test(r_lba2$rt[r_lba2$response==1], normalised_n1CDF, A=0.5, b=1, t0 = 0.5, st0 = 0.1, mean_v=c(1.2, 1.0), sd_v=0.2)$p.value, 0.0001)
  
  expect_more_than(ks.test(r_lba1$rt[r_lba1$response==1], normalised_n1CDF, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1.0), sd_v=0.2)$p.value, 0.05)
  expect_more_than(ks.test(r_lba2$rt[r_lba2$response==1], normalised_n1CDF, A=0.5, b=1, t0 = 0.5, st0 = 0.2, mean_v=c(1.2, 1.0), sd_v=0.2)$p.value, 0.05)
})

