
require(rtdists)

# data with differences:

fd <- list.files(pattern = "n1CDF_diff_example_[[:digit:]]\\.RData")

strs <- vector("list", length(fd))

for (i in seq_along(strs)) {
  load(fd[[i]])  
  diffs <- n1CDF(r_lba1$rt[ r_lba1$response==1 ], A = A, b = b, t0 = t0, mean_v = v1, sd_v = v2) - n1CDF(r_lba1$rt[ r_lba1$response==1 ] - t0, A = A, b = b, t0 = 0, mean_v = v1, sd_v = v2)
  strs[[i]] <- data.frame(A = A, b = b, t0 = t0, mean_v_1 = v1[1], mean_v_2 = v1[2], sd_v_1 = v2[1], sd_v_2 = v2[2], diff_mean = mean(diffs), sd_diffs = sd(diffs))
}
(diffs <- do.call(rbind, strs))

# data without differences

fnd <- list.files(pattern = "n1CDF_no_diff_example_[[:digit:]]\\.RData")

strs2 <- vector("list", length(fnd))

for (i in seq_along(strs2)) {
  load(fnd[[i]])  
  diffs <- n1CDF(r_lba1$rt[ r_lba1$response==1 ], A = A, b = b, t0 = t0, mean_v = v1, sd_v = v2) - n1CDF(r_lba1$rt[ r_lba1$response==1 ] - t0, A = A, b = b, t0 = 0, mean_v = v1, sd_v = v2)
  strs2[[i]] <- data.frame(A = A, b = b, t0 = t0, mean_v_1 = v1[1], mean_v_2 = v1[2], sd_v_1 = v2[1], sd_v_2 = v2[2], diff_mean = mean(diffs), sd_diffs = sd(diffs))
}
(diffs2 <- do.call(rbind, strs2))

