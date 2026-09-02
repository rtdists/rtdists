# Examples where ddiffusion needs higher precision than the default (3)
# to match WienR's WienerPDF output.
#
# These cases arise when variability parameters (st0, sz) are non-zero and
# the starting point is near a boundary, making the numerical integration
# over the variability distributions challenging for rtdists' fast-dm C code.

library(rtdists)
library(WienR)

# Example 1: st0 with starting point near lower boundary
# z = 0.15 (absolute), a = 0.5 => w = 0.3 (close to lower boundary)
# st0 = 0.16 => non-decision time uniform on [0.14, 0.30]
# The lower-boundary density is large and rapidly changing, requiring
# precision >= 5 to match WienR.
ex1 <- list(a = 0.5, v = 0.5, t0 = 0.14, z = 0.15, sv = 0, sz = 0, st0 = 0.16)
# Lower-boundary density at rt = 0.3
d_wien <- WienerPDF(t = 0.3, response = "lower", a = ex1$a, v = ex1$v,
                    w = ex1$z / ex1$a, t0 = ex1$t0, st0 = ex1$st0)$value
for (p in c(3, 4, 5, 6, 7, 8)) {
  d_rtd <- ddiffusion(0.3, "lower", a = ex1$a, v = ex1$v, t0 = ex1$t0,
                      z = ex1$z, st0 = ex1$st0, precision = p)
  cat(sprintf("    precision=%d: rtdists=%.10f  diff=%.6e\n", p, d_rtd, abs(d_rtd - d_wien)))
}
cat(sprintf("    WienR (target):           %.10f\n\n", d_wien))

# Example 2: sz + st0 with starting point near lower boundary
# z = 0.15 (absolute), a = 0.5 => w = 0.3
# sz = 0.25 => sw = 0.5, starting point range U(0.05, 0.55)
# st0 = 0.16 => non-decision time uniform on [0.14, 0.30]
# Double integration over both sz and st0 near the boundary requires
# precision >= 8 to fully converge.
ex2 <- list(a = 0.5, v = 0.5, t0 = 0.14, z = 0.15, sv = 0, sz = 0.25, st0 = 0.16)
# Lower-boundary density at rt = 0.3
d_wien <- WienerPDF(t = 0.3, response = "lower", a = ex2$a, v = ex2$v,
                    w = ex2$z / ex2$a, t0 = ex2$t0, sw = ex2$sz / ex2$a, st0 = ex2$st0)$value
for (p in c(3, 4, 5, 6, 7, 8, 9, 10)) {
  d_rtd <- ddiffusion(0.3, "lower", a = ex2$a, v = ex2$v, t0 = ex2$t0,
                      z = ex2$z, sz = ex2$sz, st0 = ex2$st0, precision = p)
  cat(sprintf("    precision=%2d: rtdists=%.10f  diff=%.6e\n", p, d_rtd, abs(d_rtd - d_wien)))
}
cat(sprintf("    WienR (target):           %.10f\n\n", d_wien))

# Example 3: sz + st0, upper boundary (same params, for contrast)
# Same parameters as Example 2, but upper boundary.
# The upper-boundary density is smoother and converges at precision = 5.
# Upper-boundary density at rt = 0.3
d_wien <- WienerPDF(t = 0.3, response = "upper", a = ex2$a, v = ex2$v,
                    w = ex2$z / ex2$a, t0 = ex2$t0, sw = ex2$sz / ex2$a, st0 = ex2$st0)$value
for (p in c(3, 4, 5, 6, 7, 8)) {
  d_rtd <- ddiffusion(0.3, "upper", a = ex2$a, v = ex2$v, t0 = ex2$t0,
                      z = ex2$z, sz = ex2$sz, st0 = ex2$st0, precision = p)
  cat(sprintf("    precision=%2d: rtdists=%.10f  diff=%.6e\n", p, d_rtd, abs(d_rtd - d_wien)))
}
cat(sprintf("    WienR (target):           %.10f\n\n", d_wien))

# Example 4: sv + st0 with starting point near lower boundary
# z = 0.15 (absolute), a = 0.5 => w = 0.3
# sv = 0.3, st0 = 0.16
# Drift variability plus non-decision time variability near the boundary.
ex4 <- list(a = 0.5, v = 0.5, t0 = 0.14, z = 0.15, sv = 0.3, sz = 0, st0 = 0.16)
# Lower-boundary density at rt = 0.3
d_wien <- WienerPDF(t = 0.3, response = "lower", a = ex4$a, v = ex4$v,
                    w = ex4$z / ex4$a, t0 = ex4$t0, sv = ex4$sv, st0 = ex4$st0)$value
for (p in c(3, 4, 5, 6, 7, 8)) {
  d_rtd <- ddiffusion(0.3, "lower", a = ex4$a, v = ex4$v, t0 = ex4$t0,
                      z = ex4$z, sv = ex4$sv, st0 = ex4$st0, precision = p)
  cat(sprintf("    precision=%2d: rtdists=%.10f  diff=%.6e\n", p, d_rtd, abs(d_rtd - d_wien)))
}
cat(sprintf("    WienR (target):           %.10f\n\n", d_wien))

# Example 5: st0 with starting point near lower boundary, small effective rt
# z = 0.15 (absolute), a = 0.5 => w = 0.3 (close to lower boundary)
# st0 = 0.1 => non-decision time uniform on [0.275, 0.375]
# At rt = 0.3, the effective decision time is only ~0.025 (near the low end
# of the st0 range), where the lower-boundary density is large and rapidly
# changing, requiring precision >= 4 to match WienR.
ex5 <- list(a = 0.5, v = 0.5, t0 = 0.275, z = 0.15, sv = 0, sz = 0, st0 = 0.1)
# Lower-boundary density at rt = 0.3
d_wien <- WienerPDF(t = 0.3, response = "lower", a = ex5$a, v = ex5$v,
                    w = ex5$z / ex5$a, t0 = ex5$t0, st0 = ex5$st0)$value
for (p in c(3, 4, 5, 6, 7, 8)) {
  d_rtd <- ddiffusion(0.3, "lower", a = ex5$a, v = ex5$v, t0 = ex5$t0,
                      z = ex5$z, st0 = ex5$st0, precision = p)
  cat(sprintf("    precision=%d: rtdists=%.10f  diff=%.6e\n", p, d_rtd, abs(d_rtd - d_wien)))
}
cat(sprintf("    WienR (target):           %.10f\n\n", d_wien))

# Summary
# - Default precision = 3 is sufficient for densities without variability
#   or when the starting point is near the center.
# - precision = 5 handles most cases with a single variability parameter.
# - precision = 8 may be needed when sz and st0 are both non-zero AND
#   the starting point range extends close to a boundary.
# - The lower boundary is more sensitive than the upper boundary when
#   the starting point is near the lower threshold, and vice versa.
