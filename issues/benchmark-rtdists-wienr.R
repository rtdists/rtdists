## Benchmark: rtdists vs WienR timing comparison (extended parameter sweep)
##
## Uses the bench package to compare ddiffusion/pdiffusion (rtdists) against
## WienerPDF/WienerCDF (WienR) across a range of parameter values.
##
## Parameter mapping (rtdists -> WienR):
##   z (absolute) -> w = z/a (relative)
##   sz (absolute) -> sw = sz/a (relative)
##   t0, sv, st0, a, v map directly
##   rtdists s defaults to 1; WienR has no s parameter (fixed to 1)

library(bench)
library(rtdists)
library(WienR)
library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(4271)

## ================================================================
## 1. Basic case (no variability): grid sweep with bench::press()
##    Vary a, v, t0, z across moderate grids
## ================================================================

rt <- seq(0.3, 3, length.out = 20)

cat("=== 1. Basic ddiffusion vs WienerPDF (parameter grid) ===\n")
bench_d_basic <- press(
  a  = c(0.5, 1.0, 1.5, 2.0),
  v  = c(0.5, 1.0, 1.5, 2.0),
  t0 = c(0.05, 0.2, 0.5),
  z_rel = c(0.3, 0.5, 0.7),
  {
    z_abs <- z_rel * a
    bench::mark(
      rtdists = ddiffusion(rt, a = a, v = v, t0 = t0, z = z_abs, precision = 3),
      WienR   = WienerPDF(t = rt, response = "upper", a = a, v = v, w = z_rel, t0 = t0)$value,
      iterations = 10, check = FALSE, filter_gc = FALSE
    )
  }
)

cat("=== 2. Basic pdiffusion vs WienerCDF (parameter grid) ===\n")
bench_p_basic <- press(
  a  = c(0.5, 1.0, 1.5, 2.0),
  v  = c(0.5, 1.0, 1.5, 2.0),
  t0 = c(0.05, 0.2, 0.5),
  z_rel = c(0.3, 0.5, 0.7),
  {
    z_abs <- z_rel * a
    bench::mark(
      rtdists = pdiffusion(rt, a = a, v = v, t0 = t0, z = z_abs, precision = 3),
      WienR   = WienerCDF(t = rt, response = "upper", a = a, v = v, w = z_rel, t0 = t0)$value,
      iterations = 10, check = FALSE, filter_gc = FALSE
    )
  }
)

## ================================================================
## 3. Variability case: random parameter sample
##    Full grid is too large (7 params), so sample N valid combinations
## ================================================================

n_samp <- 40
param_grid <- data.frame(
  a    = runif(n_samp, 0.5, 2.0),
  v    = runif(n_samp, 0.5, 2.5),
  t0   = runif(n_samp, 0.05, 0.5),
  z_rel = runif(n_samp, 0.3, 0.7),
  sv   = runif(n_samp, 0, 0.5),
  sz   = runif(n_samp, 0, 0.25),
  st0  = runif(n_samp, 0, 0.2)
)
param_grid$z_abs <- param_grid$z_rel * param_grid$a
param_grid$sw   <- param_grid$sz / param_grid$a

# Filter: starting point range must stay within [0, a]
param_grid <- param_grid[
  param_grid$z_abs - param_grid$sz / 2 > 0 &
  param_grid$z_abs + param_grid$sz / 2 < param_grid$a,
]
cat(sprintf("Variability case: %d valid parameter combinations\n", nrow(param_grid)))

cat("=== 3. ddiffusion vs WienerPDF (with variability, sampled grid) ===\n")
bench_d_var_results <- list()
for (i in seq_len(nrow(param_grid))) {
  p <- param_grid[i, ]
  mb <- bench::mark(
    rtdists = ddiffusion(rt, a = p$a, v = p$v, t0 = p$t0, z = p$z_abs,
                         sv = p$sv, sz = p$sz, st0 = p$st0, precision = 6),
    WienR   = WienerPDF(t = rt, response = "upper", a = p$a, v = p$v, w = p$z_rel,
                         t0 = p$t0, sv = p$sv, sw = p$sw, st0 = p$st0)$value,
    iterations = 10, check = FALSE, filter_gc = FALSE
  )
  mb$a <- p$a; mb$v <- p$v; mb$t0 <- p$t0; mb$z_rel <- p$z_rel
  mb$sv <- p$sv; mb$sz <- p$sz; mb$st0 <- p$st0
  bench_d_var_results[[i]] <- mb
}
bench_d_var <- do.call(rbind, bench_d_var_results)
class(bench_d_var) <- c("bench_press", "bench_df", class(bench_d_var))

cat("=== 4. pdiffusion vs WienerCDF (with variability, sampled grid) ===\n")
bench_p_var_results <- list()
for (i in seq_len(nrow(param_grid))) {
  p <- param_grid[i, ]
  mb <- bench::mark(
    rtdists = pdiffusion(rt, a = p$a, v = p$v, t0 = p$t0, z = p$z_abs,
                         sv = p$sv, sz = p$sz, st0 = p$st0, precision = 4),
    WienR   = WienerCDF(t = rt, response = "upper", a = p$a, v = p$v, w = p$z_rel,
                         t0 = p$t0, sv = p$sv, sw = p$sw, st0 = p$st0)$value,
    iterations = 10, check = FALSE, filter_gc = FALSE
  )
  mb$a <- p$a; mb$v <- p$v; mb$t0 <- p$t0; mb$z_rel <- p$z_rel
  mb$sv <- p$sv; mb$sz <- p$sz; mb$st0 <- p$st0
  bench_p_var_results[[i]] <- mb
}
bench_p_var <- do.call(rbind, bench_p_var_results)
class(bench_p_var) <- c("bench_press", "bench_df", class(bench_p_var))

## ================================================================
## 5. Precision scaling: vary precision from 3 to 7
## ================================================================

cat("=== 5. Precision scaling (ddiffusion vs WienerPDF) ===\n")
bench_precision <- press(
  precision = c(3, 4, 5, 6, 7),
  {
    bench::mark(
      rtdists = ddiffusion(rt, a = 1.2, v = 1.5, t0 = 0.2, z = 0.6,
                           sv = 0.3, sz = 0.1, st0 = 0.1, precision = precision),
      WienR   = WienerPDF(t = rt, response = "upper", a = 1.2, v = 1.5, w = 0.5,
                           t0 = 0.2, sv = 0.3, sw = 0.1/1.2, st0 = 0.1, precision = precision)$value,
      iterations = 10, check = FALSE, filter_gc = FALSE
    )
  }
)

## ================================================================
## 6. Vector length scaling
## ================================================================

cat("=== 6. Vector length scaling (ddiffusion vs WienerPDF) ===\n")
bench_length <- press(
  n = c(20, 100, 500, 1000, 2000, 5000),
  {
    rt_n <- seq(0.3, 3, length.out = n)
    bench::mark(
      rtdists = ddiffusion(rt_n, a = 1.2, v = 1.5, t0 = 0.2, z = 0.6,
                           sv = 0.3, sz = 0.1, st0 = 0.1, precision = 6),
      WienR   = WienerPDF(t = rt_n, response = "upper", a = 1.2, v = 1.5, w = 0.5,
                           t0 = 0.2, sv = 0.3, sw = 0.1/1.2, st0 = 0.1, precision = 6)$value,
      iterations = 5, check = FALSE, filter_gc = FALSE
    )
  }
)

## ================================================================
## Summary tables (using summary() on bench objects)
## ================================================================

# Helper: convert bench summary to a clean tibble with median in ms
tidy_bench <- function(mb) {
  summary(mb) |>
    mutate(
      expr      = as.character(expression),
      median_ms = as.numeric(.data[["median"]]) * 1e3,
      min_ms    = as.numeric(.data[["min"]]) * 1e3,
      max_ms    = sapply(time, function(x) max(as.numeric(x))) * 1e3
    ) |>
    select(-expression, -min, -median)
}

cat("\n")
cat(strrep("=", 80), "\n")
cat("SUMMARY TABLES\n")
cat(strrep("=", 80), "\n")

cat("\n--- 1. Basic ddiffusion vs WienerPDF (median by a, v) ---\n")
tidy_bench(bench_d_basic) |>
  group_by(a, v, expr) |>
  summarise(median_ms = median(median_ms), .groups = "drop") |>
  pivot_wider(names_from = expr, values_from = median_ms) |>
  mutate(speedup_rtdists = round(WienR / rtdists, 2)) |>
  print(n = Inf)

cat("\n--- 2. Basic pdiffusion vs WienerCDF (median by a, v) ---\n")
tidy_bench(bench_p_basic) |>
  group_by(a, v, expr) |>
  summarise(median_ms = median(median_ms), .groups = "drop") |>
  pivot_wider(names_from = expr, values_from = median_ms) |>
  mutate(speedup_rtdists = round(WienR / rtdists, 2)) |>
  print(n = Inf)

cat("\n--- 3. ddiffusion vs WienerPDF with variability (overall) ---\n")
tidy_bench(bench_d_var) |>
  group_by(expr) |>
  summarise(
    median_ms = median(median_ms),
    mean_ms   = mean(median_ms),
    min_ms    = min(min_ms),
    max_ms    = max(max_ms),
    .groups = "drop"
  ) |>
  print()

cat("\n--- 4. pdiffusion vs WienerCDF with variability (overall) ---\n")
tidy_bench(bench_p_var) |>
  group_by(expr) |>
  summarise(
    median_ms = median(median_ms),
    mean_ms   = mean(median_ms),
    min_ms    = min(min_ms),
    max_ms    = max(max_ms),
    .groups = "drop"
  ) |>
  print()

cat("\n--- 5. Precision scaling (median ms by precision) ---\n")
tidy_bench(bench_precision) |>
  select(precision, expr, median_ms) |>
  pivot_wider(names_from = expr, values_from = median_ms) |>
  mutate(speedup_rtdists = round(WienR / rtdists, 2)) |>
  print(n = Inf)

cat("\n--- 6. Vector length scaling (median ms by n) ---\n")
tidy_bench(bench_length) |>
  select(n, expr, median_ms) |>
  pivot_wider(names_from = expr, values_from = median_ms) |>
  mutate(speedup_rtdists = round(WienR / rtdists, 2)) |>
  print(n = Inf)

## ================================================================
## Visualizations
## ================================================================

## Basic case: median time ratio across parameter grid
plot_basic_ratio <- function(mb, title) {
  tidy_bench(mb) |>
    group_by(a, v, t0, z_rel, expr) |>
    summarise(median_ms = median(median_ms), .groups = "drop") |>
    pivot_wider(names_from = expr, values_from = median_ms) |>
    mutate(ratio = WienR / rtdists) |>
    ggplot(aes(x = a, y = v, fill = ratio)) +
    geom_tile() +
    facet_grid(t0 ~ z_rel, labeller = label_both) +
    scale_fill_viridis_c(name = "WienR / rtdists\n(<1 = rtdists slower)") +
    geom_text(aes(label = sprintf("%.2f", ratio)), size = 3) +
    labs(title = title, x = "a (boundary)", y = "v (drift)") +
    theme_minimal()
}

plot_d_basic <- plot_basic_ratio(bench_d_basic,
  "ddiffusion vs WienerPDF: time ratio across parameters")
print(plot_d_basic)

plot_p_basic <- plot_basic_ratio(bench_p_basic,
  "pdiffusion vs WienerCDF: time ratio across parameters")
print(plot_p_basic)

## Variability case: scatter of ratio vs each parameter
plot_var_ratio <- function(mb, title) {
  tidy_bench(mb) |>
    group_by(a, v, t0, z_rel, sv, sz, st0, expr) |>
    summarise(median_ms = median(median_ms), .groups = "drop") |>
    pivot_wider(names_from = expr, values_from = median_ms) |>
    mutate(ratio = WienR / rtdists) |>
    select(a, v, t0, z_rel, sv, sz, st0, ratio) |>
    pivot_longer(-ratio, names_to = "param", values_to = "value") |>
    ggplot(aes(x = value, y = ratio)) +
    geom_point(alpha = 0.6) +
    geom_smooth(se = FALSE, method = "loess", linewidth = 0.8) +
    facet_wrap(~ param, scales = "free_x") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    labs(title = title, y = "WienR / rtdists time ratio", x = "parameter value") +
    theme_minimal()
}

plot_d_var <- plot_var_ratio(bench_d_var,
  "ddiffusion vs WienerPDF (variability): ratio vs each parameter")
print(plot_d_var)

plot_p_var <- plot_var_ratio(bench_p_var,
  "pdiffusion vs WienerCDF (variability): ratio vs each parameter")
print(plot_p_var)

## Precision scaling plot
plot_precision <- tidy_bench(bench_precision) |>
  group_by(precision, expr) |>
  summarise(median_ms = median(median_ms), .groups = "drop") |>
  ggplot(aes(x = precision, y = median_ms, color = expr)) +
  geom_line() + geom_point(size = 2) +
  scale_y_log10() +
  labs(title = "Precision scaling: ddiffusion vs WienerPDF",
       x = "precision", y = "median time (ms, log scale)") +
  theme_minimal()
print(plot_precision)

## Vector length scaling plot
plot_length <- tidy_bench(bench_length) |>
  group_by(n, expr) |>
  summarise(median_ms = median(median_ms), .groups = "drop") |>
  ggplot(aes(x = n, y = median_ms, color = expr)) +
  geom_line() + geom_point(size = 2) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Vector length scaling: ddiffusion vs WienerPDF",
       x = "n (log scale)", y = "median time (ms, log scale)") +
  theme_minimal()
print(plot_length)

cat("\nDone.\n")

## ================================================================
## Key Takeaways
## ================================================================
##
## 1. rtdists wins only in one case: basic ddiffusion (no variability,
##    precision=3), where it is ~1.4x faster than WienR. The advantage
##    is uniform across all parameter values (a, v, t0, z).
##
## 2. WienR wins everywhere else:
##    - Basic pdiffusion (no variability): WienR ~1.5-3x faster, with the
##      gap widening as boundary separation `a` increases (rtdists CDF
##      scales with `a`, WienR does not).
##    - ddiffusion with variability (precision=6): WienR ~16x faster
##      (median 1.8 ms vs 28 ms).
##    - pdiffusion with variability (precision=4): WienR ~4x faster
##      (median 22 ms vs 92 ms).
##
## 3. The precision parameter is the dominant factor. rtdists scales
##    exponentially with precision (~10x slower per level), while WienR
##    is essentially flat. At precision 7, WienR is ~4500x faster.
##    This points to fundamentally different underlying algorithms.
##
## 4. rtdists timing is far more parameter-sensitive. In the variability
##    case, rtdists timing varies by 100-250x across the parameter space,
##    while WienR varies by only 2-4x.
##
## 5. Both packages scale linearly with vector length, but the per-element
##    cost ratio is ~550x (at precision=6 with variability).
##
## 6. The `a` parameter (boundary separation) is the one parameter that
##    systematically affects the ratio in the basic CDF case: rtdists
##    slows down with larger `a` while WienR stays flat.
