## ---- fig.height=4, fig.width=7, message=FALSE, warning=FALSE------------
require(rtdists)
require(dplyr)   # for data manipulations and looping
require(tidyr)   # for data manipulations
require(purrr)   # for data manipulations
require(lattice) # for plotting and corresponding themes
require(latticeExtra)
lattice.options(default.theme = standard.theme(color = FALSE))
lattice.options(default.args = list(as.table = TRUE))
options(digits = 3) # only three decimal digits
require(binom)  # for binomial confidence intervals

data(rr98)
rr98 <- rr98[!rr98$outlier,]  #remove outliers

# aggregate data for first plot:
agg_rr98 <- rr98  %>% group_by(id, instruction, strength) %>% 
  summarise(prop = mean(response == "dark"), mean_rt = mean(rt), median_rt = mean(rt)) %>% 
  ungroup()

xyplot(prop ~ strength|id, agg_rr98, group = instruction, type = "b", 
       auto.key = list(lines = TRUE), ylab = "Proportion of 'dark' responses")



## ---- fig.height=6, fig.width=7------------------------------------------

quantiles <- c(0.1, 0.3, 0.5, 0.7, 0.9)
## aggregate data for quantile plot
quantiles_rr98 <- rr98  %>% 
  group_by(id, instruction, strength) %>% 
  nest() %>% 
  mutate(quantiles = map(data, ~ as.data.frame(t(quantile(.x$rt, probs = quantiles))))) %>% 
  unnest(quantiles) %>% 
  gather("quantile", "rt",`10%`:`90%`) %>% 
  arrange(id, instruction, strength)

xyplot(rt ~ strength|id + instruction, quantiles_rr98, group = quantile, type = "b", 
       auto.key = list(lines = TRUE), ylab = "RT (in seconds)", subset = instruction == "speed")

xyplot(rt ~ strength|id + instruction, quantiles_rr98, group = quantile, type = "b", 
       auto.key = FALSE, ylab = "RT (in seconds)", subset = instruction == "accuracy")



## ---- fig.height=4, fig.width=7------------------------------------------

#bins <- c(-0.5, 5.5, 10.5, 13.5, 16.5, 19.5, 25.5, 32.5) # seven bins like RR98
bins <- c(-0.5, 10.5, 13.5, 16.5, 19.5, 32.5)
rr98$strength_bin <- cut(rr98$strength, breaks = bins, include.lowest = TRUE)
levels(rr98$strength_bin) <- as.character(1:7)

# aggregate data for response probability plot:
agg_rr98_bin <- rr98 %>% 
  group_by(id, instruction, strength_bin) %>%
  summarise(n = n(), 
            dark = sum(response == "dark"),
            light = sum(response == "light")) %>%
  ungroup() %>%
  mutate(prop = map2(dark, n, ~ binom.confint(.x, .y, methods = "agresti-coull"))) %>% 
  unnest(prop)
  

knitr::kable(
  rr98 %>% group_by(id, instruction, strength_bin, response) %>%
    summarise(n = n()) %>%
    spread(strength_bin, n)
)



## ---- fig.height=4, fig.width=7------------------------------------------
xyplot(mean ~ strength_bin|id, agg_rr98_bin, group = instruction, type = "b", 
       auto.key = list(lines = TRUE), ylab = "Proportion of 'dark' responses")


## ---- fig.height=6, fig.width=7------------------------------------------

## aggregate data for quantile plot
quantiles_rr98_bin <- rr98  %>% 
  group_by(id, instruction, strength_bin) %>% 
  nest() %>% 
  mutate(quantiles = map(data, ~ as.data.frame(t(quantile(.x$rt, probs = quantiles))))) %>% 
  unnest(quantiles) %>% 
  gather("quantile", "rt",`10%`:`90%`) %>% 
  arrange(id, instruction, strength_bin)

xyplot(rt ~ strength_bin|id + instruction, quantiles_rr98_bin, group = quantile, type = "b", 
       auto.key = list(lines = TRUE), ylab = "RT (in seconds)", subset = instruction == "speed")

xyplot(rt ~ strength_bin|id + instruction, quantiles_rr98_bin, group = quantile, type = "b", 
       auto.key = FALSE, ylab = "RT (in seconds)", subset = instruction == "accuracy")


## ---- fig.height=6, fig.width=7------------------------------------------

agg2_rr98_response <- rr98  %>% 
  group_by(id, instruction, strength_bin, response) %>% 
  nest() %>% 
  mutate(quantiles = map(data, ~ as.data.frame(t(quantile(.x$rt, probs = quantiles))))) %>% 
  unnest(quantiles) %>% 
  gather("quantile", "rt",`10%`:`90%`) %>% 
  arrange(id, instruction, response, strength_bin)

p1 <- xyplot(rt ~ strength_bin|id, agg2_rr98_response, group = quantile, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "speed" & response == "dark", layout = c(3,1))
p2 <- xyplot(rt ~ strength_bin|id, agg2_rr98_response, group = quantile, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "speed" & response == "light", col = "grey")
p1 + as.layer(p2)


p1 <- xyplot(rt ~ strength_bin|id, agg2_rr98_response, group = quantile, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "accuracy" & response == "dark", layout = c(3,1))
p2 <- xyplot(rt ~ strength_bin|id, agg2_rr98_response, group = quantile, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "accuracy" & response == "light", col = "grey")
p1 + as.layer(p2)



## ------------------------------------------------------------------------
d_nested <- rr98 %>% 
  group_by(id, instruction) %>% # we loop across both, id and instruction
  nest()
d_nested

## ------------------------------------------------------------------------
# objective function for diffusion with 1 a. loops over drift to assign drift rates to strength
objective_diffusion_separate <- function(pars, rt, response, drift, ...) {
  non_v_pars <- grep("^v", names(pars), invert = TRUE, value = TRUE)
  base_par <- length(non_v_pars)  # number of non-drift parameters
  densities <- vector("numeric", length(rt))
  for (i in seq_along(levels(drift))) {
    densities[drift == levels(drift)[i]] <- 
      ddiffusion(rt[drift == levels(drift)[i]], response=response[drift == levels(drift)[i]], 
                 a=pars["a"], t0=pars["t0"],  
                 sv=pars["sv"],
                 sz=if ("sz" %in% non_v_pars) pars["sz"] else 0.1,
                 z=if ("z" %in% non_v_pars) pars["z"]*pars["a"] else 0.5*pars["a"],
                 st0=if ("st0" %in% non_v_pars) pars["st0"] else 0, 
                 v=pars[base_par+i])
  }
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}


## ------------------------------------------------------------------------

# function that creates random start values, also 
get_start <- function(base_par, n_drift = 5) {
  start1 <- c(
    a = runif(1, 0.5, 3),
    a_1 = runif(1, 0.5, 3), 
    a_2 = runif(1, 0.5, 3),
    t0 = runif(1, 0, 0.5), 
    z = runif(1, 0.4, 0.6), 
    sz = runif(1, 0, 0.5),
    sv = runif(1, 0, 0.5),
    st0 = runif(1, 0, 0.5),
    d = rnorm(1, 0, 0.05)
  )
  start2 <- sort(rnorm(n_drift), decreasing = FALSE)
  names(start2) <- paste0("v_", seq_len(n_drift))
  c(start1[base_par], start2)
}

# function that tries different random start values until it works:
ensure_fit <- 
  function(data, start_function, objective_function, 
           base_pars, n_drift = 5, n_fits = 1, 
           lower = c(rep(0, length(base_pars)), -Inf,
                     rep(-Inf,length(start_function(base_pars))-length(base_pars)))) {
    best_fit <- list(objective = 1e+06)
  for (i in seq_len(n_fits)) {
    start_ll <- 1e+06
    #browser()
    while(start_ll == 1e+06) {
      start <- start_function(base_pars, n_drift=n_drift)
      start_ll <- objective_function(start, 
                                     rt = data$rt, response = data$response_num, 
                                     drift = factor(data$strength_bin, seq_len(n_drift)), 
                                     instruction = data$instruction)
    }
    cat("\nstart fitting.\n") # just for information to see if it is stuck
    
    fit <- nlminb(start, objective_function, 
                  rt = data$rt, response = data$response_num, 
                  drift = factor(data$strength_bin, seq_len(n_drift)), 
                  instruction = data$instruction,
                  lower = lower)
    
    if (fit$objective < best_fit$objective) best_fit <- fit
  }
  out <- as.data.frame(t(unlist(best_fit[1:3])))
  colnames(out) <- sub("par.", "", colnames(out))
  out
}


## ---- echo=FALSE---------------------------------------------------------
load("rr98_full-diffusion_fits.rda")
load("rr98_full-lba_fits.rda")


## ---- eval = FALSE-------------------------------------------------------
#  fit_diffusion <- d_nested %>%
#    mutate(fit =
#             map(data,
#                 ~ensure_fit(data = ., start_function = get_start,
#                              objective_function = objective_diffusion_separate,
#                              base_pars = c("a", "t0", "sv", "sz", "z")))) %>%
#    unnest(fit)

## ---- eval = FALSE-------------------------------------------------------
#  require(parallel)
#  
#  fit_diffusion <- d_nested
#  fit_diffusion$fit <-
#    mclapply(d_nested$data, function(x)
#      ensure_fit(data = x, start_function = get_start,
#                 objective_function = objective_diffusion_separate,
#                 base_pars = c("a", "t0", "sv", "sz", "z")),
#      mc.cores = 2)
#  fit_diffusion <- unnest(fit_diffusion, fit)

## ------------------------------------------------------------------------

fit_diffusion$data <- NULL
if (!("st0" %in% colnames(fit_diffusion))) fit_diffusion$st0 <- 0
if (!("z" %in% colnames(fit_diffusion))) fit_diffusion$z <- 0.5
if (!("sz" %in% colnames(fit_diffusion))) fit_diffusion$sz <- 0.1
knitr::kable(fit_diffusion)


## ----obtain_fits_not_run, eval = FALSE, include = FALSE------------------
#  
#  require(parallel)
#  
#  fit_diffusion <- d_nested %>%
#    mutate(fit =
#             map(data,
#                 ~ensure_fit(data = ., start_function = get_start,
#                              objective_function = objective_diffusion_separate,
#                              base_pars = c("a", "t0", "sv", "sz", "z")))) %>%
#    unnest(fit)
#  fit_diffusion$data <- NULL
#  
#  fit_diffusion2 <- d_nested
#  fit_diffusion2$fit <-
#    mclapply(d_nested$data, function(x)
#      ensure_fit(data = x, start_function = get_start,
#                 objective_function = objective_diffusion_separate,
#                 base_pars = c("a", "t0", "sv", "sz", "z")),
#      mc.cores = 3)
#  fit_diffusion2 <- unnest(fit_diffusion2, fit)
#  fit_diffusion2$data <- NULL
#  
#  all.equal(as.data.frame(fit_diffusion), as.data.frame(fit_diffusion2), tolerance = 0.01)
#  
#  save(fit_diffusion, fit_diffusion2, file = "rr98_full-diffusion_fits.rda")
#  
#  

## ---- fig.height=5, fig.width=7, message=FALSE---------------------------


# get predicted response proportions
pars_separate_l <- fit_diffusion %>% gather("strength_bin", "v", starts_with("v"))
pars_separate_l$strength_bin <- factor(substr(pars_separate_l$strength_bin, 3,3), 
                                       levels = as.character(seq_len(length(bins)-1)))
#pars_separate_l <- inner_join(pars_separate_l, agg_rr98_bin)
pars_separate_l <- pars_separate_l  %>% group_by(id, instruction, strength_bin) %>%
  mutate(resp_prop = pdiffusion(rt=Inf, response="lower", 
                                a=a, v=v, t0=t0, sz = sz, z=a*z, sv=sv, st0=st0))

p1 <- xyplot(mean ~ strength_bin|id + instruction, agg_rr98_bin, type = "b", auto.key = 
               list(lines = TRUE), ylab = "Proportion of 'dark' responses", col = "grey")
p2 <- segplot(strength_bin ~ upper+lower|id + instruction, agg_rr98_bin, 
              auto.key = list(lines = TRUE), ylab = "Proportion of 'dark' responses", 
              col = "grey", horizontal = FALSE, segments.fun = panel.arrows,  
              draw.bands = FALSE, angle = 90, length = 0.05, ends = "both")
p3 <- xyplot(resp_prop ~ strength_bin|id + instruction, pars_separate_l, type = "b", 
             auto.key = list(lines = TRUE), ylab = "Proportion of 'dark' responses", 
             col = "black")
p2 + as.layer(p1) + as.layer(p3)


## ---- fig.height=6, fig.width=7, message=FALSE---------------------------

# get predicted quantiles (uses predicted response proportions)
separate_pred_dark <- pars_separate_l %>% do(as.data.frame(t(
  qdiffusion(quantiles*.$resp_prop, response="lower", 
             a=.$a, v=.$v, t0=.$t0, sz = .$sz, z = .$z*.$a, sv=.$sv, st0=.$st0)))) %>% 
  ungroup() %>% gather("quantiles", "dark", V1:V5)
separate_pred_light <- pars_separate_l %>% do(as.data.frame(t(
  qdiffusion(quantiles*(1-.$resp_prop), response="upper", 
             a=.$a, v=.$v, t0=.$t0, sz = .$sz, z = .$z*.$a, sv=.$sv, st0=.$st0)))) %>% 
  ungroup() %>% gather("quantiles", "light", V1:V5)

#separate_pred_light %>% filter(is.na(light))
separate_pred <- inner_join(separate_pred_dark, separate_pred_light)
separate_pred$quantiles <- factor(separate_pred$quantiles, 
                                  levels = c("V5", "V4", "V3", "V2", "V1"), 
                                  labels = c("90%", "70%", "50%", "30%", "10%"))
separate_pred <- separate_pred %>% gather("response", "rt", dark, light)

# get SE for observed quantiles
agg2_rr98_response_se <- rr98  %>% group_by(id, instruction, strength_bin, response) %>% 
  summarise(se_median = sqrt(pi/2)*(sd(rt)/sqrt(n()))) %>%
  ungroup()

# calculate error bars for quantiles.
agg2_rr98_response <- left_join(agg2_rr98_response, agg2_rr98_response_se)
agg2_rr98_response <- agg2_rr98_response %>%
  mutate(lower = rt-se_median, upper = rt+se_median)


p1 <- xyplot(rt ~ strength_bin|id+response, agg2_rr98_response, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "speed" & quantile == "50%", 
             layout = c(3,2), col = "grey")
p1e <- segplot(strength_bin ~ upper+lower|id+response, agg2_rr98_response, 
               auto.key = list(lines = TRUE), ylab = "Proportion of 'dark' responses", 
               col = "grey", horizontal = FALSE, segments.fun = panel.arrows,  
               draw.bands = FALSE, angle = 90, length = 0.05, ends = "both", 
               subset = instruction == "speed" & quantile == "50%", layout = c(3,2))
p2 <- xyplot(rt ~ strength_bin|id + response, separate_pred, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "speed" & quantiles == "50%", 
             scales = list(y = list(limits = c(0.25, 0.5))))
p2 + as.layer(p1) + as.layer(p1e)


## ---- fig.height=6, fig.width=7------------------------------------------

p1 <- xyplot(rt ~ strength_bin|id+response, agg2_rr98_response, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "accuracy" & quantile == "50%", 
             layout = c(3,2), col = "grey")
p1e <- segplot(strength_bin ~ upper+lower|id+response, agg2_rr98_response, 
               auto.key = list(lines = TRUE), ylab = "Proportion of 'dark' responses", 
               col = "grey", horizontal = FALSE, segments.fun = panel.arrows,  
               draw.bands = FALSE, angle = 90, length = 0.05, ends = "both", 
               subset = instruction == "accuracy" & quantile == "50%", layout = c(3,2))
p2 <- xyplot(rt ~ strength_bin|id + response, separate_pred, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "accuracy" & quantiles == "50%", 
             scales = list(y = list(limits = c(0.2, 1.5))))
p2 + as.layer(p1) + as.layer(p1e)


## ---- fig.height=7, fig.width=7------------------------------------------

p1 <- xyplot(rt ~ strength_bin|id+response, agg2_rr98_response, group = quantile, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "speed", layout = c(3,2), col = "grey")
p1e <- segplot(strength_bin ~ upper+lower|id+response, agg2_rr98_response, 
               auto.key = list(lines = TRUE), ylab = "Proportion of 'dark' responses", 
               col = "grey", horizontal = FALSE, segments.fun = panel.arrows,  
               draw.bands = FALSE, angle = 90, length = 0.05, ends = "both", 
               subset = instruction == "speed")
p2 <- xyplot(rt ~ strength_bin|id + response, separate_pred, group = quantiles, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "speed", scales = list(y = list(limits = c(0.2, 0.9))))
p2 + as.layer(p1) + as.layer(p1e)


## ---- fig.height=7, fig.width=7------------------------------------------

p1 <- xyplot(rt ~ strength_bin|id+response, agg2_rr98_response, group = quantile, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "accuracy", layout = c(3,2), col = "grey")
p1e <- segplot(strength_bin ~ upper+lower|id+response, agg2_rr98_response, 
               auto.key = list(lines = TRUE), ylab = "Proportion of 'dark' responses", 
               col = "grey", horizontal = FALSE, segments.fun = panel.arrows,  
               draw.bands = FALSE, angle = 90, length = 0.05, ends = "both", 
               subset = instruction == "accuracy")
p2 <- xyplot(rt ~ strength_bin|id + response, separate_pred, group = quantiles, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "accuracy", scales = list(y = list(limits = c(0.1, 3.0))))
p2 + as.layer(p1) + as.layer(p1e)


## ------------------------------------------------------------------------

# objective function for diffusion with 1 a. loops over drift to assign drift rates to strength
objective_lba_separate <- function(pars, rt, response, drift, ...) {
  non_v_pars <- grep("^v", names(pars), invert = TRUE, value = TRUE)
  base_par <- length(non_v_pars)  # number of non-drift parameters
  densities <- vector("numeric", length(rt))
  for (i in seq_along(levels(drift))) {
    if (sum(drift == levels(drift)[i]) == 0) next
    densities[drift == levels(drift)[i]] <- dLBA(
      rt[drift == levels(drift)[i]], 
      response=response[drift == levels(drift)[i]],
      A = list(pars["a_1"], pars["a_2"]), 
      b = max(pars["a_1"], pars["a_2"])+pars["b"], 
      t0 = pars["t0"], 
      mean_v = c(pars[i], 1-pars[i]), 
      sd_v = pars["sv"], silent=TRUE)
  }
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}

# function that creates random start values
get_start_lba <- function(base_par, n_drift = 10) {
  start1 <- c(
    a = runif(1, 0.5, 3),
    a_1 = runif(1, 0.5, 3), 
    a_2 = runif(1, 0.5, 3),
    t0 = runif(1, 0, 0.5), 
    b = runif(1, 0, 0.5), 
    sv = runif(1, 0.5, 1.5),
    st0 = runif(1, 0, 0.5)
  )
  start2 <- sort(rnorm(n_drift), decreasing = FALSE)
  names(start2) <- paste0("v_", seq_len(n_drift))
  c(start2, start1[base_par])
}


## ---- eval=FALSE---------------------------------------------------------
#  
#  fit_lba <- d_nested %>%
#    mutate(fit =
#             map(data,
#                 ~ensure_fit(data = ., start_function = get_start_lba,
#                        objective_function = objective_lba_separate,
#                        base_pars = c("a_1", "a_2", "t0", "b", "sv"),
#                        lower = c(rep(-Inf, 5), rep(0, 5)),
#                        n_drift = 5, n_fits = 10))) %>%
#    unnest(fit)
#  

## ---- eval = FALSE-------------------------------------------------------
#  require(parallel)
#  
#  fit_lba <- d_nested
#  fit_lba$fit <-
#    mclapply(d_nested$data, function(x)
#      ensure_fit(data = x, start_function = get_start_lba,
#                        objective_function = objective_lba_separate,
#                        base_pars = c("a_1", "a_2", "t0", "b", "sv"),
#                        lower = c(rep(-Inf, 5), rep(0, 5)),
#                        n_drift = 5, n_fits = 10),
#      mc.cores = 2)
#  fit_lba <- unnest(fit_lba, fit)

## ------------------------------------------------------------------------
knitr::kable(fit_lba)


## ----obtain_fits_lba, eval = FALSE, include = FALSE----------------------
#  fit_lba <- d_nested %>%
#    mutate(fit =
#             map(data,
#                 ~ensure_fit(data = ., start_function = get_start_lba,
#                        objective_function = objective_lba_separate,
#                        base_pars = c("a_1", "a_2", "t0", "b", "sv"),
#                        lower = c(rep(-Inf, 5), rep(0, 5)),
#                        n_drift = 5, n_fits = 10))) %>%
#    unnest(fit)
#  fit_lba$data <- NULL
#  
#  fit_lba2 <- d_nested
#  fit_lba2$fit <-
#    mclapply(d_nested$data, function(x)
#      ensure_fit(data = x, start_function = get_start_lba,
#                        objective_function = objective_lba_separate,
#                        base_pars = c("a_1", "a_2", "t0", "b", "sv"),
#                        lower = c(rep(-Inf, 5), rep(0, 5)),
#                        n_drift = 5, n_fits = 10),
#      mc.cores = 2)
#  fit_lba2 <- unnest(fit_lba2, fit)
#  fit_lba2$data <- NULL
#  
#  all.equal(as.data.frame(fit_lba), as.data.frame(fit_lba2), tolerance = 0.03)
#  save(fit_lba, fit_lba2, file = "rr98_full-lba_fits.rda")
#  
#  
#  
#  # objective function for LBA with 1 a. loops over drift to assign drift rates to strength
#  objective_lba_separate <- function(pars, rt, response, drift, ...) {
#    non_v_pars <- grep("^v", names(pars), invert = TRUE, value = TRUE)
#    base_par <- length(non_v_pars)  # number of non-drift parameters
#    densities <- vector("numeric", length(rt))
#    for (i in seq_along(levels(drift))) {
#      if (sum(drift == levels(drift)[i]) == 0) next
#      densities[drift == levels(drift)[i]] <- dLBA(
#        rt[drift == levels(drift)[i]],
#        response=response[drift == levels(drift)[i]],
#        A = pars["a"], b = pars["a"]+pars["b"],
#        t0 = pars["t0"],
#        mean_v = pars[((i-1)*2+1):((i-1)*2+2)],
#        sd_v = c(1, pars["sv"]), silent=TRUE)
#    }
#    if (any(densities == 0)) return(1e6)
#    return(-sum(log(densities)))
#  }
#  
#  # objective function for diffusion with 1 a. loops over drift to assign drift rates to strength
#  objective_lba_separate <- function(pars, rt, response, drift, ...) {
#    non_v_pars <- grep("^v", names(pars), invert = TRUE, value = TRUE)
#    base_par <- length(non_v_pars)  # number of non-drift parameters
#    densities <- vector("numeric", length(rt))
#    for (i in seq_along(levels(drift))) {
#      if (sum(drift == levels(drift)[i]) == 0) next
#      densities[drift == levels(drift)[i]] <- dLBA(
#        rt[drift == levels(drift)[i]],
#        response=response[drift == levels(drift)[i]],
#        A = list(pars["a_1"], pars["a_2"]), b = list(pars["a_1"]+pars["b"], pars["a_2"]+pars["b"]),
#        t0 = pars["t0"],
#        mean_v = c(pars[i], 1-pars[i]),
#        sd_v = pars["sv"], silent=TRUE)
#    }
#    if (any(densities == 0)) return(1e6)
#    return(-sum(log(densities)))
#  }
#  

## ---- fig.height=5, fig.width=7, message=FALSE---------------------------
# get predicted response proportions
lba_pars_separate_l <- fit_lba %>% gather("strength_bin", "v", starts_with("v"))
lba_pars_separate_l$strength_bin <- factor(substr(lba_pars_separate_l$strength_bin, 3,3), 
                                       levels = as.character(seq_len(length(bins)-1)))
#pars_separate_l <- inner_join(pars_separate_l, agg_rr98_bin)
lba_pars_separate_l <- lba_pars_separate_l  %>% group_by(id, instruction, strength_bin) %>%
  mutate(resp_prop = pLBA(rt=Inf, response=1, A=list(a_1, a_2), sd_v=sv,
                          mean_v=c(v, 1-v), t0=t0, b=max(a_1, a_2)+b, silent=TRUE))

p1 <- xyplot(mean ~ strength_bin|id + instruction, agg_rr98_bin, type = "b", auto.key = 
               list(lines = TRUE), ylab = "Proportion of 'dark' responses", col = "grey")
p2 <- segplot(strength_bin ~ upper+lower|id + instruction, agg_rr98_bin, 
              auto.key = list(lines = TRUE), ylab = "Proportion of 'dark' responses", 
              col = "grey", horizontal = FALSE, segments.fun = panel.arrows,  
              draw.bands = FALSE, angle = 90, length = 0.05, ends = "both")
p3 <- xyplot(resp_prop ~ strength_bin|id + instruction, lba_pars_separate_l, type = "b", 
             auto.key = list(lines = TRUE), ylab = "Proportion of 'dark' responses", 
             col = "black")
p2 + as.layer(p1) + as.layer(p3)


## ---- fig.height=6, fig.width=7, message=FALSE---------------------------

# get predicted quantiles (uses predicted response proportions)
lba_separate_pred_dark <- lba_pars_separate_l %>% do(as.data.frame(t(
  qLBA(quantiles*.$resp_prop, response=1, A=list(.$a_1, .$a_2), sd_v=.$sv,
       mean_v=c(.$v, 1-.$v), t0=.$t0, b=max(.$a_1, .$a_2)+.$b, silent=TRUE)))) %>% 
  ungroup() %>% gather("quantiles", "dark", V1:V5)
lba_separate_pred_light <- lba_pars_separate_l %>% do(as.data.frame(t(
  qLBA(quantiles*(1-.$resp_prop), response=2, A=list(.$a_1, .$a_2), sd_v=.$sv,
       mean_v=c(.$v, 1-.$v), t0=.$t0, b=max(.$a_1, .$a_2)+.$b, silent=TRUE)))) %>% 
  ungroup() %>% gather("quantiles", "light", V1:V5)

#separate_pred_light %>% filter(is.na(light))
lba_separate_pred <- inner_join(lba_separate_pred_dark, lba_separate_pred_light)
lba_separate_pred$quantiles <- factor(lba_separate_pred$quantiles, 
                                  levels = c("V5", "V4", "V3", "V2", "V1"), 
                                  labels = c("90%", "70%", "50%", "30%", "10%"))
lba_separate_pred <- lba_separate_pred %>% gather("response", "rt", dark, light)

p1 <- xyplot(rt ~ strength_bin|id+response, agg2_rr98_response, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "speed" & quantile == "50%", 
             layout = c(3,2), col = "grey")
p1e <- segplot(strength_bin ~ upper+lower|id+response, agg2_rr98_response, 
               auto.key = list(lines = TRUE), ylab = "Proportion of 'dark' responses", 
               col = "grey", horizontal = FALSE, segments.fun = panel.arrows,  
               draw.bands = FALSE, angle = 90, length = 0.05, ends = "both", 
               subset = instruction == "speed" & quantile == "50%", layout = c(3,2))
p2 <- xyplot(rt ~ strength_bin|id + response, lba_separate_pred, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "speed" & quantiles == "50%", 
             scales = list(y = list(limits = c(0.25, 0.5))))
p2 + as.layer(p1) + as.layer(p1e)


## ---- fig.height=6, fig.width=7------------------------------------------

p1 <- xyplot(rt ~ strength_bin|id+response, agg2_rr98_response, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "accuracy" & quantile == "50%", 
             layout = c(3,2), col = "grey")
p1e <- segplot(strength_bin ~ upper+lower|id+response, agg2_rr98_response, 
               auto.key = list(lines = TRUE), ylab = "Proportion of 'dark' responses", 
               col = "grey", horizontal = FALSE, segments.fun = panel.arrows,  
               draw.bands = FALSE, angle = 90, length = 0.05, ends = "both", 
               subset = instruction == "accuracy" & quantile == "50%", layout = c(3,2))
p2 <- xyplot(rt ~ strength_bin|id + response, lba_separate_pred, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "accuracy" & quantiles == "50%", 
             scales = list(y = list(limits = c(0.2, 1.5))))
p2 + as.layer(p1) + as.layer(p1e)


## ---- fig.height=7, fig.width=7------------------------------------------

p1 <- xyplot(rt ~ strength_bin|id+response, agg2_rr98_response, group = quantile, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "speed", layout = c(3,2), col = "grey")
p1e <- segplot(strength_bin ~ upper+lower|id+response, agg2_rr98_response, 
               auto.key = list(lines = TRUE), ylab = "Proportion of 'dark' responses", 
               col = "grey", horizontal = FALSE, segments.fun = panel.arrows,  
               draw.bands = FALSE, angle = 90, length = 0.05, ends = "both", 
               subset = instruction == "speed")
p2 <- xyplot(rt ~ strength_bin|id + response, lba_separate_pred, group = quantiles, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "speed", scales = list(y = list(limits = c(0.2, 0.6))))
p2 + as.layer(p1) + as.layer(p1e)


## ---- fig.height=7, fig.width=7------------------------------------------

p1 <- xyplot(rt ~ strength_bin|id+response, agg2_rr98_response, group = quantile, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "accuracy", layout = c(3,2), col = "grey")
p1e <- segplot(strength_bin ~ upper+lower|id+response, agg2_rr98_response, 
               auto.key = list(lines = TRUE), ylab = "Proportion of 'dark' responses", 
               col = "grey", horizontal = FALSE, segments.fun = panel.arrows,  
               draw.bands = FALSE, angle = 90, length = 0.05, ends = "both", 
               subset = instruction == "accuracy")
p2 <- xyplot(rt ~ strength_bin|id + response, lba_separate_pred, group = quantiles, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "accuracy", scales = list(y = list(limits = c(0.1, 3.3))))
p2 + as.layer(p1) + as.layer(p1e)


## ---- fig.height=6.5, fig.width=7, message=FALSE-------------------------

key <- simpleKey(text = c("data", "LBA", "Diffusion"), lines = TRUE)
key$lines$col <- c("grey", "black", "black")
key$lines$lty <- c(1, 1, 2)
key$points$col <- c("grey", "black", "black")
key$points$pch <- c(1, 0, 4)

p1 <- xyplot(mean ~ strength_bin|id + instruction, agg_rr98_bin, type = "b", auto.key = 
               list(lines = TRUE), ylab = "Proportion of 'dark' responses", col = "grey")
p2 <- segplot(strength_bin ~ upper+lower|id + instruction, agg_rr98_bin, 
              auto.key = list(lines = TRUE), ylab = "Proportion of 'dark' responses", 
              col = "grey", horizontal = FALSE, segments.fun = panel.arrows,  
              draw.bands = FALSE, angle = 90, length = 0.05, ends = "both")
p3 <- xyplot(resp_prop ~ strength_bin|id + instruction, lba_pars_separate_l, type = "b", 
             auto.key = list(lines = TRUE), ylab = "Proportion of 'dark' responses", 
             col = "black", pch = 0)
p4 <- xyplot(resp_prop ~ strength_bin|id + instruction, pars_separate_l, type = "b", 
             auto.key = list(lines = TRUE), ylab = "Proportion of 'dark' responses", 
             col = "black", lty = 2, key = key, pch=4)
p4 + as.layer(p2) + as.layer(p1) + as.layer(p3)


## ---- fig.height=6.5, fig.width=7, message=FALSE-------------------------

p1 <- xyplot(rt ~ strength_bin|id+response, agg2_rr98_response, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "speed" & quantile == "50%", 
             layout = c(3,2), col = "grey")
p1e <- segplot(strength_bin ~ upper+lower|id+response, agg2_rr98_response, 
               auto.key = list(lines = TRUE), ylab = "Proportion of 'dark' responses", 
               col = "grey", horizontal = FALSE, segments.fun = panel.arrows,  
               draw.bands = FALSE, angle = 90, length = 0.05, ends = "both", 
               subset = instruction == "speed" & quantile == "50%", layout = c(3,2))
p2 <- xyplot(rt ~ strength_bin|id + response, lba_separate_pred, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "speed" & quantiles == "50%", 
             scales = list(y = list(limits = c(0.25, 0.5))), pch = 0)
p3 <- xyplot(rt ~ strength_bin|id + response, separate_pred, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "speed" & quantiles == "50%", 
             scales = list(y = list(limits = c(0.25, 0.5))),
             col = "black", lty = 2, key = key, pch=4)

p3 + as.layer(p2) + as.layer(p1) + as.layer(p1e)


## ---- fig.height=6.5, fig.width=7, message=FALSE-------------------------

p1 <- xyplot(rt ~ strength_bin|id+response, agg2_rr98_response, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "accuracy" & quantile == "50%", 
             layout = c(3,2), col = "grey")
p1e <- segplot(strength_bin ~ upper+lower|id+response, agg2_rr98_response, 
               auto.key = list(lines = TRUE), ylab = "Proportion of 'dark' responses", 
               col = "grey", horizontal = FALSE, segments.fun = panel.arrows,  
               draw.bands = FALSE, angle = 90, length = 0.05, ends = "both", 
               subset = instruction == "accuracy" & quantile == "50%", layout = c(3,2))
p2 <- xyplot(rt ~ strength_bin|id + response, lba_separate_pred, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "accuracy" & quantiles == "50%", 
             pch = 0)
p3 <- xyplot(rt ~ strength_bin|id + response, separate_pred, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "accuracy" & quantiles == "50%", 
             scales = list(y = list(limits = c(0.2, 1.5))),
             col = "black", lty = 2, key = key, pch=4)

p3 + as.layer(p2) + as.layer(p1) + as.layer(p1e)


