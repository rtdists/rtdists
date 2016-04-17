## ---- fig.height=4, fig.width=7------------------------------------------
require(rtdists)
require(dplyr)   # for data manipulations and looping
require(tidyr)   # for data manipulations
require(lattice) # for plotting and corresponding themes
require(latticeExtra)
lattice.options(default.theme = standard.theme(color = FALSE))
lattice.options(default.args = list(as.table = TRUE))
options(digits = 3) # only three decimal digits
require(binom)  # for binomial confidence intervals

data(rr98)
rr98 <- rr98[!rr98$outlier,]  #remove outliers

#bins <- c(-0.5, 5.5, 10.5, 13.5, 16.5, 19.5, 25.5, 32.5) # seven bins like RR98
bins <- c(-0.5, 10.5, 13.5, 16.5, 19.5, 32.5)
rr98$strength_bin <- cut(rr98$strength, breaks = bins, include.lowest = TRUE)
levels(rr98$strength_bin) <- as.character(1:7)

# aggregate data for response probability plot:
agg_rr98_bin <- rr98 %>% group_by(id, instruction, strength_bin) %>%
  summarise(n = n(), 
            dark = sum(response == "dark"),
            light = sum(response == "light")) %>%
  ungroup() %>%
   mutate(prop = binom.confint(dark, n, methods = "agresti-coull")[,"mean"],
     lower = binom.confint(dark, n, methods = "agresti-coull")$lower,
     upper = binom.confint(dark, n, methods = "agresti-coull")$upper)
  

knitr::kable(
  rr98 %>% group_by(id, instruction, strength_bin, response) %>%
    summarise(n = n()) %>%
    spread(strength_bin, n)
)


## ---- fig.height=4, fig.width=7------------------------------------------
xyplot(prop ~ strength_bin|id, agg_rr98_bin, group = instruction, type = "b", 
       auto.key = list(lines = TRUE), ylab = "Proportion of 'dark' responses")


## ---- fig.height=6, fig.width=7------------------------------------------

quantiles <- c(0.1, 0.3, 0.5, 0.7, 0.9)

## aggregate data for quantile plot
quantiles_rr98_bin <- rr98  %>% group_by(id, instruction, strength_bin) %>% 
  do(as.data.frame(t(quantile(.$rt, probs = quantiles)))) %>%
  ungroup() %>%
  gather(quantile, rt,  -id, -instruction, -strength_bin)
quantiles_rr98_bin$quantile <- factor(quantiles_rr98_bin$quantile, 
                                      levels = c("90%", "70%", "50%", "30%", "10%"))

xyplot(rt ~ strength_bin|id + instruction, quantiles_rr98_bin, group = quantile, type = "b", 
       auto.key = list(lines = TRUE), ylab = "RT (in seconds)", subset = instruction == "speed")

xyplot(rt ~ strength_bin|id + instruction, quantiles_rr98_bin, group = quantile, type = "b", 
       auto.key = FALSE, ylab = "RT (in seconds)", subset = instruction == "accuracy")




## ---- fig.height=6, fig.width=7------------------------------------------

agg2_rr98_response <- rr98  %>% group_by(id, instruction, strength_bin, response) %>% 
 do(as.data.frame(t(quantile(.$rt, probs = c(0.1, 0.3, 0.5, 0.7, 0.9))))) %>%
  ungroup() %>%
  gather(quantile, rt,  -id, - instruction, - strength_bin, -response)
agg2_rr98_response$quantile <- factor(agg2_rr98_response$quantile, 
                                      levels = c("90%", "70%", "50%", "30%", "10%"))

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
# objective function for diffusion with 1 a. loops over drift to assign drift rates to strength
objective_diffusion_separate <- function(pars, rt, boundary, drift, ...) {
  non_v_pars <- grep("^v", names(pars), invert = TRUE, value = TRUE)
  base_par <- length(non_v_pars)  # number of non-drift parameters
  densities <- vector("numeric", length(rt))
  for (i in seq_along(levels(drift))) {
    densities[drift == levels(drift)[i]] <- tryCatch(
      ddiffusion(rt[drift == levels(drift)[i]], boundary=boundary[drift == levels(drift)[i]], 
                 a=pars["a"], t0=pars["t0"],  
                 z=if ("z" %in% non_v_pars) pars["z"] else 0.5,
                 v=pars[base_par+i]), 
      error = function(e) 0)  
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
    t0_1 = runif(1, 0, 0.5),
    t0_2 = runif(1, 0, 0.5),
    z = runif(1, 0.4, 0.6),
    z = runif(1, 0.4, 0.6),
    z = runif(1, 0.4, 0.6)
  )
  start2 <- sort(rnorm(n_drift), decreasing = FALSE)
  names(start2) <- paste0("v_", seq_len(n_drift))
  c(start1[base_par], start2)
}

# function that tries different random start values until it works:
ensure_fit <- function(data, start_function, objective_function, base_pars, n_drift = 5) {
  
  start_ll <- 1e+06
  #browser()
  while(start_ll == 1e+06) {
    start <- start_function(base_pars)
    start_ll <- objective_function(start, 
                                   rt = data$rt, boundary = data$response_num, 
                                   drift = factor(data$strength_bin, seq_len(n_drift)), 
                                   instruction = data$instruction)
  }
  cat("\nstart fitting.\n") # just for information to see if it is stuck
  
  fit <- nlminb(start, objective_function, 
                rt = data$rt, boundary = data$response_num, 
                drift = factor(data$strength_bin, seq_len(n_drift)), 
                instruction = data$instruction,
                lower = c(rep(0, length(base_pars)), -Inf,
                          rep(-Inf, length(start_function(base_pars))-length(base_pars))))
  
  fit
}


## ---- echo=FALSE---------------------------------------------------------
load("rr98_wiener_fits.rda")



## ---- eval = FALSE-------------------------------------------------------
#  
#  fits_separate <- rr98 %>%
#    group_by(id, instruction) %>% # we loop across both, id and instruction
#    do(diffusion = ensure_fit(data = ., start_function = get_start,
#                              objective_function = objective_diffusion_separate,
#                              base_pars = c("a", "t0", "z"))) %>% ungroup()

## ------------------------------------------------------------------------
pars_separate <- fits_separate %>% group_by(id, instruction) %>% 
  do(as.data.frame(t(.$diffusion[[1]][["par"]]))) %>% ungroup() %>%
  as.data.frame()
pars_separate$ll <- (fits_separate %>% group_by(id, instruction) %>% 
                       do(ll = .$diffusion[[1]][["objective"]]) %>%  
                       summarize(ll2 = mean(ll[[1]])) %>% as.data.frame())[[1]]
if (!("z" %in% colnames(pars_separate))) pars_separate$z <- 0.5
knitr::kable(pars_separate)


## ----obtain_fits_not_run, eval = FALSE, include = FALSE------------------
#  fits_separate <- rr98 %>%
#    group_by(id, instruction) %>% # we loop across both, id and instruction
#    do(diffusion = ensure_fit(data = ., start_function = get_start,
#                              objective_function = objective_diffusion_separate,
#                              base_pars = c("a", "t0", "z"))) %>% ungroup()
#  
#  fits_separate_b <- rr98 %>%
#    group_by(id, instruction) %>% # we loop across both, id and instruction
#    do(diffusion = ensure_fit(data = ., start_function = get_start,
#                              objective_function = objective_diffusion_separate,
#                              base_pars = c("a", "t0", "z"))) %>% ungroup()
#  
#  
#  pars_separate_b <- fits_separate_b %>% group_by(id, instruction) %>%
#    do(as.data.frame(t(.$diffusion[[1]][["par"]]))) %>% ungroup() %>%
#    as.data.frame()
#  pars_separate_b$ll <- (fits_separate_b %>% group_by(id, instruction) %>%
#                         do(ll = .$diffusion[[1]][["objective"]]) %>%
#                         summarize(ll2 = mean(ll[[1]])) %>% as.data.frame())[[1]]
#  
#  
#  all.equal(pars_separate, pars_separate_b, tolerance = 0.001)
#  
#  fits_joint <- rr98 %>%
#    group_by(id) %>% # we loop only across id
#    do(diffusion = ensure_fit(data = ., start_function = get_start,
#                              objective_function = objective_diffusion_joint,
#                              base_pars = c("a_1", "a_2", "t0_1", "t0_2", "z"))) %>% ungroup()
#  
#  fits_joint_b <- rr98 %>%
#    group_by(id) %>% # we loop only across id
#    do(diffusion = ensure_fit(data = ., start_function = get_start,
#                              objective_function = objective_diffusion_joint,
#                              base_pars = c("a_1", "a_2", "t0_1", "t0_2", "z"))) %>% ungroup()
#  
#  pars_joint_b <- fits_joint_b %>% group_by(id) %>%
#    do(as.data.frame(t(.$diffusion[[1]][["par"]]))) %>%
#    ungroup() %>% as.data.frame()
#  
#  pars_joint_b$ll <- (fits_joint_b %>% group_by(id) %>%
#                      do(ll = .$diffusion[[1]][["objective"]]) %>%
#                      summarize(ll2 = mean(ll[[1]])) %>% as.data.frame())[[1]]
#  
#  all.equal(pars_joint, pars_joint_b, tolerance = 0.0001)
#  
#  save(fits_separate, fits_separate_b, fits_joint, fits_joint_b, file = "rr98_wiener_fits.rda")
#  

## ---- fig.height=5, fig.width=7, message=FALSE---------------------------


# get predicted response proportions
pars_separate_l <- pars_separate %>% gather("strength_bin", "v", starts_with("v"))
pars_separate_l$strength_bin <- factor(substr(pars_separate_l$strength_bin, 3,3), 
                                       levels = as.character(seq_len(length(bins)-1)))
#pars_separate_l <- inner_join(pars_separate_l, agg_rr98_bin)
pars_separate_l <- pars_separate_l  %>% group_by(id, instruction, strength_bin) %>%
  mutate(resp_prop = pdiffusion(rt=20, boundary="lower", a=a, v=v, t0=t0, z=z)) 

p1 <- xyplot(prop ~ strength_bin|id + instruction, agg_rr98_bin, type = "b", auto.key = 
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
  qdiffusion(quantiles*.$resp_prop, boundary="lower", 
             a=.$a, v=.$v, t0=.$t0, z=.$z)))) %>% 
  ungroup() %>% gather("quantiles", "dark", V1:V5)
separate_pred_light <- pars_separate_l %>% do(as.data.frame(t(
  qdiffusion(quantiles*(1-.$resp_prop), boundary="upper", 
             a=.$a, v=.$v, t0=.$t0, z=.$z)))) %>% 
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
objective_diffusion_joint <- function(pars, rt, boundary, drift, instruction) {
  non_v_pars <- grep("^v", names(pars), invert = TRUE, value = TRUE)
  base_par <- length(non_v_pars)  # number of non-drift parameters
  densities <- vector("numeric", length(rt))
  as <- c(pars["a_1"], pars["a_2"])
  ts <- c(pars["t0_1"], pars["t0_2"])
  for (j in seq_along(levels(instruction))) {
    for (i in seq_along(levels(drift))) {
      densities[drift == levels(drift)[i] & instruction == levels(instruction)[j]] <- tryCatch(
        ddiffusion(rt[drift == levels(drift)[i] & instruction == levels(instruction)[j]], 
                   boundary=boundary[drift==levels(drift)[i]&instruction==levels(instruction)[j]],
                   a=as[j], t0=ts[j], z=pars["z"], 
                   v=pars[base_par+i]), 
        error = function(e) 0)  
    }  
  }
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}



## ---- include=FALSE, eval=FALSE------------------------------------------
#  
#  fits_joint <- rr98 %>%
#    group_by(id) %>% # we loop only across id
#    do(diffusion = ensure_fit(data = ., start_function = get_start,
#                              objective_function = objective_diffusion_joint,
#                              base_pars = c("a_1", "a_2", "t0_1", "t0_2", "z"))) %>% ungroup()
#  

## ------------------------------------------------------------------------
pars_joint <- fits_joint %>% group_by(id) %>% 
  do(as.data.frame(t(.$diffusion[[1]][["par"]]))) %>% 
  ungroup() %>% as.data.frame()

pars_joint$ll <- (fits_joint %>% group_by(id) %>% 
                    do(ll = .$diffusion[[1]][["objective"]]) %>%  
                    summarize(ll2 = mean(ll[[1]])) %>% as.data.frame())[[1]]

knitr::kable(pars_joint)


## ---- message=FALSE------------------------------------------------------

ll_tab <- left_join(pars_joint[,c("id", "ll")], 
                    pars_separate %>% group_by(id) %>% summarise(ll_sep = sum(ll))) %>% 
  mutate(ll_diff_2 = 2*(ll-ll_sep), 
         p = round(pchisq(ll_diff_2, df = 5, lower.tail = FALSE), 4))
#rr98 %>% group_by(id) %>% summarise(n())
knitr::kable(ll_tab)


## ---- fig.height=5, fig.width=7, message=FALSE---------------------------


# get predicted response proportions
pars_joint_l <- pars_joint %>% gather("strength_bin", "v", starts_with("v")) %>% 
  gather("instruction_tmp", "pars", a_1, a_2, t0_1, t0_2) %>%
  separate(instruction_tmp, c("para", "instruction")) %>%
  spread(para, pars)
pars_joint_l$instruction <- factor(pars_joint_l$instruction, 
                                   levels = c("1", "2"), labels = c("speed", "accuracy"))
pars_joint_l$strength_bin <- factor(substr(pars_joint_l$strength_bin, 3,3), 
                                       levels = as.character(seq_len(length(bins)-1)))
#pars_separate_l <- inner_join(pars_separate_l, agg_rr98_bin)
pars_joint_l <- pars_joint_l  %>% group_by(id, instruction, strength_bin) %>%
  mutate(resp_prop = pdiffusion(rt=20, boundary="lower", a=a, v=v, t0=t0, z=z)) 

p1 <- xyplot(prop ~ strength_bin|id + instruction, agg_rr98_bin, type = "b", auto.key = 
               list(lines = TRUE), ylab = "Proportion of 'dark' responses", col = "grey")
p1e <- segplot(strength_bin ~ upper+lower|id + instruction, agg_rr98_bin, 
              auto.key = list(lines = TRUE), ylab = "Proportion of 'dark' responses", 
              col = "grey", horizontal = FALSE, segments.fun = panel.arrows,  
              draw.bands = FALSE, angle = 90, length = 0.05, ends = "both")
p2 <- xyplot(resp_prop ~ strength_bin|id + instruction, pars_separate_l, type = "b", 
             auto.key = list(lines = TRUE), ylab = "Proportion of 'dark' responses", 
             col = "darkgrey", lty = 3, pch = 0)
p3 <- xyplot(resp_prop ~ strength_bin|id + instruction, pars_joint_l, type = "b", 
             auto.key = list(lines = TRUE), ylab = "Proportion of 'dark' responses", 
             col = "black")
p2 + as.layer(p3) + as.layer(p1) + as.layer(p1e)


## ---- fig.height=6, fig.width=7, message=FALSE---------------------------

# get predicted quantiles (uses predicted response proportions)
joint_pred_dark <- pars_joint_l %>% do(as.data.frame(t(
  qdiffusion(quantiles*.$resp_prop, boundary="lower", 
             a=.$a, v=.$v, t0=.$t0, z=.$z)))) %>% 
  ungroup() %>% gather("quantiles", "dark", V1:V5)
joint_pred_light <- pars_joint_l %>% do(as.data.frame(t(
  qdiffusion(quantiles*(1-.$resp_prop), boundary="upper", 
             a=.$a, v=.$v, t0=.$t0, z=.$z)))) %>% 
  ungroup() %>% gather("quantiles", "light", V1:V5)

#joint_pred_light %>% filter(is.na(light))
joint_pred <- inner_join(joint_pred_dark, joint_pred_light)
joint_pred$quantiles <- factor(joint_pred$quantiles, 
                                  levels = c("V5", "V4", "V3", "V2", "V1"), 
                                  labels = c("90%", "70%", "50%", "30%", "10%"))
joint_pred <- joint_pred %>% gather("response", "rt", dark, light)

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
             scales = list(y = list(limits = c(0.25, 0.5))),
             col = "darkgrey", lty = 3, pch = 0)
p3 <- xyplot(rt ~ strength_bin|id + response, joint_pred, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "speed" & quantiles == "50%", 
             scales = list(y = list(limits = c(0.25, 0.5))))
p2 + as.layer(p3) + as.layer(p1) + as.layer(p1e)


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
             scales = list(y = list(limits = c(0.2, 1.5))),
             col = "darkgrey", lty = 3, pch = 0)
p3 <- xyplot(rt ~ strength_bin|id + response, joint_pred, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "accuracy" & quantiles == "50%", 
             scales = list(y = list(limits = c(0.25, 0.5))))
p2 + as.layer(p3) + as.layer(p1) + as.layer(p1e)


## ---- fig.height=7, fig.width=7------------------------------------------

p1 <- xyplot(rt ~ strength_bin|id+response, agg2_rr98_response, group = quantile, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "speed", layout = c(3,2), col = "grey")
p1e <- segplot(strength_bin ~ upper+lower|id+response, agg2_rr98_response, 
               auto.key = list(lines = TRUE), ylab = "Proportion of 'dark' responses", 
               col = "grey", horizontal = FALSE, segments.fun = panel.arrows,  
               draw.bands = FALSE, angle = 90, length = 0.05, ends = "both", 
               subset = instruction == "speed")
p2 <- xyplot(rt ~ strength_bin|id + response, joint_pred, group = quantiles, type = "b", 
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
p2 <- xyplot(rt ~ strength_bin|id + response, joint_pred, group = quantiles, type = "b", 
             auto.key = list(lines = TRUE), ylab = "RT (in seconds)", 
             subset = instruction == "accuracy", scales = list(y = list(limits = c(0.1, 3.0))))
p2 + as.layer(p1) + as.layer(p1e)


