## ---- fig.height=4, fig.width=7------------------------------------------
require(rtdists)
require(dplyr)   # for data manipulations and looping
require(tidyr)   # for data manipulations
require(lattice) # for plotting and corresponding themes
require(latticeExtra)
lattice.options(default.theme = standard.theme(color = FALSE))
lattice.options(default.args = list(as.table = TRUE))
options(digits = 3) # only three decimal digits


data(rr98)
rr98 <- rr98[!rr98$outlier,]  #remove outliers

# aggregate data for first plot:
agg_rr98 <- rr98  %>% group_by(id, instruction, strength) %>% 
  summarise(prop = mean(response == "dark"), mean_rt = mean(rt), median_rt = mean(rt)) %>% 
  ungroup()

xyplot(prop ~ strength|id, prop_rr98, group = instruction, type = "b", auto.key = list(lines = TRUE), ylab = "Proportion of 'dark' responses")



## ---- fig.height=6, fig.width=7------------------------------------------

## aggregate data for quantile plot
quantiles_rr98 <- rr98  %>% group_by(id, instruction, strength) %>% 
  do(as.data.frame(t(quantile(.$rt, probs = c(0.1, 0.3, 0.5, 0.7, 0.9))))) %>%
  ungroup() %>%
  gather(quantile, rt,  -id, - instruction, - strength)
quantiles_rr98$quantile <- factor(quantiles_rr98$quantile, levels = c("90%", "70%", "50%", "30%", "10%"))

xyplot(rt ~ strength|id + instruction, quantiles_rr98, group = quantile, type = "b", auto.key = list(lines = TRUE), ylab = "RT (in seconds)", subset = instruction == "speed")

xyplot(rt ~ strength|id + instruction, quantiles_rr98, group = quantile, type = "b", auto.key = FALSE, ylab = "RT (in seconds)", subset = instruction == "accuracy")



## ---- fig.height=4, fig.width=7------------------------------------------

bins <- c(-0.5, 5.5, 10.5, 13.5, 16.5, 19.5, 25.5, 32.5)
rr98$strength_bin <- cut(rr98$strength, breaks = bins, include.lowest = TRUE)
levels(rr98$strength_bin) <- as.character(1:7)

# aggregate data for first plot:
agg_rr98_bin <- rr98  %>% group_by(id, instruction, strength_bin) %>% 
  summarise(prop = mean(response == "dark"), mean_rt = mean(rt), median_rt = mean(rt)) %>% 
  ungroup()

xyplot(prop ~ strength_bin|id, agg_rr98_bin, group = instruction, type = "b", auto.key = list(lines = TRUE), ylab = "Proportion of 'dark' responses")



## ---- fig.height=6, fig.width=7------------------------------------------

## aggregate data for quantile plot
quantiles_rr98_bin <- rr98  %>% group_by(id, instruction, strength_bin) %>% 
  do(as.data.frame(t(quantile(.$rt, probs = c(0.1, 0.3, 0.5, 0.7, 0.9))))) %>%
  ungroup() %>%
  gather(quantile, rt,  -id, - instruction, - strength_bin)
quantiles_rr98_bin$quantile <- factor(quantiles_rr98_bin$quantile, levels = c("90%", "70%", "50%", "30%", "10%"))

xyplot(rt ~ strength_bin|id + instruction, quantiles_rr98_bin, group = quantile, type = "b", auto.key = list(lines = TRUE), ylab = "RT (in seconds)", subset = instruction == "speed")

xyplot(rt ~ strength_bin|id + instruction, quantiles_rr98_bin, group = quantile, type = "b", auto.key = FALSE, ylab = "RT (in seconds)", subset = instruction == "accuracy")



## ---- fig.height=6, fig.width=7------------------------------------------

agg2_rr98_response <- rr98  %>% group_by(id, instruction, strength_bin, response) %>% 
 do(as.data.frame(t(quantile(.$rt, probs = c(0.1, 0.3, 0.5, 0.7, 0.9))))) %>%
  ungroup() %>%
  gather(quantile, rt,  -id, - instruction, - strength_bin, -response)
agg2_rr98_response$quantile <- factor(agg2_rr98_response$quantile, levels = c("90%", "70%", "50%", "30%", "10%"))

p1 <- xyplot(rt ~ strength_bin|id, agg2_rr98_response, group = quantile, type = "b", auto.key = list(lines = TRUE), ylab = "RT (in seconds)", subset = instruction == "speed" & response == "dark", layout = c(3,1))
p2 <- xyplot(rt ~ strength_bin|id, agg2_rr98_response, group = quantile, type = "b", auto.key = list(lines = TRUE), ylab = "RT (in seconds)", subset = instruction == "speed" & response == "light", col = "grey")
p1 + as.layer(p2)


p1 <- xyplot(rt ~ strength_bin|id, agg2_rr98_response, group = quantile, type = "b", auto.key = list(lines = TRUE), ylab = "RT (in seconds)", subset = instruction == "accuracy" & response == "dark", layout = c(3,1))
p2 <- xyplot(rt ~ strength_bin|id, agg2_rr98_response, group = quantile, type = "b", auto.key = list(lines = TRUE), ylab = "RT (in seconds)", subset = instruction == "accuracy" & response == "light", col = "grey")
p1 + as.layer(p2)



## ------------------------------------------------------------------------
# objective function for diffusion with one a. loops over drift to assign drift rates to strength
objective_diffusion_separate <- function(pars, rt, boundary, drift, ...) {
  base_par <- 3  # number of non-drift parameters
  densities <- vector("numeric", length(rt))
  for (i in seq_along(levels(drift))) {
    densities[drift == levels(drift)[i]] <- tryCatch(
      ddiffusion(rt[drift == levels(drift)[i]], boundary=boundary[drift == levels(drift)[i]], 
                 a=pars["a"], t0=pars["t0"], z=0.5, 
                 sz=0.1, sv=pars["sv"], #st0=pars["st0"], 
                 v=pars[base_par+i]), 
      error = function(e) 0)  
  }
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}


## ------------------------------------------------------------------------

# function that creates random start values
get_start <- function(base_par) {
  start1 <- c(
    a = runif(1, 0.5, 3),
    a_1 = runif(1, 0.5, 3), 
    a_2 = runif(1, 0.5, 3),
    t0 = runif(1, 0, 0.5), 
    z = runif(1, 0.4, 0.6), 
    sz = runif(1, 0, 0.5),
    sv = runif(1, 0, 0.5)
  )
  start2 <- sort(rnorm(7), decreasing = FALSE)
  names(start2) <- paste0("v_", 1:7)
  c(start1[base_par], start2)
}

# function that tries different random start values until it works:
ensure_fit <- function(data, start_function, objective_function, base_pars) {
  
  start_ll <- 1e+06
  while(start_ll == 1e+06) {
    start <- start_function(base_pars)
    start_ll <- objective_function(start, 
                                   rt = data$rt, boundary = data$response_num, 
                                   drift = factor(data$strength_bin, 1:7), instruction = data$instruction)
  }
  cat("\nstart fitting.\n") # just for information to see if it is stuck
  
  fit <- nlminb(start, objective_function, 
                rt = data$rt, boundary = data$response_num, 
                drift = factor(data$strength_bin, 1:7), instruction = data$instruction,
                lower = c(rep(0, length(base_pars)), 
                          rep(-Inf, length(start_function(base_pars))-length(base_pars))))
  
  fit
}


## ---- echo=FALSE---------------------------------------------------------
load("rr98_full-diffusion_fits.rda")


## ---- eval = FALSE-------------------------------------------------------
#  
#  fits_separate <- rr98 %>%
#    group_by(id, instruction) %>% # we loop across both, id and instruction
#    do(diffusion = ensure_fit(data = ., start_function = get_start,
#                              objective_function = objective_diffusion_separate,
#                              base_pars = c("a", "t0", "sv"))) %>% ungroup()
#  

## ------------------------------------------------------------------------
pars_separate <- as.data.frame(fits_separate %>% group_by(id, instruction) %>% do(as.data.frame(t(.$diffusion[[1]][["par"]]))) %>% ungroup())
pars_separate$ll <- (fits_separate %>% group_by(id, instruction) %>% do(ll = .$diffusion[[1]][["objective"]]) %>%  summarize(ll2 = mean(ll[[1]])) %>% as.data.frame())[[1]]

knitr::kable(pars_separate)


## ------------------------------------------------------------------------
objective_diffusion_joint <- function(pars, rt, boundary, drift, instruction) {
  base_par <- 4
  densities <- vector("numeric", length(rt))
  as <- c(pars["a_1"], pars["a_2"])
  for (j in seq_along(levels(instruction))) {
    for (i in seq_along(levels(drift))) {
      densities[drift == levels(drift)[i] & instruction == levels(instruction)[j]] <- tryCatch(
        ddiffusion(rt[drift == levels(drift)[i] & instruction == levels(instruction)[j]], 
                   boundary=boundary[drift==levels(drift)[i]&instruction==levels(instruction)[j]],
                   a=as[j], t0=pars["t0"], z=0.5, 
                   sz=0.1, sv=pars["sv"], #st0=pars["st0"], 
                   v=pars[base_par+i]), 
        error = function(e) 0)  
    }  
  }
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}



## ---- include=FALSE, eval=FALSE------------------------------------------
#  
#  fits_separate <- rr98 %>%
#    group_by(id, instruction) %>% # we loop across both, id and instruction
#    do(diffusion = ensure_fit(data = ., start_function = get_start,
#                              objective_function = objective_diffusion_separate,
#                              base_pars = c("a", "t0", "sv"))) %>% ungroup()
#  

## ------------------------------------------------------------------------
pars_joint <- as.data.frame(fits_joint %>% group_by(id) %>% do(as.data.frame(t(.$diffusion[[1]][["par"]]))) %>% ungroup())

pars_joint$ll <- (fits_joint %>% group_by(id) %>% do(ll = .$diffusion[[1]][["objective"]]) %>%  summarize(ll2 = mean(ll[[1]])) %>% as.data.frame())[[1]]

knitr::kable(pars_joint)


## ----obtain_fits_not_run, eval = FALSE, include = FALSE------------------
#  
#  fits_separate <- rr98 %>%
#    group_by(id, instruction) %>% # we loop across both, id and instruction
#    do(diffusion = ensure_fit(data = ., start_function = get_start,
#                              objective_function = objective_diffusion_separate,
#                              base_pars = c("a", "t0", "sv"))) %>% ungroup()
#  
#  
#  fits_joint <- rr98 %>%
#    group_by(id) %>% # we loop only across instruction
#    do(diffusion = ensure_fit(data = ., start_function = get_start,
#                              objective_function = objective_diffusion_joint,
#                              base_pars = c("a_1", "a_2", "t0", "sv"))) %>% ungroup()
#  
#  
#  save(fits_separate, fits_joint, file = "rr98_full-diffusion_fits.rda")
#  

## ---- message=FALSE------------------------------------------------------

ll_tab <- left_join(pars_joint[,c("id", "ll")], pars_separate %>% 
                      group_by(id) %>% summarise(ll_sep = sum(ll))) %>% 
  mutate(ll_diff_2 = 2*(ll-ll_sep), p = afex::round_ps(pchisq(ll_diff_2, df = 9, lower.tail = FALSE)))
#rr98 %>% group_by(id) %>% summarise(n())
knitr::kable(ll_tab)


