
data(rr98)
rr98 <- rr98[!rr98$outlier,]  #remove outliers
head(rr98)
#   id session block trial instruction source strength response response_num correct    rt outlier
# 1 jf       2     1    21    accuracy   dark        8     dark            1    TRUE 0.801   FALSE
# 2 jf       2     1    22    accuracy   dark        7     dark            1    TRUE 0.680   FALSE
# 3 jf       2     1    23    accuracy  light       19    light            2    TRUE 0.694   FALSE
# 4 jf       2     1    24    accuracy   dark       21    light            2   FALSE 0.582   FALSE
# 5 jf       2     1    25    accuracy  light       19     dark            1   FALSE 0.925   FALSE
# 6 jf       2     1    26    accuracy   dark       10     dark            1    TRUE 0.605   FALSE


\dontrun{}

## fit original model: separately for each condition
require(dplyr) # for convenient looping and data.frame manipulation

rr98$v <- paste0("v_", rr98$strength)  # create variable for different vs.

# objective function for diffusion with one a. loops over drift to assign drift rates to strength
objective_diffusion_separate <- function(pars, rt, boundary, drift, ...) {
  base_par <- 3
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
  start2 <- sort(rnorm(33), decreasing = FALSE)
  names(start2) <- paste0("v_", 0:32)
  c(start1[base_par], start2)
}

# function that tries different random start values until it works:
ensure_fit <- function(data, start_function, objective_function, base_pars) {
  fit <- list(objective = 1e+06)
  try <- 1
  while(fit$objective == 1e+06) {
    cat("try: ", try, "\n") # just for information to see if it is stuck
    try <- try+1
    fit <- nlminb(start_function(base_pars), objective_function, 
                  rt = data$rt, boundary = data$response_num, 
                  drift = factor(data$strength, 0:32), instruction = data$instruction,
                  lower = c(rep(0, length(base_pars)), 
                            rep(-Inf, length(start_function(base_pars))-length(base_pars))))
  }
  fit
}

fits_separate <- rr98 %>% 
  group_by(id, instruction) %>% # we loop across both, id and instruction
  do(diffusion = ensure_fit(data = ., start_function = get_start, 
                            objective_function = objective_diffusion_separate, 
                            base_pars = c("a", "t0", "sv"))) %>% ungroup()


fits_joint <- rr98 %>% 
  group_by(id) %>% # we loop only across instruction
  do(diffusion = ensure_fit(data = ., start_function = get_start, 
                            objective_function = objective_diffusion_joint, 
                            base_pars = c("a_1", "a_2", "t0", "sv"))) %>% ungroup()

pars_separate <- as.data.frame(fits_separate %>% group_by(id, instruction) %>% do(as.data.frame(t(.$diffusion[[1]][["par"]]))) %>% ungroup())
ll_separate <- fits_separate %>% group_by(id, instruction) %>% do(as.data.frame(.$diffusion[[1]][["objective"]])) %>% ungroup() %>% as.data.frame() 
colnames(ll_separate)[3] <- "ll"
ll_separate <- ll_separate %>% group_by(id) %>% summarise(lls = sum(ll))


ll_joint <- as.data.frame(fits_joint %>% group_by(id) %>% do(as.data.frame(.$diffusion[[1]][["objective"]])) %>% ungroup())
colnames(ll_joint)[2] <- "lls"

(ll_diff <- ll_joint$lls - ll_separate$lls)

rr98 %>% group_by(id) %>% summarise(n())

pars_joint <- as.data.frame(fits_joint %>% group_by(id) %>% do(as.data.frame(t(.$diffusion[[1]][["par"]]))) %>% ungroup())

## check if fit can be recovered:

# fits_separate_b <- rr98 %>% 
#   group_by(id, instruction) %>% # we loop across both, id and instruction
#   do(diffusion = ensure_fit(data = ., start_function = get_start, 
#                             objective_function = objective_diffusion_separate, 
#                             base_pars = c("a", "t0", "sv"))) %>% ungroup()
# 
# pars_separate_a <- as.data.frame(fits_separate %>% group_by(id, instruction) %>% do(as.data.frame(t(.$diffusion[[1]][["par"]]))) %>% ungroup())
# pars_separate_b <- as.data.frame(fits_separate_b %>% group_by(id, instruction) %>% do(as.data.frame(t(.$diffusion[[1]][["par"]]))) %>% ungroup())
# 
# ll_separate_a <- as.data.frame(fits_separate %>% group_by(id, instruction) %>% do(ll = .$diffusion[[1]][["objective"]]) %>% ungroup())
# ll_separate_b <- as.data.frame(fits_separate_b %>% group_by(id, instruction) %>% do(ll = .$diffusion[[1]][["objective"]]) %>% ungroup())
# 
# all.equal(pars_separate_a, pars_separate_b, tolerance = 0.0001)  # TRUE
# 
# all.equal(ll_separate_a, ll_separate_b) # TRUE


# fits_joint_b <- rr98 %>% 
#   group_by(id) %>% # we loop only across instruction
#   do(diffusion = ensure_fit(data = ., start_function = get_start, 
#                             objective_function = objective_diffusion_joint, 
#                             base_pars = c("a_1", "a_2", "t0", "sv"))) %>% ungroup()
# 
# pars_joint_a <- as.data.frame(fits_joint %>% group_by(id) %>% do(as.data.frame(t(.$diffusion[[1]][["par"]]))) %>% ungroup())
# pars_joint_b <- as.data.frame(fits_joint_b %>% group_by(id) %>% do(as.data.frame(t(.$diffusion[[1]][["par"]]))) %>% ungroup())
# 
# ll_joint_a <- as.data.frame(fits_joint %>% group_by(id) %>% do(ll = .$diffusion[[1]][["objective"]]) %>% ungroup())
# ll_joint_b <- as.data.frame(fits_joint_b %>% group_by(id) %>% do(ll = .$diffusion[[1]][["objective"]]) %>% ungroup())
# 
# all.equal(pars_joint_a, pars_joint_b, tolerance = 0.0001)  # [1] "Component “v_32”: Mean relative difference: 0.2833519"
# all.equal(ll_joint_a, ll_joint_b) # TRUE


objective_diffusion_orig <- function(pars, rt, boundary, v, ...) {
  densities <- tryCatch(
    ddiffusion(rt, boundary=boundary, 
               a=pars["a_1"], t0=pars["t0"], z=0.5, 
               sz=0.1, sv=pars["sv"], #st0=pars["st0"], 
               v=pars[v]), 
    error = function(e) 0)  
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}


with(rr98, table(correct, strength, instruction, id))

sort(unique(rr98$strength))
round(seq(0, 32, length.out = 8))
cut

if (requireNamespace("dplyr")) {
  rr98 %>% group_by(id, instruction, correct) %>% do(as.data.frame(t(as.matrix(quantile(.$rt, c(0.1, 0.3, 0.5, 0.7, 0.9))))))
  
  rr98 %>% group_by(id, instruction, correct) %>% do(as.data.frame(t(quantile(.$rt, c(0.1, 0.3, 0.5, 0.7, 0.9)))))
  
  require(tidyr)
  rr98 %>% group_by(id, instruction, correct, strength) %>% summarise(counts = n()) %>% spread(strength, counts)
}

objective_diffusion <- function(pars, rt, boundary, drift) 
{
  base_par <- 5
  densities <- vector("numeric", length(rt))
  for (i in seq_along(levels(drift))) {
    densities[drift == levels(drift)[i]] <- tryCatch(
      ddiffusion(rt[drift == levels(drift)[i]], boundary=boundary[drift == levels(drift)[i]], 
                 a=pars["a"], t0=pars["t0"], z=pars["z"], 
                 sz=pars["sz"], sv=pars["sv"], #st0=pars["st0"], 
                 v=pars[base_par+i]), 
      error = function(e) 0)  
  }
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}

#require(optimx)

get_start_1 <- function(base_par) {
  start1 <- c(a = runif(1, 0.5, 3), 
              t0 = runif(1, 0, 0.5), 
              z = runif(1, 0.4, 0.6), 
              sz = runif(1, 0, 0.5),
              sv = runif(1, 0, 0.5)
  )
  start2 <- sort(rnorm(33), decreasing = FALSE)
  names(start2) <- paste0("v_", 0:32)
  c(start1[base_par], start2)
}





ensure_fit <- function(data, start_function, objective_function, base_pars) {
  fit <- list(objective = 1e+06)
  try <- 1
  while(fit$objective == 1e+06) {
    print(try)
    try <- try+1
    fit <- nlminb(start_function(base_pars), objective_function, 
                        rt = data$rt, boundary = data$boundary, drift = factor(data$strength, 0:32),
                        lower = c(rep(0, length(base_pars)), rep(-Inf, length(start-length(base_pars)))))
  }
  fit
}


objective_diffusion_orig <- function(pars, rt, boundary, drift) 
{
  base_par <- 3
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

fits_1_a <- rr98 %>% 
  group_by(id, instruction) %>% 
  do(diffusion = ensure_fit(data = ., start_function = get_start_1, objective_function = objective_diffusion, 
                            base_pars = c("a", "t0", "z", "sz", "sv"))) %>% ungroup()

fits_1_b <- rr98 %>% 
  group_by(id, instruction) %>% 
  do(diffusion = ensure_fit(data = ., start_function = get_start_1, objective_function = objective_diffusion, 
                            base_pars = c("a", "t0", "z", "sz", "sv"))) %>% ungroup()

pars_1_a <- as.data.frame(fits_1_a %>% group_by(id, instruction) %>% do(as.data.frame(t(.$diffusion[[1]][["par"]]))) %>% ungroup())
pars_1_b <- as.data.frame(fits_1_b %>% group_by(id, instruction) %>% do(as.data.frame(t(.$diffusion[[1]][["par"]]))) %>% ungroup())

ll_1_a <- as.data.frame(fits_1_a %>% group_by(id, instruction) %>% do(ll = .$diffusion[[1]][["objective"]]) %>% ungroup())
ll_1_b <- as.data.frame(fits_1_b %>% group_by(id, instruction) %>% do(ll = .$diffusion[[1]][["objective"]]) %>% ungroup())

all.equal(pars_1_a, pars_1_b, tolerance = 0.01)

round(abs(pars_1_a[,-(1:2)] - pars_1_b[,-(1:2)]), 3)

fits_orig_a <- rr98 %>% 
  group_by(id, instruction) %>% 
  do(diffusion = ensure_fit(data = ., start_function = get_start_1, objective_function = objective_diffusion_orig, 
                            base_pars = c("a", "t0", "sv"))) %>% ungroup()

fits_orig_b <- rr98 %>% 
  group_by(id, instruction) %>% 
  do(diffusion = ensure_fit(data = ., start_function = get_start_1, objective_function = objective_diffusion_orig, 
                            base_pars = c("a", "t0", "sv"))) %>% ungroup()

pars_orig_a <- as.data.frame(fits_orig_a %>% group_by(id, instruction) %>% do(as.data.frame(t(.$diffusion[[1]][["par"]]))) %>% ungroup())
pars_orig_b <- as.data.frame(fits_orig_b %>% group_by(id, instruction) %>% do(as.data.frame(t(.$diffusion[[1]][["par"]]))) %>% ungroup())
all.equal(pars_orig_a, pars_orig_b, tolerance = 0.0001)


ll_orig_a <- as.data.frame(fits_orig_a %>% group_by(id, instruction) %>% do(ll = .$diffusion[[1]][["objective"]]) %>% ungroup())
ll_orig_b <- as.data.frame(fits_orig_b %>% group_by(id, instruction) %>% do(ll = .$diffusion[[1]][["objective"]]) %>% ungroup())

save(fits_orig_a, fits_orig_a, fits_1_a, fits_1_b, file = "rr_fits.rda")

# prep combined fit



### old

ensure_fit <- function(data, start_function, objective_function) {
  fit <- list(objective = 1e+06)
  try <- 1
  while(fit$objective == 1e+06) {
    print(try)
    try <- try+1
    fit <- nlminb(start_function(), objective_diffusion, 
                        rt = data$rt, boundary = data$boundary, drift = factor(data$strength, 0:32),
                        lower = c(rep(0, 5), rep(-Inf, length(start-5))))
  }
  fit
}

fits_1 <- rr98 %>% 
  group_by(id, instruction) %>% 
  do(diffusion = ensure_fit(data = ., start_function = get_start_1, objective_function = objective_diffusion))

pars <- fits_1 %>% group_by(id, instruction) %>% do(as.data.frame(t(.$diffusion[[1]][["par"]]))) %>% ungroup()


plot(1, 1, xlim = c(0, 32), ylim = c(-6, 6), type = "n")
ltys <- rep(c(1,2,3), each = 2)
pchs <- rep(c(0,8,18), each = 2)
for (i in seq_len(nrow(pars))) {
  lines(0:32, as.data.frame(pars %>% select(v_0:v_32))[i,], type = "o", lty = ltys[i], pch = pchs[i])
}

