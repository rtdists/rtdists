
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




rr98$v <- paste0("v_", rr98$strength)  # create variable for different vs.
rr98$a <- "a"  # one a variable

data <- rr98[rr98$id == "jf" & rr98$instruction == "speed",]

get_start <- function(base_par) {
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

# base_pars <- c("a", "t0", "sv")
# start <- get_start(base_pars)

start <- structure(c(0.604158850735985, 0.0458056061761454, 0.00925190048292279, 
-1.79185778153834, -1.06804427986316, -0.959099580750479, -0.752494006684124, 
-0.522666513065963, -0.442051777948383, -0.428041055056563, -0.392072403455471, 
-0.346123007678541, -0.259611586137565, -0.256023174564367, -0.225951166978148, 
-0.211362873684466, -0.205349428025187, -0.0880896039957638, 
-0.0433688042588261, 0.0405097283346233, 0.0447979702838185, 
0.0510138357198057, 0.0711055601913217, 0.145967577771945, 0.14680359507143, 
0.260736666421871, 0.286148885051787, 0.361230531385396, 0.421465015940815, 
0.707683960295186, 0.805551977975223, 0.813160237905473, 0.89062047456893, 
0.985388061208959, 2.09987520415512, 2.56456995332312), .Names = c("a", 
"t0", "sv", "v_0", "v_1", "v_2", "v_3", "v_4", "v_5", "v_6", 
"v_7", "v_8", "v_9", "v_10", "v_11", "v_12", "v_13", "v_14", 
"v_15", "v_16", "v_17", "v_18", "v_19", "v_20", "v_21", "v_22", 
"v_23", "v_24", "v_25", "v_26", "v_27", "v_28", "v_29", "v_30", 
"v_31", "v_32"))

## vectorized version:

# objective function for diffusion. names of a and v need to be passed as vectors of length(rt)
objective_diffusion <- function(pars, rt, boundary, v, a) {
  densities <- tryCatch(
    ddiffusion(rt, boundary=boundary, 
               a=pars[a], t0=pars["t0"], z=0.5, 
               sz=0.1, sv=pars["sv"], #st0=pars["st0"], 
               v=pars[v]), 
    error = function(e) 0)  
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}

# base_pars <- c("a", "t0", "sv")
# start <- get_start(base_pars)
# 
# objective <- 1e6
# while (objective == 1e6) {
#   base_pars <- c("a", "t0", "sv")
#   start <- get_start(base_pars)
#   objective <- objective_diffusion(start, rt = data$rt, boundary = data$response_num, v=data$v, a=data$a)
# }


objective_diffusion_2 <- function(pars, rt, boundary, v, a) {
  densities <- tryCatch(
    ddiffusion(rt, boundary=boundary, 
               a=pars["a"], t0=pars["t0"], z=0.5, 
               sz=0.1, sv=pars["sv"], #st0=pars["st0"], 
               v=pars[v]), 
    error = function(e) 0)  
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}


objective_diffusion_3 <- function(pars, rt, boundary, drift) 
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

# function that creates random start values
system.time(
  out_vector <- nlminb(start, objective_diffusion, 
                        rt = data$rt, boundary = data$response_num, v=data$v, a=data$a, 
                        lower = c(rep(0, length(base_pars)), rep(-Inf, length(start)-length(base_pars))))
)
#       User      System verstrichen 
#    2792.09        1.32     2806.85 

system.time(
  out_vector_2 <- nlminb(start, objective_diffusion_2, 
                        rt = data$rt, boundary = data$response_num, v=data$v, a=data$a, 
                        lower = c(rep(0, length(base_pars)), rep(-Inf, length(start)-length(base_pars))))
)
#       User      System verstrichen 
#    2568.96        0.32     2576.87 

system.time(
  out_loop_1 <- nlminb(start, objective_diffusion_3, 
                        rt = data$rt, boundary = data$response_num, drift = factor(data$strength), 
                        lower = c(rep(0, length(base_pars)), rep(-Inf, length(start)-length(base_pars))))
)
#       User      System verstrichen 
#     411.64        0.03      412.45 

all.equal(out_vector, out_vector_2) #TRUE
all.equal(out_loop_1, out_vector_2) #TRUE


