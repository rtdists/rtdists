


objective_fun <- function(par, rt, response, distribution = "norm") {
  # simple parameters
  spar <- par[!grepl("[12]$", names(par))]  
  
  # distribution parameters:
  dist_par_names <- unique(sub("[12]$", "", grep("[12]$" ,names(par), value = TRUE)))
  dist_par <- vector("list", length = length(dist_par_names))
  names(dist_par) <- dist_par_names
  for (i in dist_par_names) dist_par[[i]] <- as.list(unname(par[grep(i, names(par))]))
  dist_par$sd_v <- c(1, dist_par$sd_v) # fix sd_v[1] to 1 for identification

  # get summed log-likelihood:
  d <- do.call(diLBA, args = c(rt=list(rt), response=list(response), spar, dist_par, 
                               distribution=distribution, silent=TRUE))
  if (any(d < 0e-10)) return(1e6)
  else return(-sum(log(d)))
}


# recovery of simple model
rt1 <- riLBA(500, A=0.5, b=1, t0 = 0.5, mean_v=c(2.4, 1.6), sd_v=c(1,1.2))
head(rt1)
prop.table(table(rt1$response))
objective_fun(c(A=0.5, b=1, t0=0.5, mean_v1=2.4, mean_v2=1.6, sd_v1=1.2), 
              rt=rt1$rt, response=rt1$response)

init_par <- runif(6)
names(init_par) <- c("A", "b", "t0", "mean_v1", "mean_v2", "sd_v1")
nlminb(objective_fun, start = init_par, rt=rt1$rt, response=rt1$response, lower = 0)

nlminb(objective_fun, start = c(A=0.5, b=1, t0=0.5, mean_v1=2.4, mean_v2=1.6, sd_v1=1.2), rt=rt1$rt, response=rt1$response, lower = 0)


# test recovery of t0 per accunmulator
rt1 <- riLBA(5000, A=0.5, b=1, t0 = list(0.5, 0.4), mean_v=c(2.8, 1.6), sd_v=c(1.0,1.2), distribution = "norm")
rt1 <- rlba_norm(5000, A=0.5, b=1, t0 = 0.5, mean_v=c(2.8, 1.6), sd_v=c(1.0,1.2))
head(rt1)
prop.table(table(rt1$response))

objective_fun(c(A=0.5, b=1, t01=0.5, t02=0.4, mean_v1=2.8, mean_v2=1.6, sd_v2=1.2), 
              rt=rt1$rt, response=rt1$response)


init_par <- c(runif(2), runif(2, 0, 0.3), runif(3))
names(init_par) <- c("A",  "b", "t01", "t02", "mean_v1", "mean_v2", "sd_v2")
(out1 <- nlminb(objective_fun, start = init_par, rt=rt1$rt, response=rt1$response, lower = 0))

nlminb(objective_fun, start = c(A=0.5, b=1, t01=0.5, t02=0.4, mean_v1=2.8, mean_v2=1.6, sd_v2=1.2), rt=rt1$rt, response=rt1$response, lower = 0)


## all params per accumulator:
rt1 <- riLBA(5000, A=list(0.5, 0.6), b=list(1, 1.2), t0 = list(0.5, 0.4), mean_v=c(2.8, 1.4), sd_v=c(1.0,1.2), distribution = "norm")
head(rt1)
prop.table(table(rt1$response))

objective_fun(c(A1=0.5, A2=0.6, b1=1, b2=1.2, t01=0.5, t02=0.4, mean_v1=2.8, mean_v2=1.4, sd_v2=1.2), 
              rt=rt1$rt, response=rt1$response)

init_par <- c(runif(4), runif(2, 0, 0.3), runif(3))
names(init_par) <- c("A1", "A2", "b1", "b2", "t01", "t02", "mean_v1", "mean_v2", "sd_v2")
nlminb(objective_fun, start = init_par, rt=rt1$rt, response=rt1$response, lower = 0)

nlminb(objective_fun, start = c(A1=0.5, A2=0.6, b1=1, b2=1.2, t01=0.5, t02=0.4, mean_v1=2.8, mean_v2=1.4, sd_v2=1.2), rt=rt1$rt, response=rt1$response, lower = 0)



