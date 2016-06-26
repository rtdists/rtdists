
## generate random LBA data:
rt1 <- rLBA(500, A=0.5, b=1, t0 = 0.5, mean_v=c(2.4, 1.6), sd_v=c(1,1.2))
head(rt1)
prop.table(table(rt1$response))

# original parameters have 'high' log-likelihood:
sum(log(dLBA(rt1$rt, rt1$response, A=0.5, b=1, t0 = 0.5, mean_v=c(2.4, 1.6), sd_v=c(1,1.2))))

# data can also be passed as data.frame (same is true for pLBA):
sum(log(dLBA(rt1, A=0.5, b=1, t0 = 0.5, mean_v=c(2.4, 1.6), sd_v=c(1,1.2))))

objective_fun <- function(par, rt, response, distribution = "norm") {
  # simple parameters
  spar <- par[!grepl("[12]$", names(par))]  
  
  # distribution parameters:
  dist_par_names <- unique(sub("[12]$", "", grep("[12]$" ,names(par), value = TRUE)))
  dist_par <- vector("list", length = length(dist_par_names))
  names(dist_par) <- dist_par_names
  for (i in dist_par_names) dist_par[[i]] <- as.list(unname(par[grep(i, names(par))]))
  dist_par$sd_v <- c(1, dist_par$sd_v) # fix first sd to 1

  # get summed log-likelihood:
  d <- do.call(dLBA, args = c(rt=list(rt), response=list(response), spar, dist_par, 
                               distribution=distribution, silent=TRUE))
  if (any(d < 0e-10)) return(1e6) 
  else return(-sum(log(d)))
}

# gives same value as manual calculation above:
objective_fun(c(A=0.5, b=1, t0=0.5, mean_v1=2.4, mean_v2=1.6, sd_v2=1.2), 
              rt=rt1$rt, response=rt1$response)

\dontrun{
# can we recover the parameters?
init_par <- runif(6)
names(init_par) <- c("A", "b", "t0", "mean_v1", "mean_v2",  "sd_v2")
nlminb(objective_fun, start = init_par, rt=rt1$rt, response=rt1$response, lower = 0)
}


# plot cdf (2 accumulators):
curve(pLBA(x, response = 1, A=0.5, b=1, t0 = 0.5, mean_v=c(2.4, 1.6), sd_v=c(1,1.2)), 
     xlim = c(0, 2), ylim = c(0,1), 
     ylab = "cumulative probability", xlab = "response time",
     main = "Defective CDFs of LBA")
curve(pLBA(x, response = 2, A=0.5, b=1, t0 = 0.5, mean_v=c(2.4, 1.6), sd_v=c(1,1.2)), 
     add=TRUE, lty = 2)
legend("topleft", legend=c("1", "2"), title="Response", lty=1:2)


# plot cdf (3 accumulators):
curve(pLBA(x, response = 1, A=0.5, b=1, t0 = 0.5, mean_v=c(2.4, 1.6, 1.0), sd_v=c(1,1.2, 2.0)), 
     xlim = c(0, 2), ylim = c(0,1), 
     ylab = "cumulative probability", xlab = "response time",
     main = "Defective CDFs of LBA")
curve(pLBA(x, response = 2, A=0.5, b=1, t0 = 0.5,  mean_v=c(2.4, 1.6, 1.0), sd_v=c(1,1.2, 2.0)), 
     add=TRUE, lty = 2)
curve(pLBA(x, response = 3, A=0.5, b=1, t0 = 0.5,  mean_v=c(2.4, 1.6, 1.0), sd_v=c(1,1.2, 2.0)), 
     add=TRUE, lty = 3)
legend("topleft", legend=c("1", "2", "3"), title="Response", lty=1:2)


## qLBA can only return values up to maximal predicted probability:
pLBA(Inf, response = 1, A=0.5, b=1, t0 = 0.5, mean_v=c(2.4, 1.6), sd_v=c(1,1.2))
# [1] 0.6604696

qLBA(0.66, response = 1, A=0.5, b=1, t0 = 0.5, mean_v=c(2.4, 1.6), sd_v=c(1,1.2))
# 2.559532

qLBA(0.67, response = 1, A=0.5, b=1, t0 = 0.5, mean_v=c(2.4, 1.6), sd_v=c(1,1.2))
# NA

# qLBA also accepts a data.frame as first argument:
t <- data.frame(p = rep(c(0.05, 0.1, 0.66), 2), response = rep(1:2, each = 3))
#      p response
# 1 0.05        1
# 2 0.10        1
# 3 0.87        1
# 4 0.05        2
# 5 0.10        2
# 6 0.87        2
qLBA(t, A=0.5, b=1, t0 = 0.5, mean_v=c(2.4, 1.6), sd_v=c(1,1.2))
