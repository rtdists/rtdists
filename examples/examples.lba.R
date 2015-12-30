
## generate random LBA data:
rt1 <- riLBA(500, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3))
head(rt1)
prop.table(table(rt1$response))

# original parameters have 'high' log-likelihood:
sum(log(diLBA(rt1$rt, rt1$response, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3))))
sum(log(diLBA(rt1$rt, rt1$response, A=0.5, b=1, t0 = 0.5, mean_v=c(1.5, 0.6), sd_v=c(0.2,0.3))))

\dontrun{
# can we recover the parameters?
objective_fun <- function(par, rt, response, distribution = "norm") {
  # simple parameters
  spar <- par[!grepl("[12]$", names(par))]  
  
  # distribution parameters:
  dist_par_names <- unique(sub("[12]$", "", grep("[12]$" ,names(par), value = TRUE)))
  dist_par <- vector("list", length = length(dist_par_names))
  names(dist_par) <- dist_par_names
  for (i in dist_par_names) dist_par[[i]] <- unname(par[grep(i, names(par))])
  
  # get summed log-likelihood:
  d <- do.call(diLBA, args = c(rt=list(rt), response=list(response), spar, dist_par, 
                               distribution=distribution))
  if (any(d == 0)) return(1e6)
  else return(-sum(log(d)))
}

objective_fun(c(A=0.5, b=1, t0=0.5, mean_v1=1.2, mean_v2=1.0, sd_v1=0.2, sd_v2=0.3), 
              rt=rt1$rt, response=rt1$response)

init_par <- runif(7)
names(init_par) <- c("A", "b", "t0", "mean_v1", "mean_v2", "sd_v1", "sd_v2")
nlminb(objective_fun, start = init_par, rt=rt1$rt, response=rt1$response, lower = 0)

}
