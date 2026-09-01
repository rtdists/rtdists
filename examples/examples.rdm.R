
## generate random RDM data:
rt1 <- rRDM(500, A=0.5, b=1, t0 = 0.5, v=c(2.4, 1.6))
head(rt1)
prop.table(table(rt1$response))

# original parameters have 'high' log-likelihood:
sum(log(dRDM(rt1$rt, rt1$response, A=0.5, b=1, t0 = 0.5, v=c(2.4, 1.6))))

# data can also be passed as data.frame (same is true for pRDM):
sum(log(dRDM(rt1, A=0.5, b=1, t0 = 0.5, v=c(2.4, 1.6))))

objective_fun <- function(par, rt, response) {
  d <- dRDM(rt, response, A=par["A"], b=par["b"], t0=par["t0"],
            v=c(par["v1"], par["v2"]), silent=TRUE)
  if (any(d < 0e-10)) return(1e6)
  else return(-sum(log(d)))
}

# gives same value as manual calculation above:
objective_fun(c(A=0.5, b=1, t0=0.5, v1=2.4, v2=1.6),
              rt=rt1$rt, response=rt1$response)

\dontrun{
# can we recover the parameters?
# should be run several times with different random values of init_par
init_par <- runif(5)
init_par[2] <- sum(init_par[1:2]) # ensures b is larger than A
init_par[3] <- runif(1, 0, min(rt1$rt)) # ensures t0 is not too large
names(init_par) <- c("A", "b", "t0", "v1", "v2")
nlminb(objective_fun, start = init_par, rt=rt1$rt, response=rt1$response, lower = 0)
}

# plot defective cdfs (2 accumulators):
curve(pRDM(x, response = 1, A=0.5, b=1, t0 = 0.5, v=c(2.4, 1.6), silent=TRUE),
      xlim = c(0, 2), ylim = c(0, 1),
      ylab = "cumulative probability", xlab = "response time",
      main = "Defective CDFs of the RDM")
curve(pRDM(x, response = 2, A=0.5, b=1, t0 = 0.5, v=c(2.4, 1.6), silent=TRUE),
      add=TRUE, lty = 2)
legend("topleft", legend=c("1", "2"), title="Response", lty=1:2)

# quantiles of the first response:
qRDM(c(0.1, 0.5, 0.9), response = 1, A=0.5, b=1, t0 = 0.5, v=c(2.4, 1.6),
     scale_p = TRUE)
