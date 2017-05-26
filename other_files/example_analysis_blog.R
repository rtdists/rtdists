require(rtdists)

# Exp. 1; Wagenmakers, Ratcliff, Gomez, & McKoon (2008, JML)
data(speed_acc)   
# remove excluded trials:
speed_acc <- droplevels(speed_acc[!speed_acc$censor,]) 
# create numeric response variable with 1 = error and 2 = correct response: 
speed_acc$corr <- with(speed_acc, as.numeric(stim_cat == response))+1 
# select participant 11, accuracy condition, non-word trials only
p11 <- speed_acc[speed_acc$id == 11 & 
                   speed_acc$condition == "accuracy" & 
                   speed_acc$stim_cat == "nonword",] 
prop.table(table(p11$corr))
#          1          2 
# 0.04166667 0.95833333 

ll_lba <- function(pars, rt, response) {
  d <- dLBA(rt = rt, response = response, 
            A = pars["A"], 
            b = pars["A"]+pars["b"], 
            t0 = pars["t0"], 
            mean_v = pars[c("v1", "v2")], 
            sd_v = c(1, pars["sv"]), 
            silent=TRUE)
  if (any(d == 0)) return(1e6)
  else return(-sum(log(d)))
}

start <- c(runif(3, 0.5, 3), runif(2, 0, 0.2), runif(1))
names(start) <- c("A", "v1", "v2", "b", "t0", "sv")
p11_norm <- nlminb(start, ll_lba, lower = c(0, -Inf, 0, 0, 0, 0), 
                   rt=p11$rt, response=p11$corr)
p11_norm[1:3]
# $par
#          A         v1         v2          b         t0         sv 
#  0.1182940 -2.7409230  1.0449963  0.4513604  0.1243441  0.2609968 
# 
# $objective
# [1] -211.4202
# 
# $convergence
# [1] 0


ll_diffusion <- function(pars, rt, response) 
{
  densities <- ddiffusion(rt, response=response, 
                          a=pars["a"], 
                          v=pars["v"], 
                          t0=pars["t0"], 
                          sz=pars["sz"], 
                          st0=pars["st0"],
                          sv=pars["sv"])
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}

start <- c(runif(2, 0.5, 3), 0.1, runif(3, 0, 0.5))
names(start) <- c("a", "v", "t0", "sz", "st0", "sv")
p11_diff <- nlminb(start, ll_diffusion, lower = 0, 
                   rt=p11$rt, response=p11$corr)
p11_diff[1:3]
# $par
#         a         v        t0        sz       st0        sv 
# 1.3206011 3.2727202 0.3385602 0.4621645 0.2017950 1.0551706 
# 
# $objective
# [1] -207.5487
# 
# $convergence
# [1] 0


# quantiles:
q <- c(0.1, 0.3, 0.5, 0.7, 0.9)

## observed data:
(p11_q_c <- quantile(p11[p11$corr == 2, "rt"], probs = q))
#    10%    30%    50%    70%    90% 
# 0.4900 0.5557 0.6060 0.6773 0.8231 
(p11_q_e <- quantile(p11[p11$corr == 1, "rt"], probs = q))
#    10%    30%    50%    70%    90% 
# 0.4908 0.5391 0.5905 0.6413 1.0653 

### LBA:
# predicted error rate  
(pred_prop_correct_lba <- pLBA(Inf, 2, 
                               A = p11_norm$par["A"], 
                               b = p11_norm$par["A"]+p11_norm$par["b"], 
                               t0 = p11_norm$par["t0"], 
                               mean_v = c(p11_norm$par["v1"], p11_norm$par["v2"]), 
                               sd_v = c(1, p11_norm$par["sv"])))
# [1] 0.9581342

(pred_correct_lba <- qLBA(q*pred_prop_correct_lba, response = 2, 
                          A = p11_norm$par["A"], 
                          b = p11_norm$par["A"]+p11_norm$par["b"], 
                          t0 = p11_norm$par["t0"], 
                          mean_v = c(p11_norm$par["v1"], p11_norm$par["v2"]), 
                          sd_v = c(1, p11_norm$par["sv"])))
# [1] 0.4871710 0.5510265 0.6081855 0.6809796 0.8301286
(pred_error_lba <- qLBA(q*(1-pred_prop_correct_lba), response = 1, 
                        A = p11_norm$par["A"], 
                        b = p11_norm$par["A"]+p11_norm$par["b"], 
                        t0 = p11_norm$par["t0"], 
                        mean_v = c(p11_norm$par["v1"], p11_norm$par["v2"]), 
                        sd_v = c(1, p11_norm$par["sv"])))
# [1] 0.4684374 0.5529575 0.6273737 0.7233961 0.9314820


### diffusion:
# same result as when using Inf, but faster:
(pred_prop_correct_diffusion <- pdiffusion(rt = 20,  response = "upper",
                                      a=p11_diff$par["a"], 
                                      v=p11_diff$par["v"], 
                                      t0=p11_diff$par["t0"], 
                                      sz=p11_diff$par["sz"], 
                                      st0=p11_diff$par["st0"], 
                                      sv=p11_diff$par["sv"]))  
# [1] 0.964723

(pred_correct_diffusion <- qdiffusion(q*pred_prop_correct_diffusion, 
                                      response = "upper",
                                      a=p11_diff$par["a"], 
                                      v=p11_diff$par["v"], 
                                      t0=p11_diff$par["t0"], 
                                      sz=p11_diff$par["sz"], 
                                      st0=p11_diff$par["st0"], 
                                      sv=p11_diff$par["sv"]))
# [1] 0.4748271 0.5489903 0.6081182 0.6821927 0.8444566
(pred_error_diffusion <- qdiffusion(q*(1-pred_prop_correct_diffusion), 
                                    response = "lower",
                                    a=p11_diff$par["a"], 
                                    v=p11_diff$par["v"], 
                                    t0=p11_diff$par["t0"], 
                                    sz=p11_diff$par["sz"], 
                                    st0=p11_diff$par["st0"], 
                                    sv=p11_diff$par["sv"]))
# [1] 0.4776565 0.5598018 0.6305120 0.7336275 0.9770047


### plot predictions

par(mfrow=c(1,2), cex=1.2)
plot(p11_q_c, q*prop.table(table(p11$corr))[2], pch = 2, ylim=c(0, 1), xlim = c(0.4, 1.3), ylab = "Cumulative Probability", xlab = "Response Time (sec)", main = "LBA")
points(p11_q_e, q*prop.table(table(p11$corr))[1], pch = 2)
lines(pred_correct_lba, q*pred_prop_correct_lba, type = "b")
lines(pred_error_lba, q*(1-pred_prop_correct_lba), type = "b")
legend("right", legend = c("data", "predictions"), pch = c(2, 1), lty = c(0, 1))

plot(p11_q_c, q*prop.table(table(p11$corr))[2], pch = 2, ylim=c(0, 1), xlim = c(0.4, 1.3), ylab = "Cumulative Probability", xlab = "Response Time (sec)", main = "Diffusion")
points(p11_q_e, q*prop.table(table(p11$corr))[1], pch = 2)
lines(pred_correct_diffusion, q*pred_prop_correct_diffusion, type = "b")
lines(pred_error_diffusion, q*(1-pred_prop_correct_diffusion), type = "b")


png("p11_predictions_correct.png", width = 600, height = 400, units = "px")
par(mfrow=c(1,2), cex=1.2)
plot(p11_q_c, q*prop.table(table(p11$corr))[2], pch = 2, ylim=c(0, 1), xlim = c(0.4, 1.3), ylab = "Cumulative Probability", xlab = "Response Time (sec)", main = "LBA")
points(p11_q_e, q*prop.table(table(p11$corr))[1], pch = 2)
lines(pred_correct_lba, q*pred_prop_correct_lba, type = "b")
lines(pred_error_lba, q*(1-pred_prop_correct_lba), type = "b")
legend("right", legend = c("data", "predictions"), pch = c(2, 1), lty = c(0, 1))

plot(p11_q_c, q*prop.table(table(p11$corr))[2], pch = 2, ylim=c(0, 1), xlim = c(0.4, 1.3), ylab = "Cumulative Probability", xlab = "Response Time (sec)", main = "Diffusion")
points(p11_q_e, q*prop.table(table(p11$corr))[1], pch = 2)
lines(pred_correct_diffusion, q*pred_prop_correct_diffusion, type = "b")
lines(pred_error_diffusion, q*(1-pred_prop_correct_diffusion), type = "b")
dev.off()

png("p11_CDF_check.png", width = 600, height = 400, units = "px")
par(cex=1.5)
rand_rts <- rdiffusion(1e5, a=p11_diff$par["a"], 
                            v=p11_diff$par["v"], 
                            t0=p11_diff$par["t0"], 
                            sz=p11_diff$par["sz"], 
                            st0=p11_diff$par["st0"], 
                            sv=p11_diff$par["sv"])
plot(ecdf(rand_rts[rand_rts$response == "upper","rt"]), lwd = 2)

normalised_pdiffusion = function(rt,...) pdiffusion(rt,...)/pdiffusion(rt=Inf,...) 
curve(normalised_pdiffusion(x, response = "upper",
                            a=p11_diff$par["a"], 
                            v=p11_diff$par["v"], 
                            t0=p11_diff$par["t0"], 
                            sz=p11_diff$par["sz"], 
                            st0=p11_diff$par["st0"], 
                            sv=p11_diff$par["sv"]), 
      add=TRUE, col = "yellow", lty = 2, lwd = 2)
dev.off()

