
data(speed_acc)
str(speed_acc)

# remove excluded trials:
speed_acc <- droplevels(speed_acc[!speed_acc$censor,])

# new factors for obtaining values as in Table 1, Wagenmakers et al. (2008, p. 152)
speed_acc$freq <- with(speed_acc, 
                       factor(ifelse(stim_cat == "nonword", "nonword", 
                                     as.character(frequency)), 
                              levels = c("high", "low", "very_low", "nonword")))
# corr = correct (0 = correct, 1 = error)
speed_acc$corr <- with(speed_acc, 1-as.numeric(stim_cat == response))

str(speed_acc)

## aggregated RTs:
aggregate(rt ~ condition + freq + corr, speed_acc, mean)
## Error Rate:
aggregate(corr ~ condition + freq + corr, speed_acc, mean)

