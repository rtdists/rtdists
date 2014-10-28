
data(speed_acc)
str(speed_acc)

# remove excluded trials:
speed_acc <- droplevels(speed_acc[!speed_acc$censor,])
str(speed_acc)
