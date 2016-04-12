
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

