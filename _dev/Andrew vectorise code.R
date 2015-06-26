system.time ({
  #### Stepwise parameter compare
  # 100000 rows, ~90%+ Unique: 214.25
  # 100000 rows, ~25%+ Unique:  64.36
  # 100000 rows, ~<1%  Unique:   1.23
key <- apply(params,1,paste,collapse=".")
ukey <- unique(key)
#out <- numeric(dim(params)[1])
for (i in ukey) {
  is.in <- key==i
  p <- params[is.in,,drop=FALSE][1,]

  cat ("i = ", i, "; is.in = ", is.in, "; p = ", p, "\n")
  
  #  print(p)
#  out[is.in] <- rrd()
  
}

})