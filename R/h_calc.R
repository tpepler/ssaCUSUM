#********************************************************************************
# Function to calculate the alert threshold for SSA CUSUM charts

h_calc <- function(resids, MRL0 = 200, reps = 10000){
  # reps: number of replications used to estimate the value of h
  # MRL0: average in-control run length

  #print(mean(resids))
  #resids<-resids-mean(resids)

  hvals <- rep(NA, times = reps)
  for(r in 1:reps){
    data0 <- resids[sample(1:length(resids), size = MRL0, replace = TRUE)] * sample(c(-1, 1), size = MRL0, replace = TRUE)
    cusumvals <- rep(NA, times = MRL0)
    cusumvals[1] <- max(0, data0[1])
    for(i in 2:MRL0){
      cusumvals[i] <- max(0, cusumvals[i - 1] + data0[i])
    }
    hvals[r] <- max(cusumvals)
  }
  return(median(hvals))
  #return(mean(hvals))
}

#********************************************************************************
