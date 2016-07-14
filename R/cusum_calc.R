#********************************************************************************
# Function to calculate the CUSUM values

cusum_calc <- function(ssa.obj, h = h_calc(resids = ssa.obj$obs - ssa.obj$exp)){
  n <- length(ssa.obj$obs)
  cusumvals <- rep(NA, times = n)
  alarms <- rep(FALSE, times = n)
  outbreak <- rep(FALSE, times = n)
  outbreak.current <- FALSE
  cusum.current <- 0
  for(i in 1:n){
    cusumvals[i] <- max(0, cusum.current + ssa.obj$obs[i] - ssa.obj$exp[i])
    #cusumvals[i]<-max(0,(cusum.current+obs[i]-exp[i])/sqrt(exp[i])) # deviation in standard deviations
    cusum.current <- cusumvals[i]

    # If threshold is breached
    if(cusumvals[i] >= h){
      cusum.current <- max(0, ssa.obj$obs[i] - ssa.obj$exp[i]) # reset the CUSUM after the breach
      #cusum.current<-max(0,(cusum.current+obs[i]-exp[i])/sqrt(exp[i])) # reset the CUSUM after the breach
      # Raise alarm if not currently in outbreak scenario
      if(outbreak.current == FALSE){
        alarms[i] <- TRUE
        outbreak.current <- TRUE
        #cusum.current<-max(0,obs[i]-exp[i]) # reset the CUSUM after alarm
        # Indicate outbreak from last time CUSUM was zero
        j <- i
        while(j > 0 & if(j > 0){cusumvals[j] > 0}else{FALSE}){
          outbreak[j] <- TRUE
          j <- j - 1
        }
      }
    }

    # If threshold was not breached
    else{
      # Outbreak indicator ceases if CUSUM reaches zero
      if(cusumvals[i] == 0 & outbreak.current == TRUE){
        outbreak.current <- FALSE
      }
    }
    outbreak[i] <- outbreak.current
  }
  return(list(dates = ssa.obj$dates, obs = ssa.obj$obs, exp = ssa.obj$exp, q = ssa.obj$q,
              L = ssa.obj$L, h = h, cusumvals = cusumvals, alarms = alarms, outbreak = outbreak))
}
