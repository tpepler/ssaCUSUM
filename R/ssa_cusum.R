#********************************************************************************
# Function to perform a SSA-CUSUM

ssa_cusum <- function(x, dates = NULL, L = length(x) %/% 4, q = 1, verbose = TRUE, diagplots = FALSE,
                    MRL0 = 200, h.reps = 10000,
                    #sim = FALSE,
                    main = NULL, eps.out = NULL){

  # Check to see if dates are NULL -- if so, simply number the observations
  if(is.null(dates)){
    dates <- c(1:length(x))
  }

  # Perform SSA on the time series
  ssa.obj <- ssa_calc(x = x, dates = dates, L = L, q = q, diagplots = diagplots)

  # Calculate alert threshold, h
  #h<-h.calc(resids=(ssa.obj$obs-ssa.obj$exp)[1:L],MRL0=MRL0,reps=h.reps)
  h <- h_calc(resids = (ssa.obj$obs - ssa.obj$exp)[1:(L*2)], MRL0 = MRL0, reps = h.reps)
  #h<-h.calc(resids=((ssa.obj$obs-ssa.obj$exp)/sqrt(ssa.obj$exp))[1:L],MRL0=MRL0,reps=h.reps)

  # Calculate the CUSUM
  #cusum.obj<-cusum.calc(obs=ssa.obj$obs,exp=ssa.obj$exp,dates=ssa.obj$dates,q=ssa.obj$q,L=ssa.obj$L,h=h)
  cusum.obj <- cusum_calc(ssa.obj = ssa.obj, h = h)

  if(verbose){
    # Print CUSUM alarm threshold value
    cat(paste('CUSUM alarm threshold value: ', round(h, 2), '\n', sep = ''))

    # Plot the CUSUM
    #plot_cusum(cusum.obj, main = main, eps.out = eps.out)
    plot_cusum(cusum.obj, main = main)

    # Alarm dates
    cat('Alarm dates:\n')
    print(dates[cusum.obj$alarms])
  }

  # Output for simulation experiments
  #if(sim){return(list(obs = cusum.obj$obs, exp = cusum.obj$exp, alarms = cusum.obj$alarms))}
  #if(!verbose){return(cusum.obj)}
  return(invisible(cusum.obj))
}

#********************************************************************************
