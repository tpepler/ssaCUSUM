#********************************************************************************
# Function to perform SSA

ssa_calc <- function(x, dates, L = length(x) %/% 4, q = 1, diagplots = TRUE){

  n <- length(x)
  # Embed the time series (only use the first L*2 data points!)
  x.traj <- ts_stack(x[1:(L * 2)], window = L)

  # Centre the columns of the trajectory matrices
  x.meanvec <- apply(x.traj, 2, mean)
  x.traj <- scale(x.traj, center = T, scale = F)

  # SVD of the trajectory matrices multiplied by their transpose matrices
  x.svd2 <- svd(x.traj %*% t(x.traj))

  # Reconstructing the time series prospectively (from time point 2L onwards) using SSA
  B.red <- x.svd2$v[, 1:q]
  P.red <- B.red %*% t(B.red)

  x.2L.mean <- mean(x[1:(L * 2)])
  x.pred.ssa <- rep(NA, times = n)

  # Calculate predicted values for the first 2L time points (already observed)
  xi <- x[1:(2 * L)]
  ni <- length(xi)
  x.part.traj <- ts_stack(c(xi[L:2], xi, xi[(ni-1):(ni - L + 1)]), window = L) # add (L-1) neighbouring months to beginning and end of series, to improve predicted values?
  x.2L.meanvec <- rep(x.2L.mean, times = ncol(x.part.traj))
  # Centre the columns of the trajectory matrix
  x.part.traj <- scale(x.part.traj, center = T, scale = F)
  x.red.ssa <- outer(rep(1, times = nrow(x.part.traj)), x.2L.meanvec) + P.red %*% x.part.traj
  x.red.ssa <- hankelize(x.red.ssa)
  x.pred.ssa.temp <- ts_unstack(x.red.ssa)
  x.pred.ssa[1:(2 * L)] <- x.pred.ssa.temp[L:(3 * L - 1)]

  for(i in (2 * L + 1):n){
    # PROPOSAL 1: Use full data set
    # Embed the time series
    xi <- x[1:i]
    ni <- length(xi)
    x.part.traj <- ts_stack(c(xi[L:2], xi, xi[(ni - 1):(ni - L + 1)]), window = L) # add (L-1) neighbouring months to beginning and end of series, to improve predicted values?
    x.2L.meanvec <- rep(x.2L.mean, times = ncol(x.part.traj))
    # Centre the columns of the trajectory matrix
    x.part.traj <- scale(x.part.traj, center = T, scale = F)
    x.red.ssa <- outer(rep(1, times = nrow(x.part.traj)), x.2L.meanvec) + P.red %*% x.part.traj
    x.red.ssa <- hankelize(x.red.ssa)
    x.pred.ssa.temp <- ts_unstack(x.red.ssa)
    x.pred.ssa[i] <- x.pred.ssa.temp[L - 1 + i]
  }

  #x.pred.ssa[x.pred.ssa < 0] <- 0 # set negative predictive values to the minimum possible count of zero

  # Calculate the residuals
  x.resid.ssa <- x - x.pred.ssa

  if(diagplots){
    oldpar <- par()
    par(mfrow = c(2, 2), bty = 'n')
    # Scree plot of the eigenvalues
    plot(1:length(x.svd2$d), x.svd2$d, type = 'l',
         main = 'SSA: Eigenvalue scree plot',
         xlab = 'Eigenvalue number',
         ylab = 'Eigenvalue',
         col = 'black')
    points(1:length(x.svd2$d), x.svd2$d, pch = 20, cex = 0.8)

    # Histograms of the eigenvalues
    #hist(x.svd2$d, breaks = 40, col = 'darkgray', main = 'Eigenvalues distribution', xlab = 'Eigenvalue')

    # Proportion of variation accounted for by the first q eigenvectors
    cat(paste('Proportion of variation accounted for by the first q = ', q, ' eigenvector(s): ',
              round(sum(x.svd2$d[1:q]) / sum(x.svd2$d), 2), '\n', sep = ''))

    # Plot: residuals vs. time
    plot(dates, x.resid.ssa,
         type = 'p', lty = 3, col = 'black', pch = 20,
         main = paste('SSA residuals (q = ', q, ')', sep = ''),
         ylab = 'Residuals', xlab = 'Time')
    abline(h = 0, lty = 3)

    # Histograms of residuals
    hist(x.resid.ssa, col = 'darkgray', freq = FALSE,
         main = paste('SSA residuals (q = ', q, ')', sep = ''),
         xlab = 'Residuals')
    lines(density(x.resid.ssa), lwd = 2)
    abline(v = 0, lty = 3)

    # Q-Q plot of residuals
    #qqnorm(x.resid.ssa, pch = 20)
    #qqline(x.resid.ssa)

    # Correlogram of the SSA residuals
    acf(x = x.resid.ssa, lag.max = L, main = 'SSA residuals: Correlogram')

    par(mfrow = oldpar$mfrow, bty = oldpar$bty)
  }

  return(list(q = q, L = L, dates = dates, obs = x, exp = x.pred.ssa))
}

#********************************************************************************
