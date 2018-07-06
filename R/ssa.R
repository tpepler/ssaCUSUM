ssa <- function(tsdata, # univariate time series data
                window, # window length for SSA,
                train.periods, # number of time periods (at start of time series) used for training the model
                q = NULL, # number of singular spectra to retain, if known
                prop.var.exp = NULL, # minimum proportion of variance accounted for by SSA model, if required
                verbose = FALSE, # logical: print outputs of analysis?
                plot = FALSE # logical: produce plots for SSA?
  ){
  # Check that window parameter value is appropriate
  if(window < 2){stop('Error: window parameter must be at least 2.')}
  if(window > floor(length(tsdata)/2)){stop('Error: window parameter value must not exceed half the length of the time series.')}

  # Create empty lists for the SSA data objects
  ts.obj <- vector('list', 3)
  names(ts.obj) <- c('traj', 'meanvec', 'svd')

  # Embed the time series
  ts.obj$traj <- ssaCUSUM::ts_stack(x = as.numeric(tsdata[1:train.periods]), window = window) # training data
  traj.full <- ssaCUSUM::ts_stack(x = as.numeric(tsdata), window = window) # full time series data set

  # Centre the rows of the trajectory matrices
  ts.obj$meanvec <- apply(ts.obj$traj, 1, mean)
  ts.obj$traj <- ts.obj$traj - ts.obj$meanvec
  traj.full <- traj.full - ts.obj$meanvec

  # SVD of the trajectory matrices multiplied by their transpose matrices
  ts.obj$svd <- svd(ts.obj$traj %*% t(ts.obj$traj))

  # Proportion of variance accounted for by each of the SSA components
  if(verbose){
    print(round(ts.obj$svd$d / sum(ts.obj$svd$d), 3))
  }

  # Cumulative proportion of variance accounted for by each of the SSA components
  if(verbose){
    print(round(cumsum(ts.obj$svd$d / sum(ts.obj$svd$d)), 3))
  }

  if(plot){
    # Proportion variance explained by each of the eigenvectors
    plot(x = c(1:window),
         y = c(ts.obj$svd$d / sum(ts.obj$svd$d)),
         type = 'l',
         ylab = 'Proportion variance explained',
         xlab = 'Eigenvalue',
         main = 'SSA eigenvalues')
    points(x = c(1:window),
           y = c(ts.obj$svd$d / sum(ts.obj$svd$d)),
           pch = 16)
  }

  # Set value for q if prop.var.exp is specified
  if(!is.null(prop.var.exp)){
    q <- which(cumsum(ts.obj$svd$d / sum(ts.obj$svd$d)) >= prop.var.exp)[1]
  }

  # If neither q nor prop.var.exp are specified, prompt user for input
  if(is.null(q)){
    q <- as.numeric(readline(prompt = paste('Specify number of components to retain for SSA model (should be between 1 and ', window, '): ', sep = '')))
  }

  # Reconstructing the time series using SSA
  #q <- 4 # number of SSA components used for reconstruction (selected from inspection of a scree plot of the squared singular values)
  B.red <- ts.obj$svd$v[, 1:q]
  P.red <- B.red %*% t(B.red)
  res.ssa <- P.red %*% traj.full + outer(ts.obj$meanvec, rep(1, times = ncol(traj.full)))
  res.ssa <- ssaCUSUM::hankelize(res.ssa)
  pred.ssa <- ssaCUSUM::ts_unstack(res.ssa)

  # Calculate the residuals
  resid.ssa <- tsdata - pred.ssa

  # Correlogram of the residuals
  #acf(x = resid.ssa, lag.max = 36, main = 'Correlogram: SSA residuals')

  if(plot){
    # Plot of the original and reconstructed time series
    plot(x = c(1:length(tsdata)),
         y = tsdata,
         xlab = 'Time',
         ylab = 'Time series values',
         main = 'Original and reconstructed time series',
         type = 'l')
    lines(x = c(1:length(tsdata)),
          y = pred.ssa,
          lwd = 2,
          col = 'red')
    legend('topleft',
           legend = c('original', 'reconstructed'),
           lwd = c(1, 2),
           col = c('black', 'red'))

    # Plot of the residuals vs. time
    plot(x = c(1:length(tsdata)),
         y = as.numeric(resid.ssa),
         xlab = 'Time',
         ylab = 'SSA residual',
         main = 'SSA residuals vs. time',
         type = 'p',
         pch = 20)
    abline(h = 0, lty = 2, col = 'darkgray')
    abline(h = 2*sd(resid.ssa), lty = 3, col = 'red')
    abline(h = -2*sd(resid.ssa), lty = 3, col = 'red')
  }

  # Create output object to be returned to user
  output.obj <- vector('list', 7)
  names(output.obj) <- c('tsdata', 'fitted', 'values', 'vectors', 'q', 'projmat', 'rowmeans')
  output.obj$tsdata <- tsdata # original time series
  output.obj$fitted <- pred.ssa # fitted/predicted values
  output.obj$values <- ts.obj$svd$d # singular values
  output.obj$vectors <- ts.obj$svd$v # singular vectors
  output.obj$q <- q # number of components retained in SSA model
  output.obj$projmat <- P.red # projection matrix for the SSA model
  output.obj$rowmeans <- ts.obj$meanvec # row means of the training data trajectory matrix

  return(output.obj)
}
