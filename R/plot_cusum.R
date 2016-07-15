#********************************************************************************
# Function to plot a CUSUM chart

plot_cusum <- function(cusum.obj, main = NULL, eps.out = NULL){
  # cusum.obj: CUSUM object, as output from cusum.calc

  # Calculate outbreak start and end dates
  outbreak.start <- NULL
  outbreak.end <- NULL
  for(i in 2:length(cusum.obj$obs)){
    if(cusum.obj$outbreak[i - 1] == FALSE & cusum.obj$outbreak[i] == TRUE){
      outbreak.start <- append(outbreak.start, cusum.obj$dates[i])
    }
    else {
      if(cusum.obj$outbreak[i - 1] == TRUE & cusum.obj$outbreak[i] == FALSE){
        outbreak.end <- append(outbreak.end, cusum.obj$dates[i - 1])
      }
    }
  }
  # In case of outbreak spanning the full time series
  if(is.null(outbreak.start) & cusum.obj$outbreak[1] == TRUE){
    outbreak.start <- cusum.obj$dates[1]
  }
  # For outbreaks continuing to the end of the time series
  if(!is.null(outbreak.start)){
    if(length(outbreak.end) < length(outbreak.start)){
      outbreak.end <- append(outbreak.end, cusum.obj$dates[length(cusum.obj$obs)])
    }
  }

  # Plot the CUSUM chart
  oldpar <- par()
  #layout(matrix(c(1, 2), ncol = 1))
  par(#mfrow = c(2, 1),
    mar = c(3, 4.5, 4, 1), las = 1,
    bty = 'o')
  if(!is.null(eps.out)){
    setEPS()
    postscript(eps.out[1], width = 8, height = 5)
  }
  if(!is.null(main)){cusum.title <- paste(main, ': Control chart', sep = '')}
  else{cusum.title <- 'Control chart'}
  if(sum(cusum.obj$outbreak) > 0){
    plot(x = c(cusum.obj$dates), y = c(cusum.obj$cusumvals), type = 'n',
         xlab = '', ylab = 'CUSUM', main = cusum.title,
         ylim = c(0, max(cusum.obj$h, cusum.obj$cusumvals) * 1.1),
         panel.first = rect(xleft = outbreak.start, ybottom = -1e6, xright = outbreak.end, ytop = 1e6, col = 'gray80', border = NA))
  }
  else{
    plot(x = c(cusum.obj$dates), y = c(cusum.obj$cusumvals), type = 'n',
         xlab = '', ylab = 'CUSUM', main = cusum.title,
         ylim = c(0, max(cusum.obj$h, cusum.obj$cusumvals) * 1.1))
  }
  lines(x = cusum.obj$dates, y = cusum.obj$cusumvals)
  points(x = cusum.obj$dates, y = cusum.obj$cusumvals, pch = 20)
  abline(h = 0, lty = 2)
  abline(h = cusum.obj$h, lty = 2, col = 'red')
  #points(x=cusum.obj$dates[cusum.obj$outbreak],y=cusum.obj$cusumvals[cusum.obj$outbreak],pch=20,col='darkorange')
  points(x = cusum.obj$dates[cusum.obj$alarms], y = cusum.obj$cusumvals[cusum.obj$alarms], pch = 4, col = 'red', lwd = 3)
  abline(v = cusum.obj$dates[cusum.obj$alarms], lty = 3, col = 'red')
  if(!is.null(eps.out)){dev.off()}

  # Plot the time series data
  par(mar = c(5, 4.5, 4, 1))
  if(!is.null(eps.out)){
    setEPS()
    postscript(eps.out[2], width = 8, height = 5)
  }
  if(!is.null(main)){cusum.title <- paste(main, ': Time series (q = ', cusum.obj$q, ')', sep = '')}
  else{cusum.title <- paste('Time series (q = ', cusum.obj$q, ')', sep = '')}
  if(sum(cusum.obj$outbreak) > 0){
    plot(cusum.obj$dates,cusum.obj$obs,
         type = 'S', lty = 1, col = 'black',
         main = cusum.title,
         xlab = '',
         ylab = 'Number of events',
         panel.first = rect(xleft = outbreak.start, ybottom = -1e6, xright = outbreak.end, ytop = 1e6, col = 'gray80', border = NA))
  }
  else{
    plot(cusum.obj$dates, cusum.obj$obs,
         type = 'S', lty = 1, col = 'black',
         main = cusum.title,
         xlab = '',
         ylab = 'Number of events')
  }
  lines(cusum.obj$dates[1:(cusum.obj$L * 2)], cusum.obj$exp[1:(cusum.obj$L * 2)], type = 'l', col = 'blue', lwd = 2)
  lines(cusum.obj$dates[(cusum.obj$L * 2):length(cusum.obj$obs)], cusum.obj$exp[(cusum.obj$L * 2):length(cusum.obj$obs)],
        type = 'l', col = 'blue', lwd = 2, lty = 3)
  if(sum(cusum.obj$outbreak) > 0){
    abline(v = cusum.obj$dates[cusum.obj$alarms], col = 'red', lwd = 2)
  }
  #temp<-locator(1)
  #print(temp)
  legend(x = cusum.obj$dates[1],
         y = min(cusum.obj$obs) - (max(cusum.obj$obs) - min(cusum.obj$obs)) / 5,
         #y = temp$y,
         legend = c('Observed counts', 'Estimated signal', 'Predicted signal', 'Alarm', 'Out-of-control'),
         lty = c(1, 1, 3, 1, 1), lwd = c(1, 2, 2, 2, 12), col = c('black', 'blue', 'blue', 'red', 'gray80'),
         cex = 0.7,
         #fill=c(NA,NA,NA,NA,'gray80'),
         border = NA,
         bty = 'n',
         xpd = TRUE,
         horiz = TRUE)

  par(#mfrow = oldpar$mfrow,
      mar = oldpar$mar, las = oldpar$las,
      bty = oldpar$bty)
  if(!is.null(eps.out)){dev.off()}
}

#********************************************************************************
