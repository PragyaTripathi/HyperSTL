#' @title vis.fitted
#' @name vis.fitted
#' @description Creates a graph of fitted values
#' @export
vis.fitted <- function(y, y.hat) {
  plot(y, type='l', col='blue', 
        ylim=c(min(y[!is.na(y)], y.hat[!is.na(y.hat)]), 
               max(y[!is.na(y)], y.hat[!is.na(y.hat)])))
  lines(y.hat, col='red')
  legend("topright", 
         col=c("red", "blue"), 
         legend=c("Fitted", "Observed"), 
         lty=c(1,1))
}

#' @title vis.rainfall
#' @description Plots accumulative rainfall and given depths
#' @export
vis.rainfall <- function(data, depths) {
  mmin <- Inf
  mmax <- -Inf
  for(d in depths) {
    mmin <- min(mmin, data[[d]])
    mmax <- max(mmax, data[[d]])
  }
  
  par(mar=c(5,4,4,5))
  
  for(i in 1:length(depths)) {
    if(i == 1) {
      plot(data$date, data[[depths[[i]]]], 
           col=1, 
           ylim=c(mmin, mmax),
           ylab='Soil moisture',
           xlab="date",
           type='l')  
    }
    
    lines(data$date, data[[depths[[i]]]], col=i)
  }
  par(new=TRUE)
  
  plot(data$date, data$rainfall, type='l', 
       col=length(depths) + 1,
       ylab="",
       xlab="",
       axes=FALSE)
  
  axis(4)
  axis(side=4, at = pretty(range(data$rainfall)))
  mtext("rainfall",side=4,line=3)
  legend("topright",col=1:length(depths),
         lty=1, legend=depths)
}

#' @title vis.soilmoisture
#' @description Visualizes soil moisture time series
#' @export
vis.soilmoisture <- function(data, depths) {
  mmin <- Inf
  mmax <- -Inf
  for(d in depths) {
    mmin <- min(mmin, data[[d]])
    mmax <- max(mmax, data[[d]])
  }
  
  for(i in 1:length(depths)) {
    if(i == 1) {
      plot(data$date, data[[depths[[i]]]], 
           col=1, 
           ylim=c(mmin, mmax),
           ylab='Soil moisture',
           xlab="date",
           type='l')  
    }
    
    lines(data$date, data[[depths[[i]]]], col=i)
  }
  legend("topright",col=1:length(depths),
         lty=1, legend=depths)
}
