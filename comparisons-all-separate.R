# library(AWI)
#library(TSA)
library(TTR)
library(Metrics)

generateComparisonsPlot <- function(data, depth, hyper_stl_output = NULL,
  image_scale = 2.5, 
  plot_resolution = 1, legend_x_coord = 150, 
  y_axis_label =  '') {
    plot_line_widths = 3
    legend_axes_labels_font = 1.3

  if (0) {
  setupEpsFileForPlot("sma-fit.eps", image_scale);
  plotObserved(data, depth)
  plotSMA(data, depth)
  dev.off()
 
  setupEpsFileForPlot("spline-fit.eps", image_scale);
  plotObserved(data, depth)
  plotSpline(data, depth)
  dev.off()

   setupEpsFileForPlot("wma-fit.eps", image_scale);
  plotObserved(data, depth)
  plotWMA(data, depth)
  dev.off()

  setupEpsFileForPlot("loess-fit.eps", image_scale);
  plotObserved(data, depth)
  plotLOESS(data, depth)
  dev.off()
}
if (1) {
  if (!is.null(d)) {
    setupEpsFileForPlot("hyper-stl-infra.eps", image_scale);
    plotObserved(data, depth, y_axis_label = y_axis_label)
    plotHyperStl(hyper_stl_output, data)
    dev.off()
  }
  

  setupEpsFileForPlot("spline-fit-infra.eps", image_scale);
  plotObserved(data, depth, y_axis_label = y_axis_label)
  plotSpline(data, depth)
  dev.off()


  setupEpsFileForPlot("peak-preserving-fit-infra.eps", image_scale);
  plotObserved(data, depth, y_axis_label = y_axis_label)
  plotPeakPreserving(data, depth)
  dev.off()
}

 }


plotHyperStl <- function (y, data, plot_resolution = 1, plot_line_width = 3, 
  legend_axes_labels_font = 1.3, line_color = 'blue') {
  x <- data$date
  date_stamps = x[seq(1, length(x), plot_resolution)]
  y = y[seq(1, length(y), plot_resolution)]

  lines(date_stamps, y,  
    col = line_color, lty=1, lwd = plot_line_width)

  legend_x_coord = median(date_stamps[1:round(length(date_stamps))])
  legend_location = list(x = legend_x_coord, y = 3)
  legend(legend_location, legend = c("Observed", "Hyper STL"),
   bg = "transparent", col = c('burlywood3', line_color), lwd=plot_line_width,
   cex = legend_axes_labels_font, bty = "n")
}

plotObserved <- function(data, depth, plot_resolution = 1, plot_line_width = 3, 
  legend_axes_labels_font = 1.3, line_color = 'burlywood3', 
  y_axis_label = expression('VWC (' ~ mm^3 * '/' * mm^3 * ')')) {
  
  y <- data[[depth]]
  x <- data$date
  date_stamps = x[seq(1, length(x), plot_resolution)]

  plot(date_stamps, y[seq(1, length(y), plot_resolution)], type = 'l', 
    xaxt = "n", col = line_color, xlab = "Date", 
    ylab = y_axis_label, 
    lwd = plot_line_width, cex.lab=legend_axes_labels_font)
  axis.POSIXct(1, at = seq(date_stamps[1], date_stamps[length(date_stamps)], length.out = 4),
   format = "%m/%d/%y", cex.axis = 1.2)
  grid(10, 10, lwd = 1.5)
}

plotSpline <- function(data, depth, plot_resolution = 1, plot_line_width = 3, 
  legend_axes_labels_font = 1.3, line_color = 'blue') {
  y <- data[[depth]]
  x <- data$date
  date_stamps = x[seq(1, length(x), plot_resolution)]

  y.spline <- smooth.spline(y, spar=0.5)
  lines(date_stamps, y.spline$y[seq(1, length(y.spline$y), plot_resolution)],  
    col = line_color, lty=1 , lwd = plot_line_width)

  legend_x_coord = median(date_stamps[1:round(length(date_stamps))])
  legend_location = list(x = legend_x_coord, y = 3)
  legend(legend_location, legend = c("Observed", "Spline"),
   bg = "transparent", col = c('burlywood3', line_color), lwd=plot_line_width,
   cex = legend_axes_labels_font, bty = "n")
}

plotSMA <- function(data, depth, plot_resolution = 1, plot_line_width = 3, 
  legend_axes_labels_font = 1.3, line_color = 'blue') {
  y <- data[[depth]]
  x <- data$date
  date_stamps = x[seq(1, length(x), plot_resolution)]

  y.sma <- SMA(y, n=50)
  lines(date_stamps, y.sma[seq(1, length(y.sma), plot_resolution)],
    col = line_color, lwd = plot_line_width)
  
  legend_x_coord = median(date_stamps[1:round(length(date_stamps)*0.5)])
  legend_location = list(x = legend_x_coord, y = 0.15)
  legend(legend_location, legend = c("Observed", "SMA"),
   bg = "transparent", col = c('burlywood3', line_color), lwd=plot_line_width,
   cex = legend_axes_labels_font)
}

plotWMA <- function(data, depth, plot_resolution = 1, plot_line_width = 3, 
  legend_axes_labels_font = 1.3, line_color = 'blue') {
  y <- data[[depth]]
  x <- data$date
  date_stamps = x[seq(1, length(x), plot_resolution)]

  w <- rep(1, length(y))
  w[y > mean(y)] <- 1000
  y.wma <- WMA(y, w, n=50)
  lines(date_stamps, as.numeric(y.wma[seq(1, length(y.wma), plot_resolution)]), 
    col = line_color, lwd = plot_line_width)

  legend_x_coord = median(date_stamps[1:round(length(date_stamps)*0.5)])
  legend_location = list(x = legend_x_coord, y = 0.15)
  legend(legend_location, legend = c("Observed", "WMA"),
   bg = "transparent", col = c('burlywood3', line_color), lwd=plot_line_width,
   cex = legend_axes_labels_font)
}

plotLOESS <- function(data, depth, plot_resolution = 1, plot_line_width = 3, 
  legend_axes_labels_font = 1.3, line_color = 'blue') {
  y <- data[[depth]]
  x <- data$date
  date_stamps = x[seq(1, length(x), plot_resolution)]

  y.loess <- loess(y ~ x, span = 0.15, 
                   data.frame(x=1:length(y),y=y))
  y.predict <- predict(y.loess)
  lines(date_stamps, y.predict[seq(1, length(y.predict), plot_resolution)],
   col = line_color, lwd = plot_line_width)

  legend_x_coord = median(date_stamps[1:round(length(date_stamps)*0.5)])
  legend_location = list(x = legend_x_coord, y = 0.15)
  legend(legend_location, legend = c("Observed", "LOESS"),
   bg = "transparent", col = c('burlywood3', line_color), lwd=plot_line_width,
   cex = legend_axes_labels_font)
}

plotPeakPreserving <- function(data, depth, plot_resolution = 1, plot_line_width = 3, 
  legend_axes_labels_font = 1.3, line_color = 'blue') {
  y <- data[[depth]]
  x <- data$date
  date_stamps = x[seq(1, length(x), plot_resolution)]

  y.hat <- peakPreserving(data, depth, 0.02) 
                   
  lines(date_stamps, y.hat[seq(1, length(y.hat), plot_resolution)],
   col = line_color, lwd = plot_line_width, lty = 2)

  legend_x_coord = median(date_stamps[1:round(length(date_stamps)*0.87)])
  legend_location = list(x = legend_x_coord, y = 2)
  legend(legend_location, legend = c("Observed", "Peak Preserving"),
   bg = "transparent", col = c('burlywood3', line_color), lwd=plot_line_width,
   lty = c(1,2), cex = legend_axes_labels_font, bty = "n")
}


peakPreserving <- function(data, depth, threshold) {
  helper <- function(y) {
    if(length(y) == 1) {
      return(y)
    }
    y.center <- getLine(y)
    y.lower <- y.center - threshold
    y.upper <- y.center + threshold
    
    idx <- which(y > y.upper)[which.max(y[which(y > y.upper)])]
    if(length(idx) == 0) {
      idx <- which(y < y.lower)[which.min(y[which(y < y.lower)])]  
    }
    
#     plot(y, type='l')
#     lines(y.center)
#     lines(y.lower, lty=3)
#     lines(y.upper, lty=3)
  
    if(length(idx) != 0) {
#       print(idx)
      y[1:idx] <- helper(y[1:idx])
      y[idx:length(y)] <- helper(y[idx:length(y)])
    } else {
      y <- y.center
    }
    return(y)
  }
  y <- data[[depth]]
  y.hat <- helper(y)
   # generate.plot.special.purpose(y, y.hat)
   y.hat
#
}

getLine <- function(y) {
  #mx + c
  left <- matrix(c(1, 1, length(y), 1), nrow=2, byrow=TRUE)
  #y
  right <- matrix(c(y[1],y[length(y)]), nrow=2, ncol=1)
  sol <- solve(left, right)
  y.line <- 1:length(y) * sol[1,1] + sol[2,1]
  y.line
}


generate.plot.special.purpose <- function(y, y.hat, image_scale = 2.8) {
    plot_line_widths = 3

dev.off()
setEPS()
postscript("comparison-with-peak-preserving.eps", width = 3 * image_scale, height = 2 * image_scale, paper="a4")
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 2, 0, 0))

    plot(y, type='l', col='blue',
    ylim=c(min(y[!is.na(y)], y.hat[!is.na(y.hat)]),
    max(y[!is.na(y)], y.hat[!is.na(y.hat)])),  xlab = "Time Index",
    ylab = expression('VWC (' ~ mm^3 * '/' * mm^3 * ')'), lwd = plot_line_widths, cex.lab=1.5)
    lines(y.hat, col='red', lty=1 , lwd = plot_line_widths)
    
    legend_location = "topright"
    legend(legend_location, legend=c("Observed","Peak presering smoothing"),
        lty=c(1,1),
         bg="transparent",
        col=c("blue", "red"),
        lwd=plot_line_widths, cex=1.5)
    dev.off()

}

setupEpsFileForPlot <- function(filename, image_scale = 2.5) {
  setEPS()
  postscript(filename, width = 2.4 * image_scale, 
    height = 2 * image_scale, paper="a4")
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 2, 0, 0))

  }

TestErrorAWI <- function (d2, depth, y.smooth, train.end, 
  num.points.in.day, subsample.rate = 10, shift.prediction.back = FALSE) {
  num.points.in.day = num.points.in.day 
  y.smooth[is.na(y.smooth)] <- 0
  train.data <- d2[1:(train.end / subsample.rate), ]
  train.data[[depth]] = y.smooth[1:(train.end / subsample.rate)]
  params = fit.mod.AWI.model.to.data(train.data, depth, tau = num.points.in.day)
  # print(params)
  test.data <- d2[(train.end / subsample.rate + 1):nrow(d2), ]
  test.data[[depth]] <- y.smooth[(train.end / subsample.rate + 1):length(y.smooth)]
  y.predict = mod.AWI.model.prediction(params, test.data, depth, tau = num.points.in.day)

  actual.test.data <- d2[(train.end / subsample.rate + 1):nrow(d2), ]
  y = actual.test.data[[depth]]
  
  if (shift.prediction.back) {
    n = length(y.predict)    
    y.predict.shifted = y.predict[(num.points.in.day+1):n]
    y = y[1:(n-num.points.in.day)]
    y.predict = y.predict.shifted
    # print(length(y))
    # print(length(y.predict))
    # print(which.max(abs(y- y.predict)))
  }
  errors = c(rmse(y, y.predict), max(abs(y- y.predict)))

  # plot(y,type='l',col='blue')
  # lines(y.predict,col='red')
  # print(errors)
  # char = readline(prompt="Press [enter] to continue")
  # stopifnot(char == '')

  errors
}










