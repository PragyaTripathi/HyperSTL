# setwd("~/Google Drive/USGS-project/code/")
# source("compare-smooth-non-smooth-train.R")

detach(name="package:smmodeling", unload=TRUE)
install.packages("./smmodeling/", repos = NULL, type="source")
library(smmodeling)

library(DEoptim)
library(gplots)

max.iter.per.run = 10
subsample.period = 10

data.directory = "~/Google Drive/USGS-project/code/data/prepared-csv/";

d = read.csv(paste(data.directory, "pit1.sanitized.csv", sep = ""),
 		stringsAsFactors=FALSE)
d.smooth =  read.csv(paste(data.directory, "pit1.smooth.csv", sep = ""),
 	stringsAsFactors=FALSE) 		
d$date <- as.POSIXct(d$date, format = "%Y-%m-%d %H:%M:%S")
d.smooth$date <- as.POSIXct(d.smooth$date, format = "%Y-%m-%d %H:%M:%S")


d2 = d[6e3:3.6e4,]
d2.s = subsample.data(d2, subsample.period)

d2.smooth = d.smooth[6e3:3.6e4,]
d2.smooth.s = subsample.data(d2.smooth, subsample.period)

params = fit.AEAR.model.to.data(data = d2.s, depth = "X5cm", 
	tau = 720/subsample.period, max.iter = max.iter.per.run, 0)

params.smoooth = fit.AEAR.model.to.data(data = d2.smooth.s, depth = "X5cm", 
	tau = 720/subsample.period, max.iter = max.iter.per.run, 0)



d3 = d[3.6e4:6e4,]
d3.s = subsample.data(d3, subsample.period)
d3.smooth = d.smooth[3.6e4:6e4,]
d3.smooth.s = subsample.data(d3.smooth, subsample.period)

# p1 = c(-1.759354, 2.36687, 26.51044, 664.5325)
# p2 = c(-2.001128, 1.911353, 27.20462, 496.6171)
p1= params
p2 = params.smoooth

y1 = AEAR.model.prediction(p1, data = d3.s, depth = "X5cm", tau = 720/subsample.period)
y2 = AEAR.model.prediction(p2, data = d3.smooth.s, depth = "X5cm", tau = 720/subsample.period)


plot(d3.s$X5cm, type = 'l')
lines(y1, col = 'red')

plot(d3.smooth.s$X5cm, type = 'l')
lines(y2, col = 'red')

plot(d3.s$X5cm[150:1190], type = 'l')
lines(y1[150:1190], col = 'red')
plot(d3.smooth.s$X5cm[150:1190], type = 'l')
lines(y2[150:1190], col = 'red')



# ----------------------------------------
x <- d3.s$date[150:1190]
date_stamps = x[seq(1, length(x), 1)]
y_axis_label = expression('VWC (' ~ mm^3 * '/' * mm^3 * ')')

setEPS()
image_scale=2
  postscript("non-smmoth-AEAR.eps", width = 2.4 * image_scale, 
    height = 2 * image_scale, paper="a4")
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 2, 0, 0))

y = d3.s$X5cm[150:1190]
line_color1 = 'lightcoral'
line_color1 = 'lightsalmon3'

plot_line_width = 4
legend_axes_labels_font = 1.3
plot(date_stamps, y, type = 'l', 
    xaxt = "n", col = line_color1, xlab = "Date", 
    ylab = y_axis_label, 
    lwd = plot_line_width, cex.lab=legend_axes_labels_font)
  axis.POSIXct(1, at = seq(date_stamps[1], date_stamps[length(date_stamps)], length.out = 4),
   format = "%m/%d/%y", cex.axis = 1.2)
  grid(10, 10, lwd = 1.5)

y = y1[150:1190]
line_color2 = 'blue'
lines(date_stamps, y,  col = line_color2, lty=3 , lwd = plot_line_width)

legend_location = "bottomleft"
    legend(legend_location, legend=c("Observed","Predicted"),
        lty=c(1,3),
         bg="transparent",
        col=c(line_color1, line_color2),
        lwd=plot_line_width, cex=1.2, bty = "n")
dev.off()


setEPS()
image_scale=2
  postscript("smmoth-AEAR.eps", width = 2.4 * image_scale, 
    height = 2 * image_scale, paper="a4")
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 2, 0, 0))


y = d3.smooth.s$X5cm[150:1190]
line_color = 'lightcoral'
line_color = 'lightsalmon3'

plot_line_width = 4
legend_axes_labels_font = 1.3
plot(date_stamps, y, type = 'l', 
    xaxt = "n", col = line_color, xlab = "Date", 
    ylab = y_axis_label, 
    lwd = plot_line_width, cex.lab=legend_axes_labels_font)
  axis.POSIXct(1, at = seq(date_stamps[1], date_stamps[length(date_stamps)], length.out = 4),
   format = "%m/%d/%y", cex.axis = 1.2)
  grid(10, 10, lwd = 1.5)

y = y2[150:1190]
line_color = 'blue'
lines(date_stamps, y,  col = line_color, lty=3 , lwd = plot_line_width)
legend_location = "bottomleft"
    legend(legend_location, legend=c("Observed","Predicted"),
        lty=c(1,3),
         bg="transparent",
        col=c(line_color1, line_color2),
        lwd=plot_line_width, cex=1.2, bty = "n")
dev.off()












