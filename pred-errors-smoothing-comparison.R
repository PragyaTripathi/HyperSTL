# setwd("~/Google Drive/USGS-project/code/")
# source("pred-errors-smoothing-comparison.R")	

# detach(name="package:smmodeling", unload=TRUE)
# install.packages("./smmodeling/", repos = NULL, type="source")
# library(smmodeling)


# library(smmodeling)
library(DEoptim)
library(GenSA)
library(Metrics)
# source('comparisons-v2.R');
source('comparisons-all-separate.R');
# dataset = 'infra'

# dataset = 'canyon'
dataset = 'bucket'
# dataset = 'gap'

switch(dataset,
    canyon = {
		d =  read.csv("~/Google Drive/USGS-project/code/data/prepared-csv/pit1.sanitized.csv",
		 stringsAsFactors=FALSE)
		d$date <- as.POSIXct(d$date, format = "%Y-%m-%d %H:%M:%S")
		d.smooth =  read.csv("~/Google Drive/USGS-project/code/data/prepared-csv/pit1.smooth.csv",
		 stringsAsFactors=FALSE)
		d.smooth$date <- as.POSIXct(d.smooth$date, format = "%Y-%m-%d %H:%M:%S")

		# 1e4, 3.5e4, 8e4
		subsample.rate = 10
		d2 = d[seq(1e4, 8e4, subsample.rate), ]
		d2.smooth = d.smooth[seq(1e4, 8e4, subsample.rate), ]
		# data has 1 points every 2 minutes => 0.5 of 24 hours in minutes
		num.points.in.day = 24*60*0.5/subsample.rate
		train.end = 3.5e4
		depthsArray = c("X5cm","X15cm","X30cm")
	},
	# -----------------------------------
	bucket = {
		bucket.data = read.csv(paste("~/Google Drive/USGS-project/code/data/",
			"lab-experiment/data-analysis/non-leaky-bucket-wunderground-data.csv", sep = ""), 
			stringsAsFactors=FALSE)
		bucket.data$date <- as.POSIXct(bucket.data$Measurement.Time, format = "%m/%d/%Y %I:%M %p")
		d2 = bucket.data[-1,]
		colnames(d2)[which(names(d2) == "PrecipitationIn")] <- "rainfall"
		subsample.rate = 1
		# data has 1 point every hour => 24 hours in day
		num.points.in.day = 24/subsample.rate
		train.end = 400
		depthsArray = c("D..10cm","D..20cm","D..28cm")
		# overlaySelectedMethods(bucket.data, "D..28cm", smooth.stl = smooth.stl, legend_y_coord = 0.3);
	},
	# -----------------------------------
	gap = {
		gap.data = read.csv(paste("~/Google Drive/USGS-project/code/data/",
			"gap-fire/aligned-gap-fire-data.csv", sep = ""), 
			stringsAsFactors=FALSE)
		gap.data$date <- as.POSIXct(gap.data$date, format = "%Y-%m-%d %H:%M:%S")
		
		subsample.rate = 10
		d2 = gap.data[seq(1, 6.3e4, subsample.rate), ]
		# data has 1 points every 2 minutes => 0.5 of 24 hours in minutes
		num.points.in.day = 24*60/(3*subsample.rate)
		train.end = 3.2e4
		depthsArray = c("X5cm","X15cm","X30cm")
		# overlaySelectedMethods(bucket.data, "D..28cm", smooth.stl = smooth.stl, legend_y_coord = 0.3);
	},
	# -----------------------------------
	infra = {
		d =  read.csv("~/Google Drive/USGS-project/code/data/prepared-csv/InfraWatch.csv",
			stringsAsFactors=FALSE)
		subsample.rate = 1000
		d2 = d[seq(1, nrow(d), subsample.rate), ]
		d2$date <- as.POSIXct(d2$timestamp, format = "%Y-%m-%d %H:%M:%S")
		# data has 10 points per second
		num.points.in.day = 24*3600*10/subsample.rate
		y.optimized.stl = optimize.stl(d2, "value", 
			algorithm = 'sga', max.feval = 1e2,
			weights = c(5, 5, 8), data.freq = num.points.in.day, subsample.rate)
		smooth.stl = as.numeric(y.optimized.stl$time.series[,"trend"])
		overlaySelectedMethods(d2, "value", smooth.stl = smooth.stl,legend_x_coord = 6000,
		 legend_y_coord = 7, plot_resolution = 10, y_axis_label = 'Strain',
		 filename = "comparisons-few-infra.eps");

		subsample.rate = 1000
		# data has 10 points per second
		num.points.in.day = 24*3600*10/subsample.rate


	},
	# -----------------------------------
	{
		print('Incorret dataset')
	}
	)


# d2$date <- as.POSIXct(d2$timestamp, format = "%Y-%m-%d %H:%M:%S")

	rmseMatrix = matrix(0, 3, 6)
	maxAbsErrMatrix = matrix(0, 3, 6)
for (i.depth in 1:3){
	i.depth
	depth = depthsArray[i.depth]

	y <- d2[[depth]]
	y.spline <- smooth.spline(y, spar=0.5)
	y.sma <- SMA(y, n=50)
	w <- rep(1, length(y))
	w[y > mean(y)] <- 1000
	y.wma <- WMA(y, w, n=50)
	y.loess.fit <- loess(y ~ x, span = 0.15, data.frame(x=1:length(y),y=y))
	y.loess <- predict(y.loess.fit)
	y.peakPreserve <- peakPreserving(d2, depth, 0.02) 
	if (dataset == 'canyon') {
		y.hyperSTL <- d2.smooth[[depth]]
	} else {
		# optimized.stl = optimize.stl(d2, depth, weights = c(2, 10, 50),
		# 	algorithm = 'de', max.feval = 1500)
		# y.hyperSTL <- as.numeric(optimized.stl$time.series[,"trend"])		
	}
	# train.data <- d2[1:(train.end / subsample.rate), ]
	# train.data[[depth]] = y.spline$y[1:(train.end / subsample.rate)]
	# params = fit.AWI.model.to.data(train.data, depth, tau = num.points.in.day, subsample.rate)
	# test.data <- d2[(train.end / subsample.rate + 1):nrow(d2), ]
	# y.spline.predict = AWI.model.prediction(params, test.data, depth, tau = num.points.in.day, subsample.rate)
	# y = test.data[[depth]]
	# cat('RMSE = ', rmse(y, y.spline.predict), ', Max. abs. error = ', max(abs(y- y.spline.predict)))

	print('Calling AWI: None')
	e.none <- TestErrorAWI(d2, depth, y, train.end, num.points.in.day, subsample.rate)
	print(e.none)
	print('Calling AWI: sma')
	print(e.sma)
	# e.sma <- TestErrorAWI(d2, depth, y.sma, train.end, num.points.in.day, subsample.rate)
	# print('Calling AWI: wma')
	# e.wma <- TestErrorAWI(d2, depth, y.wma, train.end, num.points.in.day, subsample.rate)
	# print('Calling AWI: loess')
	# e.loess <- TestErrorAWI(d2, depth, y.loess, train.end, num.points.in.day, subsample.rate)
	# print('Calling AWI: spline')
	# e.spline <- TestErrorAWI(d2, depth, y.spline$y, train.end, num.points.in.day, subsample.rate)
	# print('Calling AWI: peakPreserve')
	# e.peakPreserve <- TestErrorAWI(d2, depth, y.peakPreserve, train.end, num.points.in.day, subsample.rate)
	# print('Calling AWI: HyperSTL')
	# e.hyperSTL <- TestErrorAWI(d2, depth, y.hyperSTL, train.end, num.points.in.day, subsample.rate, TRUE)

	# rmseMatrix[i.depth,1] = e.sma[1] 
	# rmseMatrix[i.depth,2] = e.wma[1] 
	# rmseMatrix[i.depth,3] = e.loess[1] 
	# rmseMatrix[i.depth,4] = e.spline[1] 
	# rmseMatrix[i.depth,5] = e.peakPreserve[1] 
	# # rmseMatrix[i.depth,6] = e.hyperSTL[1] 
	# maxAbsErrMatrix[i.depth,1] = e.sma[2] 
	# maxAbsErrMatrix[i.depth,2] = e.wma[2] 
	# maxAbsErrMatrix[i.depth,3] = e.loess[2] 
	# maxAbsErrMatrix[i.depth,4] = e.spline[2] 
	# maxAbsErrMatrix[i.depth,5] = e.peakPreserve[2] 
	# maxAbsErrMatrix[i.depth,6] = e.hyperSTL[2] 

}
stop()
df = data.frame(rmseMatrix)
colnames(df) = c('SMA', 'WMA', 'LOESS', 'Spline', 'Peak Preserving', 'HyperSTL')
print(xtable(df,digits = 3,caption='Standard Error'))

df = data.frame(maxAbsErrMatrix)
colnames(df) = c('SMA', 'WMA', 'LOESS', 'Spline', 'Peak Preserving', 'HyperSTL')
print(xtable(df,digits = 3,caption='Maximum Absolute Error'))
# smooth.hyper.stl = as.numeric(y.optimized.stl$time.series[,"trend"])






