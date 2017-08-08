# Apply HyperSTL and plot result
source('HyperSTL.R')


#bucket.data = read.csv("InfraWatch12DaysMean.csv", stringsAsFactors=FALSE)
bucket.data = read.csv("InfraWatch12DaysMedian.csv", stringsAsFactors=FALSE)
bucket.data$date <- as.POSIXct(bucket.data$Date, format = "%m/%d/%Y %H")
source('PeakCreator.R')
smooth_trained_data = HyperSTL(bucket.data, 'Median.Strain', data.freq = 24)
plot(smooth_trained_data)
#fit1 <- arima(smooth_trained_data$time.series[, "trend"], order = c(3,0,0))
#pred<-predict(fit1, 100)
#plot(pred$pred)


