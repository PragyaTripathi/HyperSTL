# Read data
bucket.data = read.csv("InfraWatch_S100_15-11-08_26-11-08_10Hz.csv", stringsAsFactors=FALSE)
bucket.data$date <- as.POSIXct(bucket.data$Timestamp, format = "%Y-%m-%d %H:%M:%S")

# Apply HyperSTL and plot result
source('HyperSTL.R')
y=HyperSTL(bucket.data, 'Strain', data.freq = 14400)
plot(y)