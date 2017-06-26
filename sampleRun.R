# Read data
bucket.data = read.csv("experimental-data.csv", stringsAsFactors=FALSE)
bucket.data$date <- as.POSIXct(bucket.data$Measurement.Time, format = "%m/%d/%Y %I:%M %p")
d2 = bucket.data[-1,]
colnames(d2)[which(names(d2) == "PrecipitationIn")] <- "rainfall"

# Apply HyperSTL and plot result
source('HyperSTL.R')
y=HyperSTL(d2,'D..10cm')
plot(y)
