df = read.csv("InfraWatch_S100_15-11-08_26-11-08_10Hz.csv", stringsAsFactors=FALSE)

dates<-c()
strainMeasurementsForAnHour<-c()
finalMatrix<-matrix(NA, nrow = (12*24), ncol = 36000)
jumpedDateName = FALSE
count = nrow(df)
for ( i in 1:count) {
  dateName <- sub(":.*", "", df$Timestamp[i])
  if (is.element(dateName, dates)) {
    strainMeasurementsForAnHour<-c(strainMeasurementsForAnHour, df$Strain[i])
    jumpedDateName = FALSE
  } else {
    dates<-c(dates, dateName)
    jumpedDateName = TRUE
    print(length(dates))
    print(i)
  }
  
  if (jumpedDateName || (i == count)) {
    if (i == count) { matrixIndex<-length(dates) } else { matrixIndex<-(length(dates)-1) }
    naNeeded<-(36000 - length(strainMeasurementsForAnHour))
    if (naNeeded > 0) {
      strainMeasurementsForAnHour<-c(strainMeasurementsForAnHour,rep(c("NA"),naNeeded))
    }
    finalMatrix[matrixIndex,]<-strainMeasurementsForAnHour
    strainMeasurementsForAnHour<-c(df$Strain[i])
  }
}
