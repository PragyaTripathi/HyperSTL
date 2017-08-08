finalMatrix<-as.matrix(read.csv('FinalMatrix.csv',sep=",", row.names=1))

squishMatrixRow <- function(row, type) {
  cleanRow <- na.omit(as.numeric(row))
  switch(type, mean = mean(cleanRow), median = median(cleanRow))
}

meanForRows <- c()
for (i in 1:nrow(finalMatrix)) {
  meanForRows<-c(meanForRows, squishMatrixRow(finalMatrix[i,], "median"))
}

newDF = data.frame(rownames(finalMatrix),meanForRows)

write.csv(newDF, 'InfraWatch12DaysMedian.csv')