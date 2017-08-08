source('findpeaks.R')

peaks <- findpeaks(bucket.data$Mean.Strain, nups=4)

rainfall<-c(rep(0,288))
for (i in 1:nrow(peaks)) {
  rainfall[peaks[i,2]]<-peaks[i,1]
}

bucket.data$rainfall<-rainfall