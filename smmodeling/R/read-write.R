#' @title read.soilmoisture
#' @name read.soilmoisture
#' @description reads soil moisture data from csv file, \cr \cr
#' \strong{Expected format} \cr
#' \code{date, depth1, depth2, depth3, ..., depthN}
#' \code{date should be in yyyy/mm/dd h:mm:ss format}
#' @param filepath file to read
#' @return data.frame
#' @export
read.soilmoisture <- function(filepath) {
  data <- read.csv(filepath, stringsAsFactors=FALSE)
  data$date <- as.POSIXct(data$date, format="%m/%d/%Y %H:%M")
  data
}

#' @title read.rainfall
#' @name read.rainfall
#' @description reads soil moisture data from csv file, \cr \cr
#' expected format: \cr
#' \code{date, depth1, depth2, depth3, ..., depthN}
#' @param filepath file to read
#' @return data frame
#' @export
read.rainfall <- function(filepath) {
  data <- read.csv(filePath, stringsAsFactors=FALSE)
  data$date <- as.POSIXct(data$date)
  data
}

#' @export
read.ofs <- function(filepath) {
  data <- read.csv(filePath, stringsAsFactors=FALSE)
  data$date <- as.POSIXct(data$date)
  data
}

#' @title align.rainfall
#' @description Aligns rainfall and soil moisture
#' readings based on timestamp
#' @return data.frame
#' @export
align.rainfall <- function(soilmoisture, rainfall) {
  soilm = soilmoisture
  soilm$rainfall = rep(0, nrow(soilm))
  pb <- txtProgressBar(min = 0, max = nrow(rainfall)
                       , style = 3)

  for(i in 1:nrow(rainfall)) {
    row = rainfall[i,]
    index = which(soilm$date > row$date)[1]

    setTxtProgressBar(pb, i)

    if(is.na(index)) {
      #print(i)
      next
    }
    if(index != 0) {
      # soilm[index,]$rainfall = row$mm.hr
      soilm[index,]$rainfall = row$rainfall
    }
  }
  soilm
}

#' @title write.data
#' @description Writes soil moisture data in csv
#' @export
write.data <- function(data, filename) {
  write.csv(data, filename, row.names=FALSE)
}
