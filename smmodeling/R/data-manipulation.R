#' @title subsample.data
#' @export
subsample.data <- function(data, subsampling.rate = 10) {
  d = data[seq(1, nrow(data), subsampling.rate), ]
  rain = data$rainfall
  accum.rain = cumsum(rain)
  sub.accum.rain = accum.rain[seq(1, length(accum.rain), subsampling.rate)]
  l = length(sub.accum.rain)
  sub.rain = c(0, sub.accum.rain[-1] - sub.accum.rain[-l])
  d$rain = sub.rain

  d
}


# minparams=c(7, 7, 72),
# maxparams=c(99, 99, 200),
# Range = (92, 92, 80)
# Bits needed = (7, 7, 7)
decode.params.GA <- function(binary.params) {
  BAD_INPUT_RETURN = c(55, 55, 55)
  if (length(binary.params) %% 21 != 0) {
    print(cat("BAD INPUT: ",length(binary.params), " bits "))
    return(BAD_INPUT_RETURN)
  }
  x1 = 7 + unbinary(paste(binary.params[1:7], collapse=""))
  x2 = 7 + unbinary(paste(binary.params[8:14], collapse=""))
  x3 = 72 + unbinary(paste(binary.params[15:21], collapse=""))
  integer.parameters = c(x1, x2, x3)

  return(integer.parameters)
}


#' @title obj.smoothing
#' @description Removes all non-numeric values form data
#' @description period.series should be in minutes
#' @return data.frame
#' @export
obj.smoothing <- function(data, depth, parameters, 
                          subsample = 1, period.series = 2, 
                          weights.obj.terms = c(20, 10, 30)) { 
  input.series = data[[depth]]
  y.hat.no.rain <- poly.smooth.no.rain(data, depth)

  #  The "frequency" is the number of observations per season.
  #  Season = 24 hours
  #  For hourly data: frequency = 24
  y.ts <- ts(input.series, frequency = 24 * 60 / period.series)

  y.stl <- stl(y.ts, s.window = parameters[1], t.window = parameters[2],
              l.window = parameters[3], robust=TRUE)

  r <- as.numeric(y.stl$time.series[,"remainder"])
  y.trend <- as.numeric(y.stl$time.series[,"trend"])
  mse <- sqrt(mean((y.trend - y.hat.no.rain) ^ 2))
  objective <- sum(c(sd(r), abs(max(r) - min(r)), mse) * weights.obj.terms)

  objective
}

#' @title stl.by.params
#' @description Removes all non-numeric values form data
#' @return stl object
#' @export
stl.by.params <- function(data, depth, params) {
  y <- data[[depth]]
  time.diff <- difftime(data$date[2], data$date[1], units="mins") 
  time.diff <- as.numeric(time.diff, units="mins")
  
  y.ts <- ts(y, frequency=24 * 60 / time.diff)

    y.stl <- stl(y.ts, s.window=params[1],
                 t.window=params[2],
                 l.window=params[3],
                 robust=TRUE)
    y.stl
  }
  


poly.smooth.no.rain <- function(data, depth) {
    y <- data[[depth]]
    intervals <- findInterval(data$rainfall, c(0,1))
    index <- 1
    while(index < length(intervals)) {
      startIndex <- index
      
      while(index < length(intervals) 
            && intervals[index] != 2) {
        index <- index + 1
      }
      
      ysplit <- y[startIndex:index]
      chunk <- data.frame(y=ysplit, x=1:length(ysplit))
      
      model <- lm(y ~ stats:::poly(x, 1, raw=TRUE), 
                  data=chunk)
      
      y[startIndex:index] <- fitted(model)
      index <- index + 1
      
    }
    y
  }


#' @title replace.negative
#' @description Removes all non-numeric values form data
#' @return data.frame
#' @export
clean.data <- function(data, depths) {
    for(i in depths) {
        y <- data[[i]]
        # ---- replace wired values ----
        y2 = as.numeric(y)
        y2[is.na(y2)] = 0.001
        # -------
        data[[i]] <- y2
    }
    data
}

#' @title replace.negative
#' @description Replaces negative values with 0
#' @return data.frame
#' @export
replace.negative <- function(data, depths) {
  for(i in depths) {
    y <- data[[i]]
    y[y < 0] <- 0
    data[[i]] <- y
  }
  data
}

#' @title drop.na
#' @description Removes NA values from different depths
#' @return data.frame
#' @export
drop.na <- function(data, depths) {
  for(i in depths) {
    data <- data[!is.na(data[[i]]), ]
  }
  data
}

#' @title smooth.stl.remainder
#' @name smooth.stl.remainder
#' @description \cr
#' Uses STL remainder threshold approach to 
#' calculate smoothed signal. 
#' TODO: Add equation
#' Make sure data doesn't have na \cr
#' \link{stl}
#' @return smoothed data vector
#' @export
smooth.stl.remainder <- function(data, depth) {
  y <- data[[depth]]
  ysplit <- y
  
  time.diff <- difftime(data$date[2], data$date[1], 
                        units="mins")
  
  time.diff <- as.numeric(time.diff, units="mins")
  
  ysplit.ts <- ts(ysplit, frequency=24 * 60 / time.diff)
  y.stl <- stl(ysplit.ts, s.window='per')
  
  plot(y.stl)
  
  rem <- y.stl$time.series[,"remainder"]
  threshold <- mean(rem) + 2 * sd(rem)
  above <- which(rem > threshold)
  y.new <- y.stl$time.series[,"trend"]
  y.new[above] <- y[above]
    
  data.sm <- data
  data.sm[[depth]] <-as.numeric(y.new)
  data.sm
}

#' @title optimize.stl
#' @description Optimizes STL with GA using #{equation} 
#' optimization function
#' @import GA
#' @export
#' @param data Soil moisture and rainfall data
#' @param depth Soil moisture depth to smooth
#' @param weights Weights for objective function
#' TODO: add equation
#' @return STL object
# optimize.stl <- function(data, depth, 
#                          weights=c(20, 10, 30),
#                          minparams=c(7, 7, 721),
#                          maxparams=c(99, 99, 801),...) {
optimize.stl <- function(data, depth, 
  algorithm = 'sa',
  max.feval = 2500,
  weights=c(20, 10, 30),
  minparams=c(7, 7, 72),
  maxparams=c(99, 99, 200),
  data.freq = 24,...) {  
  print("In optimize.stl")

  MinParams <- minparams
  MaxParams <- maxparams
  iters <- 50
  weights = weights / sum(weights)
  time.diff <- difftime(data$date[2], data$date[1], 
                        units="mins")
  
  time.diff <- as.numeric(time.diff, units="mins")
  
  poly.smooth <- function(data, depth) {
    y <- data[[depth]]
    intervals <- findInterval(data$rainfall, c(0,1))
    index <- 1
    while(index < length(intervals)) {
      startIndex <- index
      
      while(index < length(intervals) 
            && intervals[index] != 2) {
        index <- index + 1
      }
      
      ysplit <- y[startIndex:index]
      chunk <- data.frame(y=ysplit, x=1:length(ysplit))
      
      model <- lm(y ~ stats:::poly(x, 1, raw=TRUE), 
                  data=chunk)
      
      y[startIndex:index] <- fitted(model)
      index <- index + 1
      
    }
    y
  }
  
  y <- data[[depth]]
  
  #  The "frequency" is the number of observations per season.
  #  Season = 24 hours
  #  For hourly data: frequency = 24
  y.ts <- ts(y, frequency=data.freq)
  yhat <- poly.smooth(data, depth)
      
  # pb <- txtProgressBar(min = 0, max = iters
  #                      , style = 3)
  
  getStl <- function(params) {
    y.stl <- stl(y.ts, s.window=params[1],
                 t.window=params[2],
                 l.window=params[3],
                 robust=TRUE)
    y.stl
  }
  
  evalFunc <- function(params) {
    -evalFuncToMinimize(params)
  }

  evalFuncToMinimize <- function(params) {
    y.stl <- getStl(params)
    r <- as.numeric(y.stl$time.series[,"remainder"])
    y.trend <- as.numeric(y.stl$time.series[,"trend"])
    mse <- sqrt(mean((y.trend - yhat) ^ 2))

    # print(weights)
    objective <- sum(c(sd(r), abs(max(r) - min(r)), mse) * weights)
    if (is.nan(objective) || is.na(objective) || is.null(objective)) {
      objective = 99999999
    }
    # print(objective)
    objective
  }
  
  evalFuncForBinary <- function(binary.params) {
     decoded.params =  decode.params.GA(binary.params)
     # penalty for out of range input 
     OUT_OF_RANGE_PENALTY = 1e2;

     params = decoded.params 
     # # decoding min is set to MinParams
     # ind.small = which(decoded.params < MinParams)
     # params[ind.small] = MinParams[ind.small]
     ind.Big = which(decoded.params > MaxParams)
     params[ind.Big] = MaxParams[ind.Big]

    y.stl <- getStl(params)
    r <- as.numeric(y.stl$time.series[,"remainder"])
    y.trend <- as.numeric(y.stl$time.series[,"trend"])
    mse <- sqrt(mean((y.trend - yhat) ^ 2))
    objective <- sum(c(sd(r), abs(max(r) - min(r)), mse) * weights)

    objective <- objective 
     + OUT_OF_RANGE_PENALTY * norm(params - decoded.params, type="2")^2

    -objective
  }

  # monitorFunc <- function(object) {
  #   setTxtProgressBar(pb, object@iter)
  # }
  # ------------------------------------------------------------------
  # algorithm = 'irace'
  switch(algorithm,
    sa = {
      # Simulated Annealing
      out <- GenSA(lower = MinParams, upper = MaxParams, fn = evalFuncToMinimize,
               control=list(max.call= max.feval, verbose=TRUE))
      print(out$par)
      print(cat("Evalfunc(GenSA) = ", evalFuncToMinimize(out$par)))
      solution = out$par    
      solution.obj = evalFuncToMinimize(solution)
      },
  # ------------------------------------------------------------------      
    de = {
      # DE/rand/1/bin
      itermax = max(1,floor(max.feval/50))
      outDEoptim <- DEoptim(evalFuncToMinimize, MinParams, MaxParams, 
        DEoptim.control(NP = 50, itermax = itermax, F = 0.8))
      solution = outDEoptim$optim$bestmem
      solution.obj = evalFuncToMinimize(solution)
      },
  # ------------------------------------------------------------------
    jade = {
      # JADE
      itermax = max(1,floor(max.feval/50))
      outDEoptim <- DEoptim(evalFuncToMinimize, MinParams, MaxParams, 
        DEoptim.control(NP = 50, itermax = itermax, F = 0.8, CR = 0.5, strategy = 6, c = 0.4))
      solution = outDEoptim$optim$bestmem
      solution.obj = evalFuncToMinimize(solution)
        # print(evalFuncToMinimize(solution))
      },
  # ------------------------------------------------------------------
    sga = {
      # Real coded GA
      GAModel <- ga(type = "real-valued",
                 fitness = evalFunc,
                 monitor=NULL,
                 min = MinParams,
                 max = MaxParams,
                 maxfitness=max.feval)
      solution = GAModel@solution
      solution.obj = -evalFunc(GAModel@solution)
    
    },
  # ------------------------------------------------------------------
    bga = {
      # Binary-GA
      GAModel <- ga(
        type = "binary", 
        fitness = evalFuncForBinary, 
        nBits = 21,
        maxfitness=max.feval,
        monitor = NULL)
      solution.binary = GAModel@solution
      solution = decode.params.GA(solution.binary)
      ind.Big = which(solution > MaxParams)
      solution[ind.Big] = MaxParams[ind.Big]

      solution.obj = -1 * evalFuncForBinary(solution.binary)
      # print(GAModel@solution)
      
      # print(solution)
      print(-1 * evalFuncForBinary(solution.binary))

      # print(obj.smoothing(data, depth, solution))

      },
  # ------------------------------------------------------------------
    irace = {
      # iRace
      hook.run = function(experiment, config = list()) {
        candidate <- experiment$candidate
        x = c(candidate[["x1"]], candidate[["x2"]], candidate[["x3"]])
        # print(x)
        y = evalFuncToMinimize(x)
        y
      }
      # parameters.table <- 'x1 \"\" i \n (7, 99) x2 \"\" i (7, 99) \n x3 \"\" i (72, 200)'
      parameters.table <- 'x1 "" i (7, 99)\n x2 "" i (7, 99)\n x3 "" i (72, 200)' 
      parameters <- readParameters(text=parameters.table)
      print('Read parameters')
      irace.weights <- rnorm(200, mean = 0.9, sd = 0.02)
      result <- irace(tunerConfig = list(
        hookRun = hook.run,
        instances = irace.weights[1:100],
        maxExperiments = 4000,
        # nbIterations = 1,
        logFile = ""),
        parameters = parameters)

      print(cat("\n\n\n final result  "))
      solution = c(result[1,]$x1, result[1,]$x2, result[1,]$x3)
      solution.obj = evalFuncToMinimize(solution)
      print(evalFuncToMinimize(solution))
      },
  # ------------------------------------------------------------------
    {
      print('No algorithm specified')
      solution = c(0,0,0)
      solution.obj = 12345.12345
    }
      )
  stl.result <- getStl(solution)
  
  # print(solution.obj)
  # solution.obj

  stl.result
}

#' @title generate.smoothed
#' @description creates smoothed data from 
#' seasonality threshold
#' @param y.stl stl class obtained from optimization
#' @export
generate.smoothed <- function(y.stl) {
  mmin <- min(as.numeric(y.stl$time.series[,"seasonal"]))
  mmax <- max(as.numeric(y.stl$time.series[,"seasonal"]))
  y <- as.numeric(y.stl$time.series[,"trend"])
  r <- as.numeric(y.stl$time.series[,"remainder"])
  y[which(r > mmax)] <- y[which(r > mmax)] + r[r > mmax]
  y[which(r < mmin)] <- y[which(r < mmin)] + r[r < mmin]
  y
}
