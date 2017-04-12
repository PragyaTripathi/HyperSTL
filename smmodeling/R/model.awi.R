
#' @title calculate.awi
#' @name calculate.awi
#' @description \cr
#' \deqn{
#'    AWI_i = AWI_i + 1
#' }
#' @return data.frame
#' @export

# modified AWI model
fit.mod.AWI.model.to.data <- function(data, 
                                  depth, 
                                  tau  = 72,
                                  max.iter = 20,
                                  run.index = 0,
                                  return_parameters = FALSE) {
  fitness.to.min <- function(params) {
    kd = 10 ^ (params[1]/10)
    eta = params[3]
    y.hat = mod.awi.predict.regular(data, depth, kd, eta, tau)
    moisture = data[[depth]]
    mse <- mean(sqrt((moisture - y.hat) ^ 2))
    mse
  }
  
  MinParams = c(-60,1)
  MaxParams = c(0,4e3)

  outDEoptim <- DEoptim(fn = fitness.to.min, lower = MinParams, 
    upper = MaxParams, control = DEoptim.control(itermax = max.iter, 
    trace = FALSE))
  fitted.parameters = tail(outDEoptim$member$bestmemit,1)
  if (return_parameters) {
    return(fitted.parameters)
  }

  fitted.parameters
}
mod.awi.predict.regular <- function (data, depth, kd, eta, tau) {
  moisture = data[[depth]]
  rain = data$rainfall
  awi.predicted = moisture
  for(t in (tau+1):length(moisture)) {
    awi.predicted[[t]] = moisture[[t-tau]] * exp(- kd * tau)
    + rain[[t]] * (1 - exp(- kd * tau)) / eta
  }
  awi.predicted
}
mod.AWI.model.prediction <- function(params, data, depth, tau=72) {
  moisture <- data[[depth]]
  rain <- data$rainfall
  kd = 10 ^ params[1]
  eta = params[3]
  y.hat = mod.awi.predict.regular(data, depth, kd, eta, tau)

  y.hat
}


fit.existing.AWI.model.to.data <- function(data, 
                                  depth, 
                                  tau  = 72,
                                  max.iter = 20,
                                  run.index = 0,
                                  return_parameters = FALSE) {
  fitness.to.min <- function(params) {
    kd = 10 ^ params[1]
    y.hat = existing.awi.predict.regular(data, depth, kd, tau)
    moisture = data[[depth]]
    mse <- mean(sqrt((moisture - y.hat) ^ 2))
    mse
  }
  
  MinParams = c(-6)
  MaxParams = c(0)

  outDEoptim <- DEoptim(fn = fitness.to.min, lower = MinParams, 
    upper = MaxParams, control = DEoptim.control(itermax = max.iter, 
    trace = FALSE))
  fitted.parameters = tail(outDEoptim$member$bestmemit,1)
  if (return_parameters) {
    return(fitted.parameters)
  }

  fitted.parameters
}

existing.AWI.model.prediction <- function(params, data, depth, tau=72) {
  moisture <- data[[depth]]
  rain <- data$rainfall
  kd = 10 ^ params[1]
  y.hat = existing.awi.predict.regular(data, depth, kd, tau)

  y.hat
}
  
#' @export
fitness.awi.regular <- function (data, depth, kd, tau) {
  moisture = data[[depth]]
  predicted.awi = existing.awi.predict.regular(data, depth, kd, tau)
  mse <- mean(sqrt((moisture - predicted.awi) ^ 2))

  mse
  }

#' @title fitness.awi.irregular
#' @export
fitness.awi.irregular <- function (data, depth, kd, tau) {
  moisture = data[[depth]]
  predicted.awi = existing.awi.predict.irregular(data, depth, kd, tau)
  mse <- mean(sqrt((moisture - predicted.awi) ^ 2))

  mse
  }


#' @export
existing.awi.predict.regular <- function (data, depth, kd, tau) {
  moisture = data[[depth]]
  rain = data$rainfall
  awi.predicted = moisture
  for(t in (tau+1):length(moisture)) {
    awi.predicted[[t]] = moisture[[t-tau]] * exp(- kd * tau)
    + rain[[t]] * (1 - exp(- kd * tau)) / kd
  }
  awi.predicted
}

#' @export
existing.awi.predict.irregular <- function (data, depth, kd, tau) {
  moisture = data[[depth]]
  rain = data$rainfall
  awi.predicted = moisture
  for(t in (tau+1):length(moisture)) {
    awi.predicted[[t]] = awi.predicted[[t-tau]] * exp(- kd * tau)
    + rain[[t]] * (1 - exp(- kd * tau)) / kd
  }
  awi.predicted
}


