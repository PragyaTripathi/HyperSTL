HyperSTL
---

Accompanying code for paper "Optimizing the Decomposition of Time Series using Evolutionary Algorithms: Soil Moisture Analytics". Contains HyperSTL function to smooth the signal while preserving peaks and valleys.

# HyperSTL

## Dependencies

Luca Scrucca (2013). GA: A Package for Genetic Algorithms in R.
  Journal of Statistical Software, 53(4), 1-37. URL
  http://www.jstatsoft.org/v53/i04/.

Mullen, K.M, Ardia, D., Gil, D., Windover, D., Cline, J. (2011). DEoptim: An R Package for
    Global Optimization by Differential Evolution. Journal of Statistical Software, 40(6), 1-26. URL
    http://www.jstatsoft.org/v40/i06/.

Sylvain Gubian, Yang Xiang, Brian Suomela, Julia Hoeng, PMP SA (2016). GenSA: R Functions for Generalized Simulated Annealing
    https://cran.r-project.org/web/packages/GenSA/GenSA.pdf

## Usage

```
y.optimized.stl = HyperSTL(data, "value", 
                                algorithm = 'sga', 
                                max.feval = 1e2,
                                weights = c(5, 5, 8), 
                                data.freq = num.points.in.day, 
                                subsample.rate)

Params:
    #' @param data data.frame containing time series
    #' @param column column containing values of time series
    #' @param algorithm Optimization algorithm. Possible values are 
    #' 'sa', 'de', 'jade', 'sga', 'bga', 'irace'. Defaults to 'sa'.
    #' @param max.feval Max number of function evaluations. Defults to 2500.
    #' @param weights Weights for decomposed time series using STL. 
    #' Expects a vector containing three values. e.g. c(10,10,10) i.e. Equal weights for trend, seasonality and remainder.
    #' @param minparams Minimum params for the grid search
    #' @param maxparams Maximum params for the grid search
    #' @param data.freq 
    #' @param ... Rest of the params are passed to STL function

Returns:
    STL object (object containing trend, seasonality, and remainder).
```
