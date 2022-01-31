library(foreach)
library(survival)
library(tidyverse)
library(ggplot2)
library(ipw)

#data simulator
#helper
lagger <- function(x) {
  if(length(x) == 1) {return(0)}
  c(0, x[1:(length(x) - 1)])
}

#simulate a subject's trajectory
gen_subject <- function(max_t = 5, id = NA, beta = 1, beta2 = 1, beta_xa = 0.2, tx = 0) {
  
  # The fugutive: the subject's hazard rate 
  # (spans factor of 3 across population)
  u <- .5 / max_t + 1 / max_t * runif(1)
  
  #generate the first time step
  t <- 0
  x <- beta * max_t * (u - 1 / max_t) + rnorm(1)
  a <- runif(1) < (.3 + .4 * (x < 0))
  
  #step forward in time
  out <- c()
  while(runif(1) > (u - u*tx*(a == 1))) {
    out <- rbind(out, c(t, t+1, x, a, 0))
    t <- t + 1
    x <- beta * max_t * (u - 1/max_t) + beta2 * (a == 1) + rnorm(1) # x reflects health status and Tx
    a <- runif(1) < .5 - beta_xa * x #if x small, higher prob of treatment
  }
  out <- rbind(out, c(t, t+1, x, a, 1))
  
  #cleanup
  out <- cbind(out, cumsum(out[, 3]) / 1:nrow(out), 
               cumsum(out[, 4]) / 1:nrow(out),
               lagger(cumsum(out[, 4]) / 1:nrow(out))) #cumulative mean
  colnames(out) <- c("time1", "time2", "x", "a", "status", "x_avg", "a_avg", "a_avg_lag")
  if(!is.na(id)) {out <- cbind(id, out)} #subject id, for bookkeeping
  
  out[1:min(nrow(out), max_t),]
}

#generate synthetic data
make_data <- function(n = 1000, ...) {
  data <- as.data.frame(foreach(i = 1:n, .combine = rbind) %do% 
                          gen_subject(id = i, ...))
  #add lag to dataframe
  data <- data %>%
    group_by(id) %>%
    mutate("a_lag1" = lag(a, default = .5),
           "x_avg_lag1" = lag(x_avg, default = 0)) %>%
    ungroup()
  
  data
}
