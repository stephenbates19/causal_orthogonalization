library(tidyverse)
library(ggplot2)
library(latex2exp)

#make data from At -> Xt chain graph with phantom
make_data <- function(tsteps = 2, n = 10000, 
                      ax_effect = .5, ux_effect = .5, xa_effect = .8,
                      unoise = rnorm, xnoise = rnorm, ynoise = rnorm) {
  
  #time step 1
  U <- unoise(n)
  X <- ux_effect * U + sqrt(1 - ux_effect^2) * xnoise(n)
  A <- 2*((X > 0) - .5)
  flip_idx <- rbinom(n, 1, 1 - xa_effect)
  A[flip_idx == 1] <- (2*rbinom(sum(flip_idx), 1, .5) - 1)
  Y <- U + ynoise(n) #assymetric, heavy-tailed
  out <- data.frame("Y" = Y, "U" = U, "A_1" = A, "X_1" = X)
  
  #remaining time steps
  for(t in 2:tsteps) {
    X <- ax_effect* A + ux_effect *U + sqrt(1 - ax_effect^2 - ux_effect^2) * xnoise(n)
    A <- 2*((X > 0) - .5)
    flip_idx <- rbinom(n, 1, 1 - xa_effect)
    A[flip_idx == 1] <- (2*rbinom(sum(flip_idx), 1, .5) - 1)
    out[[paste0("A_",t)]] <- A
    out[[paste0("X_",t)]] <- X
  }
  
  return(out)
}

#fit one ortho model as a test
test <- make_data(n = 1000000, tsteps = 5, ynoise = function(n){rt(n, df = 5) + rexp(n) - 1}, xa_effect = 0)
dim(test)
head(test)

cor(test$X_1, test$A_2)

#naive regression
lm(Y ~ A_1 + X_1 + A_2 + X_2, data = test)

#ortho regression
test[["Xt_1"]] <- lm(X_1 ~ 1, data = test)$residuals
test[["Xt_2"]] <- lm(X_2 ~ A_1 + X_1, data = test)$residuals
lm(Y ~ A_1 + Xt_1 + A_2 + Xt_2, data = test)



#carry out the simulation over many settings
steps <- 5
results <- data.frame()
set.seed(10)
for(esize in c(0, .1, .2, .3 , .4, .5, .6, .7, .8)) {
  print(esize)
  dat <- make_data(n = 1000000, ax_effect = esize,
                   tsteps = steps,
                   ynoise = function(n){rt(n, df = 5) + rexp(n) - 1})
  
  #naive regression
  naive_coeffs <- coefficients(lm(Y ~ . - U, data = dat))
  worst_naive <- max(abs(naive_coeffs[1:5 * 2]))
  
  #ortho regression
  dat[["Xt_1"]] <- lm(X_1 ~ 1, data = dat)$residuals
  dat[["Xt_2"]] <- lm(X_2 ~ A_1 + X_1, data = dat)$residuals
  dat[["Xt_3"]] <- lm(X_3 ~ A_1 + X_1 + A_2 + X_2, data = dat)$residuals
  dat[["Xt_4"]] <- lm(X_4 ~ A_1 + X_1 + A_2 + X_2 + A_3 + X_3, data = dat)$residuals
  dat[["Xt_5"]] <- lm(X_5 ~ A_1 + X_1 + A_2 + X_2 + A_3 + X_3 + A_4 + X_4, data = dat)$residuals
  ortho_coeffs <- coefficients(lm(Y ~ A_1 + Xt_1 + A_2 + Xt_2 + A_3 + 
                                    Xt_3 + A_4 + Xt_4 + A_5 + Xt_5, data = dat))
  worst_ortho <- max(abs(ortho_coeffs[1:5 * 2]))
  
  results <- rbind(results, c("ax_effect" = esize, 
                              "naive_coef" = worst_naive, 
                              "ortho_coef" = worst_ortho))
}
colnames(results) <- c("ax_effect", "naive_coef", "ortho_coef")
results


plt <- ggplot(results, aes(x = ax_effect)) + 
  geom_point(aes(y = naive_coef, color = "naive", shape = "naive")) +
  geom_line(aes(y = naive_coef, color = "naive")) + 
  geom_point(aes(y = ortho_coef, color = "RWR", shape = "RWR")) +
  geom_line(aes(y = ortho_coef, color = "RWR")) + 
  geom_hline(yintercept = 0) + 
  theme_bw(base_size = 10) +
  labs(x = TeX("$$A_t \\rightarrow$$ $$X_{t+1}$$ correlation"), 
       y = TeX("Max $$A_t$$ coefficient"), color = "Method", shape = "Method") + 
  theme(aspect.ratio = 1) 
plt

ggsave(plt, file = "figures/linear_nongaussian.pdf", width = 3.5, height = 2.5)


### Compare variance to IPW
n <- 1000
steps <- 5
results <- data.frame()
set.seed(10)
xa_effect <- 0
ax_effect <- .1
ux_effect <- .9
for(xa_effect in c(0, .9)) {
  for(ux_effect in c(0, .25, .5, .7, .9, .95)) {
    print(ux_effect)
    for(i in 1:500) {
      dat <- make_data(n = n, ax_effect = ax_effect , ux_effect = ux_effect, xa_effect = xa_effect,
                       tsteps = steps,
                       ynoise = function(n){0})
      
      #naive regression
      naive_coeffs <- coefficients(lm(Y ~ . - U, data = dat))
      worst_naive <- max(abs(naive_coeffs[1:5 * 2]))
      
      #ortho regression
      dat[["Xt_1"]] <- lm(X_1 ~ 1, data = dat)$residuals
      dat[["Xt_2"]] <- lm(X_2 ~ A_1 + X_1, data = dat)$residuals
      dat[["Xt_3"]] <- lm(X_3 ~ A_1 + X_1 + A_2 + X_2, data = dat)$residuals
      dat[["Xt_4"]] <- lm(X_4 ~ A_1 + X_1 + A_2 + X_2 + A_3 + X_3 , data = dat)$residuals
      dat[["Xt_5"]] <- lm(X_5 ~ A_1 + X_1 + A_2 + X_2 + A_3 + X_3 + A_4 + X_4, data = dat)$residuals
      ortho_coeffs <- coefficients(lm(Y ~ A_1 + Xt_1 + A_2 + Xt_2 + A_3 + 
                                        Xt_3 + A_4 + Xt_4 + A_5 + Xt_5, data = dat))
      worst_ortho <- max(abs(ortho_coeffs[1:5 * 2]))
      
      #msm (marginal) fit
      n_match <- ((dat$X_1 > 0) == dat$A_1) + 
        ((dat$X_2 > 0) == dat$A_2) + 
        ((dat$X_3 > 0) == dat$A_3) + 
        ((dat$X_4 > 0) == dat$A_4) + 
        ((dat$X_5 > 0) == dat$A_5)
      w <- (.5 + xa_effect / 2)^n_match * (.5 - xa_effect / 2)^(steps - n_match) 
      ipw_coeffs <- coefficients(lm(Y ~ A_1 + A_2 +  A_3 +  A_4 + A_5, data = dat, weights = 1/w))
      
      results <- rbind(results, c(ax_effect, ux_effect, xa_effect,
                                  "naive_coef" = naive_coeffs[6], 
                                  "ortho_coef" = ortho_coeffs[6],
                                  "ipw_coef" = ipw_coeffs[4]
      ))
    }
  }
}
colnames(results) <- c("ax", "ux", "xa", "naive_coef", "ortho_coef", "ipw_coef")
#results

results2 <- results %>% 
  group_by(ax, ux, xa) %>%
  summarize(var_ortho = var(ortho_coef),
            var_ipw = var(ipw_coef),
            mean_ortho = mean(ortho_coef),
            mean_ipw = mean(ipw_coef)) %>%
  mutate(rel_var = var_ipw / var_ortho)
results2


plt <- ggplot(data = results2,
              aes(x = ux, y = rel_var, color = as.factor(xa), shape = as.factor(xa))) +
  geom_point() + 
  geom_line() +
  geom_hline(yintercept = 1, color = "black") +
  scale_y_log10() +
  scale_color_discrete(name = "Overlap", labels = c("high", "low")) +
  scale_shape_discrete(name = "Overlap", labels = c("high", "low")) +
  theme_bw(base_size = 10) +
  labs(x = TeX("$$U \\rightarrow$$ $$X_{t}$$ correlation"), y = "IPW relative variance") +
  theme(aspect.ratio = 1)
plt
ggsave(plt, file = "figures/linear_nongaussian_ipw.pdf", width = 3.5, height = 2.5)
