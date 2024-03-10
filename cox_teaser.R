n <- 40000
rho_xa <- .5
rho_ux <- .5

set.seed(1)
rhos <- c(0,.2,.4,.6,.8)
betas <- c()
betas2 <- c()
for(rho_ax in rhos) {
  
  #simulate data
  A_1 <- rnorm(n)
  U <- rnorm(n)
  X_2 <- rho_ax * A_1 + rho_ux * U + sqrt(1 - rho_ax^2 - rho_ux^2) * rnorm(n)
  A_2 <- rho_xa * X_2 + sqrt(1 - rho_xa^2) * rnorm(n)
  
  # survival outcome
  avg_surv <- 1.2 - pnorm(U)
  Y <- sapply(1:n, function(i){rexp(n = 1, rate = 1/avg_surv[i])})
  
  #orthogonalization
  Xt_2 <- lm(X_2 ~ A_1 - 1)$residuals
  
  max_t <- 5
  long_dat <- c()
  for(i in 1:n){
    t <- 0
    while(t < max_t) {
      long_dat <- rbind(long_dat, c(A_1[i], X_2[i], A_2[i], Xt_2[i], Y[i], t, t+1, Y[i] < t+1))
      t <- t + 1
      if(Y[i] < t + 1) {
        t <- max_t
      } #event happened
    }
  }
  colnames(long_dat) <- c("A_1", "X_2", "A_2", "Xt_2", "Y", "time", "time2", "event")
  
  #naive fit
  fit <- coxph(Surv(time, time2, event) ~ A_1 + X_2 + A_2, 
               data = as.data.frame(long_dat))
  betas <- c(betas, coef(fit)[1])
  
  #ortho fit
  fit_ortho <- coxph(Surv(time, time2, event) ~ A_1 + Xt_2 + A_2, 
               data = as.data.frame(long_dat))
  betas2 <- c(betas2, coef(fit_ortho)[1])
}

library(ggplot2)
library(latex2exp)
plt2 <- ggplot(data.frame("rhos" = rhos, "betas_naive" = betas, "betas_ortho" = betas2)) +
  geom_hline(yintercept = 0, color = "dark gray") +
  geom_line(aes(x = rhos, y = betas_naive, color = "naive")) + 
  geom_point(aes(x = rhos, y = betas_naive, color = "naive", shape = "naive")) +
  geom_line(aes(x = rhos, y = betas_ortho, color = "RWR")) + 
  geom_point(aes(x = rhos, y = betas_ortho, color = "RWR", shape = "RWR")) +
  labs(x = TeX("$$A_1 \\rightarrow$$ $$X_2$$ correlation"), y = TeX("$$A_1$$ coefficient"), 
       shape = "Method", color = "Method", title = "Cox Regression") +
  theme_bw(base_size = 10) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=11))
plt2

#run after linear_gaussian.R to combine plots
library(ggpubr)
double_plot <- ggarrange(plt, plt2, ncol = 2, common.legend = T, legend = "bottom")
double_plot
ggsave(double_plot, file = "figures/double_teaser.pdf", height = 2.75, width = 6.0)

