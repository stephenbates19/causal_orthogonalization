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
  
  #naive
  fit <- lm(U ~ A_1 + X_2 + A_2 - 1)
  betas <- c(betas, coef(fit)[1])
  
  #orthogonalization
  Xt_2 <- lm(X_2 ~ A_1 - 1)$residuals
  fit2 <- lm(U ~ A_1 + Xt_2 + A_2 - 1)
  betas2 <- c(betas2, coef(fit2)[1])
}

library(ggplot2)
library(latex2exp)
plt <- ggplot(data.frame("rhos" = rhos, "betas_naive" = betas, "betas_ortho" = betas2)) +
  geom_hline(yintercept = 0, color = "dark gray") +
  geom_line(aes(x = rhos, y = betas_naive, color = "naive")) + 
  geom_point(aes(x = rhos, y = betas_naive, color = "naive", shape = "naive")) +
  geom_line(aes(x = rhos, y = betas_ortho, color = "ortho")) + 
  geom_point(aes(x = rhos, y = betas_ortho, color = "ortho", shape = "ortho")) +
  labs(x = TeX("$$A_1 \\rightarrow$$ $$X_2$$ correlation"), y = TeX("$$A_1$$ coefficient"), 
       shape = "Method", color = "Method") +
  ggtitle("Linear Regression") +
  theme_bw(base_size = 10) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=11))
plt

ggsave(plt, file = "figures/linear_gaussian_simple.pdf", height = 2.25, width = 4)
