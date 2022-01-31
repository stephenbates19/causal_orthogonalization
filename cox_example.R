library(foreach)
library(survival)
library(tidyverse)
library(ggplot)
library(ipw)

#data simulator
source("cox_simulator.R")
head(gen_subject(id = 1))

#generate synthetic data
set.seed(1111)
n <- 5000
data <- as.data.frame(foreach(i = 1:n, .combine = rbind) %do% 
                        gen_subject(id = i, max_t = 5, beta = 2, beta2 = 2))
data <- data %>%
  group_by(id) %>%
  mutate("a_lag1" = lag(a, default = 0),
         "x_avg_lag1" = lag(x_avg, default = 0)) #add lag to dataframe

dim(data)
head(data, 10)

#naive regression
fit <- coxph(Surv(time1, time2, status) ~ x + a + a_lag1, data = data)
summary(fit)

#ortho regression
data[["x_ortho"]] <- lm(x ~ a_lag1, data = data)$residuals
data <- data %>%
  group_by(id) %>%
  mutate("x_ortho_avg" = cummean(x_ortho))
fit <- coxph(Surv(time1, time2, status) ~ x_ortho + a + a_lag1, data = data)
summary(fit)

reps <- 250
set.seed(333)
n <- 2000
out <- data.frame()
t1 <- Sys.time()
for(beta2 in c(0, .5, 1.1, 1.8, 3.1)) {
  print(beta2)
  for(j in 1:reps) {
    data <- as.data.frame(foreach(i = 1:n, .combine = rbind) %do% 
                            gen_subject(id = i, max_t = 5, beta = 1, beta2 = beta2))
    data <- data %>%
      group_by(id) %>%
      mutate("a_lag1" = lag(a, default = 0),
             "x_avg_lag1" = lag(x_avg, default = 0)) #add lag to dataframe
    
    #naive regression
    fit_naive <- coxph(Surv(time1, time2, status) ~ a_avg + x_avg + x + a, data = data)
    
    #ortho regression
    data[["x_ortho"]] <- lm(x ~ a_lag1 + x_avg_lag1, data = data)$residuals
    data <- data %>%
      group_by(id) %>%
      mutate("x_ortho_avg" = cummean(x_ortho))
    fit_ortho <- coxph(Surv(time1, time2, status) ~ a_avg + x_ortho_avg + x_ortho + a, data = data)
    
    out <- rbind(out, c(beta2, cor(data$x, data$a_lag1), coef(fit_naive)[1:2], coef(fit_ortho)[1:2]))
  }
}
print(Sys.time() - t1)
colnames(out) <- c("ax_effect", "ax_cor", "naive_aavg", "naive_xavg", "ortho_aavg", "ortho_xavg")
save(out, file = "data/cox_example.RData")

out2 <- out %>%
  group_by(ax_effect) %>%
  summarize(ax_cor = mean(ax_cor),
            naive_aavg = mean(naive_aavg), 
            naive_xavg = mean(naive_xavg),
            ortho_aavg = mean(ortho_aavg),
            ortho_xavg = mean(ortho_xavg))

plt <- ggplot(data = out2, aes(x = ax_cor)) + 
  geom_point(aes(y = naive_aavg, color = "naive", shape = "naive")) + 
  geom_line(aes(y = naive_aavg, color = "naive")) + 
  geom_point(aes(y = ortho_aavg, color = "ortho", shape = "ortho")) + 
  geom_line(aes(y = ortho_aavg, color = "ortho")) + 
  geom_hline(yintercept = 0) +
  labs(x = TeX("$$A_t \\rightarrow$$ $$X_{t+1}$$ correlation"), 
       y = "Regression coefficient", shape = "Method", color = "Method") +
  theme_bw() +
  theme(plot.title = element_text(size=11)) +
  theme(aspect.ratio = 1)
plt

ggsave(plt, file = "figures/cox_consistency.pdf", width = 3.5, height = 2.25)

