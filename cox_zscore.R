library(foreach)
library(survival)
library(tidyverse)
library(ggplot2)
library(ipw)
source("cox_simulator.R")

set.seed(11111)
n <- 5000
data <- make_data(n = n, tx = .5, beta = 1, beta2 = 1)
dim(data)
head(data, 10)

#naive regression
fit_naive <- coxph(Surv(time1, time2, status) ~ a + x + x_avg + a_avg, data = data)
summary(fit_naive)

data[["x_ortho"]] <- lm(x ~ a_lag1 + x_avg_lag1, data = data)$residuals
data <- data %>%
  group_by(id) %>%
  mutate("x_ortho_avg" = cummean(x_ortho))
fit_ortho <- coxph(Surv(time1, time2, status) ~ a + a_avg + x_ortho + x_ortho_avg, data = data)
summary(fit_ortho)

#ipw
ipw_fit <- ipwtm(exposure = a, 
                 family = "binomial", 
                 link = "logit",
                 numerator = ~ 1, 
                 denominator = ~ I(x > 0) + x + x_avg, 
                 id = id,  
                 tstart = time1,
                 timevar = time2,
                 type = "all",
                 data = as.data.frame(data),
                 trunc = .01)
#note: this only seems to work with very large n, otherwise this regression fails
fit_msm <- coxph(Surv(time1, time2, status) ~ a + a_avg + cluster(id), 
                 data = data,
                 weights = ipw_fit$ipw.weights)
summary(fit_msm)


#run simulations for teaser figure
set.seed(111)
n <- 2000
n_reps <- 50
out <- data.frame()
for(beta_xa in c(0, .25, .5)) {
  for(tx in c(0, .1, .2, .3)) {
    print(tx)
    for(rep_idx in 1:n_reps) {
      data <- make_data(n = n, max_t = 5, beta = 1, beta2 = 1, beta_xa = beta_xa, tx = tx)
      data <- data %>%
        group_by(id) %>%
        mutate("a_lag1" = lag(a, default = 0),
               "x_lag1" = lag(x, default = 0),
               "x_avg_lag1" = lag(x_avg, default = 0)) #add lag to dataframe
      
      #naive regression
      fit_naive <- coxph(Surv(time1, time2, status) ~ a + x + a_avg + x_avg , data = data, cluster = id)
      
      #ortho regression
      data[["x_ortho"]] <- lm(x ~ a_lag1 + x_avg_lag1, data = data)$residuals
      data <- data %>%
        group_by(id) %>%
        mutate("x_ortho_avg" = cummean(x_ortho))
      fit_ortho <- coxph(Surv(time1, time2, status) ~ a +  x_ortho + a_avg + x_ortho_avg, data = data, cluster = id)
      z_ortho <- summary(fit_ortho)$coefficients[,5]
      
      #ipw
      ipw_fit <- ipwtm(exposure = a, 
                       family = "binomial", 
                       link = "logit",
                       numerator = ~ 1, 
                       denominator = ~ I(x > 0) + x + x_avg, 
                       id = id,  
                       tstart = time1,
                       timevar = time2,
                       type = "all",
                       data = as.data.frame(data),
                       trunc = .001)
      fit_msm <- coxph(Surv(time1, time2, status) ~ a + a_avg + cluster(id), 
                       data = data,
                       weights = ipw_fit$ipw.weights)
      z_msm <- summary(fit_msm)$coefficients[,5]
      
      
      out <- rbind(out, c(tx, beta_xa, coef(fit_ortho)[1], coef(fit_msm)[1]))
    }
  }
}
colnames(out) <- c("tx_effect", "beta_xa", "RWR", "IPW")
head(out)

save(out, file = "data/cox_zscore.RData")
load(file = "data/cox_zscore.RData")

plt_dat <- pivot_longer(out, cols = c("RWR", "IPW")) %>%
  group_by(tx_effect, beta_xa, name) %>%
  mutate(zscore = value / sd(value))

levels_beta = c("None", "Moderate", "Strong")
names(levels_beta) = levels(factor(plt_dat$beta_xa))
plt <- ggplot(data = plt_dat, aes(x = as.factor(100*tx_effect))) + 
  geom_boxplot(aes(y = zscore, color = name)) + 
  geom_hline(yintercept = 0) +
  facet_wrap(~beta_xa, labeller = labeller(beta_xa = levels_beta)) +
  labs(x = "Treatment effect size (% reduced risk)", y = "z-score", color = "Method") +
  theme_bw() +
  theme(plot.title = element_text(size=11)) +
  theme(aspect.ratio = 1)
plt

ggsave(plot = plt, filename = "figures/cox_zscore.pdf", width = 6.5, height = 2.5)

