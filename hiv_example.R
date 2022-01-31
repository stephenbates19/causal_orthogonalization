library(ipw)
library(survival)
library(tidyverse)
library(rsample)

data(haartdat)
head(haartdat)

graphics.off()
ipwplot(weights = temp$ipw.weights, timevar = haartdat$fuptime,
        binwidth = 100, ylim = c(-1.5, 1.5), main = "Stabilized inverse probability weights")

summary(coxph(Surv(tstart, fuptime, event) ~ haartind + cluster(patient),
              data = haartdat, weights = temp$ipw.weights))

summary(coxph(Surv(tstart, fuptime, event) ~ haartind, data = haartdat))



#orthogonalization
hist(haartdat$cd4.sqrt)

#add lagged terms
dat2 <- haartdat %>%
  group_by(patient) %>%
  mutate(haartind_lag1 = lag(haartind, default = 0),
         haartind_lag2 = lag(haartind, n = 2, default = 0),
         haartind_lag3 = lag(haartind, n = 3, default = 0),
         cd4_lag1 = lag(cd4.sqrt, default = 22.9),
         cd4_lag2 = lag(cd4.sqrt, n = 2, default = 22.9),
         cd4.baseline = cd4.sqrt[1]) %>%
  ungroup()

head(dat2)

fit <- lm(cd4.sqrt ~ haartind_lag1 + haartind_lag2 + haartind_lag3 + cd4_lag1 + cd4_lag2 + cd4.baseline + sex + age, data = dat2)
dat2$cd4.sqrt_resid <- dat2$cd4.sqrt - predict(fit)

summary(coxph(Surv(tstart, fuptime, event) ~ haartind_lag1 + cd4.sqrt_resid + cd4_lag1 + cd4.baseline + sex + age, data = dat2))
ortho_beta_full <- summary(coxph(Surv(tstart, fuptime, event) ~ haartind_lag1 + cd4.sqrt_resid + cd4_lag1 + cd4.baseline + sex + age, data = dat2))$coeff[1]

summary(coxph(Surv(tstart, fuptime, event) ~ haartind_lag1 + cd4.sqrt + cd4.baseline + sex + age, data = dat2))


#ipw example
temp <- ipwtm(
  exposure = haartind_lag1,
  family = "survival",
  numerator = ~ sex + age,
  denominator = ~ sex + age + cd4_lag1,
  id = patient,
  tstart = tstart,
  timevar = fuptime,
  type = "first",
  trunc = .01,
  data = as.data.frame(dat2))

summary(coxph(Surv(tstart, fuptime, event) ~ haartind_lag1 + cluster(patient),
              data = dat2, weights = temp$ipw.weights))
ipw_beta_full <- summary(coxph(Surv(tstart, fuptime, event) ~ haartind + cluster(patient),
                               data = dat2, weights = temp$ipw.weights))$coeff[1]

summary(coxph(Surv(tstart, fuptime, event) ~ haartind_lag1 + cd4.baseline + cluster(patient),
              data = dat2, weights = temp$ipw.weights * temp2$ipw.weights))






#bootstrap
set.seed(1000)
n_boot <- 250

out <- c()
binf <- bootstraps(dat2, strata = patient, times = n_boot)
for(b in 1:n_boot) {
  if(b %% 25 == 0){print(b)}
  bdat <- dat2[binf$splits[[b]]$in_id, ]
  
  temp <- ipwtm(
    exposure = haartind_lag1,
    family = "survival",
    numerator = ~ sex + age,
    denominator = ~ sex + age + cd4_lag1,
    id = patient,
    tstart = tstart,
    timevar = fuptime,
    type = "first",
    data = as.data.frame(bdat))
  ipw_beta <- summary(coxph(Surv(tstart, fuptime, event) ~ haartind_lag1 + cluster(patient),
                            data = bdat, weights = temp$ipw.weights))$coeff[1]
  
  fit <- lm(cd4.sqrt ~ haartind_lag1 + haartind_lag2 + haartind_lag3 + cd4_lag1 + 
              cd4_lag2 + cd4.baseline + sex + age, data = bdat)
  bdat$cd4.sqrt_resid <- bdat$cd4.sqrt - predict(fit)
  ortho_beta <- summary(coxph(Surv(tstart, fuptime, event) ~ haartind_lag1 + cd4.sqrt_resid + 
                                cd4_lag1 + cd4.baseline + sex + age, data = bdat))$coeff[1]
  
  out <- rbind(out, c(ipw_beta, ortho_beta))
}
#save(out, file = "data/hiv_boot_4000.RData")
load(file = "data/hiv_boot_4000.RData")
out2 <- out[out[, 1] > -5, ] #clip outliers


c(ipw_beta_full, ortho_beta_full)
colMeans(out2)
apply(out2, 2, sd)

#z-scores of fitted coefficients
c(ipw_beta_full, ortho_beta_full) / apply(out2, 2, sd)

#basic bootstrap CIs
ipw_beta_full
ipw_beta_full - quantile(out2[, 1] - ipw_beta_full, probs = c(.05, .95))

ortho_beta_full
ortho_beta_full - quantile(out2[, 2] - ortho_beta_full, probs = c(.05, .95))


ipw_boot <- qplot(out2[, 1]) +
  geom_vline(xintercept = ipw_beta_full, color = "blue") +
  xlim(c(-5, 2)) +
  theme_bw(base_size = 11) + 
  labs(x = "IPW estimate") +
  theme(aspect.ratio = 1) + 
  theme(plot.title = element_text(size=11))
ipw_boot
ggsave(plot = ipw_boot, filename = "figures/hiv_ipw_boot.pdf", height = 2.5, width = 3)

ortho_boot <- qplot(out2[, 2]) +
  geom_vline(xintercept = ortho_beta_full, color = "blue") +
  theme_bw(base_size = 11) + 
  xlim(c(-5, 2)) +
  labs(x = "Ortho estimate") +
  theme(aspect.ratio = 1) + 
  theme(plot.title = element_text(size=11))
ortho_boot
ggsave(plot = ortho_boot, filename = "figures/hiv_ortho_boot.pdf", height = 2.5, width = 3)
