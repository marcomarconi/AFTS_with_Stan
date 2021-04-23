library(forecast)
library(tidyverse)
library(rstan)
library(loo)
library(xts)
library(MASS)

setwd("chapter 3/")

# fake data ARCH
set.seed(10)
N <- 1000
mu <- 0.02
alpha <- c(0.01, 1)
err <- rnorm(N, 0, 1) 
sigma_squared <- vector()
shocks <- vector()
shocks[1] <- 0
returns <- vector()
returns[1] <- 0
for(i in 2:N) {
  sigma_squared[i] <- alpha[1] + alpha[2] * (shocks[i-1]^2)
  shocks[i] <- sqrt(sigma_squared[i]) * err[i]
  returns[i] <- mu + shocks[i]
}


fit_arch <- stan(file = "ARCH.stan", data = list(T = N, r = returns, order = 1), chains = 4, cores = 4)
# test the errors are white noise
sq <- extract(fit_arch, pars="sigma_squared")[[1]] %>% colMeans()
mu <- extract(fit_arch, pars="mu")[[1]] %>% mean()
((returns - mu) / sqrt(sq)) %>% Box.test()
# model checking
preds <- extract(fit_arch, pars="r_hat")[[1]] 
preds <- cbind(
  Estimate = colMeans(preds),
  Q5 = apply(preds, 2, quantile, probs = 0.05),
  Q95 = apply(preds, 2, quantile, probs = 0.95)
)
ggplot(data.frame(year=1:N, y=returns, preds), aes(x = year, y = Estimate)) +
  geom_smooth(aes(ymin = Q5, ymax = Q95), stat = "identity", size = 0.5) +
  geom_line(aes(y = y)) 


# 3.1
intel <- read.table("m-intc7308.txt", header=T)
intel <- xts(intel$rtn, as.Date(as.character(intel$date), format = "%Y%m%d")) 
fit_intel_1 <- stan(file = "ARCH.stan", data = list(T = length(intel), r = intel %>% as.vector(), order = 1, family = 0), chains = 4, cores = 4)
fit_intel_3 <- stan(file = "ARCH.stan", data = list(T = length(intel), r = intel %>% as.vector(), order = 3, family = 0), chains = 4, cores = 4)
fit_intel_t <- stan(file = "ARCH.stan", data = list(T = length(intel), r = intel %>% as.vector(), order = 1, family = 1), chains = 4, cores = 4)
log_lik_fit_intel_1 <- extract_log_lik(fit_intel_1, merge_chains = FALSE)
log_lik_fit_intel_3 <- extract_log_lik(fit_intel_3, merge_chains = FALSE)
log_lik_fit_intel_t <- extract_log_lik(fit_intel_t, merge_chains = FALSE)
r_eff_1 <- relative_eff(exp(log_lik_fit_intel_1), cores = 2)
r_eff_3 <- relative_eff(exp(log_lik_fit_intel_3), cores = 2)
r_eff_t <- relative_eff(exp(log_lik_fit_intel_t), cores = 2)
loo_compare(loo(log_lik_fit_intel_1, r_eff = r_eff_1), loo(log_lik_fit_intel_3, r_eff = r_eff_3), loo(log_lik_fit_intel_t, r_eff = r_eff_t))

# 3.2
MRKUSD <- scan("exch-perc.txt")
fit_MRKUSD <- stan(file = "ARCH.stan", data = list(T = length(MRKUSD), r = MRKUSD %>% as.vector(), order = 3, family = 0), chains = 4, cores = 4)


#fake data CHARMA
N <- 200
mu <- 0.02
m = 2
Cov <- matrix(nrow = 2, ncol = 2)
Cov[1,1] <- 0.3
Cov[2,2] <- 0.6
Cov[2,1] <- Cov[1,2] <- 0.15
delta <- mvrnorm(1, rep(0, m), Cov)
sigma_eta <- 0.1
sigma_squared <- vector()
shocks <- vector()
shocks[1:m] <- 0
returns <- vector()
returns[1:m] <- 0
for(t in (m+1):N) {
  sigma_squared[t] <- sigma_eta^2 + shocks[(t-m):(t-1)]  %*%  Cov  %*%  shocks[(t-m):(t-1)]
  shocks[t] <- shocks[(t-m):(t-1)]  %*%  delta + rnorm(1, 0, sqrt(sigma_squared[t]))
  returns[t] <- mu + shocks[t]
}

# fake data SVM
N <- 500
mu_h <- -1.02
phi <- 0.95
sigma <- 0.25
err <- rnorm(N, 0, 1) 
delta <- rnorm(N, 0, 1) 
h <- vector()
h[1] <- rnorm(1, mu_h, sigma / sqrt(1 - phi^2))
shocks <- vector()
shocks[1] <- 0
returns <- vector()
returns[1] <- 0
for(i in 2:N) {
  h[i] <- rnorm(1, mu_h + phi*(h[i-1] - mu_h), sigma)
  returns[i] <- rnorm(1, 0, exp(h[i]/2))
}


