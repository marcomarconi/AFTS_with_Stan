

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




# fake data GARCH
set.seed(10)
N <- 1000
mu <- 0.02
alpha <- c(0.01, 0.4)
beta <- c(0.03)
err <- rnorm(N, 0, 1) 
sigma <- vector()
sigma[1:2] <- sqrt(alpha[1])
shocks <- vector()
shocks[1:2] <- 0
returns <- vector()
returns[1:2] <- 0
for(i in 2:N) {
  sigma[i] <- sqrt(alpha[1] + alpha[2] * (shocks[i-1]^2) + beta[1] * (sigma[i-1]^2))
  shocks[i] <- sigma[i] * err[i]
  returns[i] <- mu + shocks[i]
}

# 3.3
sp500 <- scan("sp500.dat.txt")
fit_sp500_garch <- stan(file = "GARCH.stan", data = list(N = length(sp500), y = sp500), chains = 4, cores = 4)
fit_sp500_ar_garch <- stan(file = "AR-GARCH.stan", data = list(N = length(sp500), y = sp500, Kar=3), chains = 4, cores = 4)
loo_compare(loo(fit_sp500_garch), loo(fit_sp500_ar_garch)) # AR-GARCH is actually better with LOO-



# fake data IGARCH
N <- 1000
mu <- 0.02
alpha <- c(0.01)
beta <- c(0.03)
err <- rnorm(N, 0, 1) 
sigma <- vector()
sigma[1] <- sqrt(alpha[1])
shocks <- vector()
shocks[1] <- 0
returns <- vector()
returns[1] <- 0
for(i in 2:N) {
  sigma[i] <- sqrt(alpha + beta * shocks[i-1]^2 + (1-beta) * (sigma[i-1]^2))
  shocks[i] <- sigma[i] * err[i]
  returns[i] <- mu + shocks[i]
}

# fake data GARCH-M
N <- 1000
mu <- 0.02
premium <- 0.5
alpha <- c(0.01, 0.4)
beta <- c(0.03)
err <- rnorm(N, 0, 1) 
sigma <- vector()
sigma[1] <- sqrt(alpha[1])
shocks <- vector()
shocks[1] <- 0
returns <- vector()
returns[1] <- 0
for(i in 2:N) {
  sigma[i] <- sqrt(alpha[1] + alpha[2] * (shocks[i-1]^2) + beta[1] * (sigma[i-1]^2))
  shocks[i] <- sigma[i] * err[i]
  returns[i] <- mu + premium * (sigma[i])^2 + shocks[i]
}

# fake data EGARCH
N <- 200
mu <- 0.0118
alpha <- c(-0.557, 0.220)
beta <- c(0.929)
gamma <- -0.264
err <- rnorm(N, 0, 1) 
sigma <- vector()
sigma[1] <- 0.06
shocks <- vector()
shocks[1] <- 0
returns <- vector()
returns[1] <- 0
for(i in 2:N) {
  sigma[i] <- sqrt(exp(alpha[1] + 
                         alpha[2] * (abs(shocks[i-1]) + gamma*shocks[i-1]) / sigma[i-1] +
                         beta[1] * log(sigma[i-1]^2)))
  shocks[i] <- sigma[i] * err[i]
  returns[i] <- mu + shocks[i]
}


# fake data TGARCH
N <- 200
mu <- 0.0118
alpha <- c(0.557, 0.220)
beta <- c(0.0929)
gamma <- 0.264
err <- rnorm(N, 0, 1) 
sigma <- vector()
sigma[1] <- 1
shocks <- vector()
shocks[1] <- 0
returns <- vector()
returns[1] <- 0
for(i in 2:N) {
  sigma[i] <- sqrt(alpha[1] + 
                     (alpha[2] + gamma*(shocks[i-1]<0))*shocks[i-1]^2  +
                     beta[1] * sigma[i-1]^2
  )
  shocks[i] <- sigma[i] * err[i]
  returns[i] <- mu + shocks[i]
}

# fake data RGARCH (requires an additional sigma proxy vector)
set.seed(10)
N <- 500
mu <- 0.02
alpha0 <- 0.01
alpha1 <- 0.4
alpha2 <- 0.66
beta1 <- 0.03
err <- rnorm(N, 0, 1) 
sigma <- vector()
sigma[1:2] <- sqrt(alpha[1])
shocks <- vector()
shocks[1:2] <- 0
#proxy <- ...
returns <- vector()
returns[1:2] <- 0
for(i in 2:N) {
  sigma[i] <- sqrt(alpha0 +  alpha1 * (shocks[i-1]^2) + alpha2 * (proxy[i-1]^2) + beta1 * (sigma[i-1]^2))
  shocks[i] <- sigma[i] * err[i]
  returns[i] <- mu + shocks[i]
}


# fake data SVM (STAN manual model)
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

# 3.4
ibm_sp500 <- read.table("chapter 3/m-ibmspln.dat", header=T)
fit_argarch_ibm2 <- stan(file = "chapter 3/AR-GARCH.stan" , data=c(list(T = length(ibm_sp500$ibm), r = ibm_sp500$ibm, ar_order=1,M = 0, m = vector())), chains=4, cores=4)
fit_argarch_I_ibm2 <- stan(file = "chapter 3/AR-GARCH-I.stan" , data=c(list(T = length(ibm_sp500$ibm), r = ibm_sp500$ibm, u=ibm_sp500_u$summer, ar_order = 1, M = 0, m = vector())), chains=4, cores=4, seed = SEED)
loo_compare(loo(fit_argarch_ibm2), loo(fit_argarch_I_ibm2))

# 3.5
fit_garch_x_sp500 <- stan(file = "chapter 3/GARCH_x.stan" , data=c(list(N = length(ibm_sp500$sp), y = ibm_sp500$sp, x=ibm_sp500$ibm, Kar = 2, M = 0, m = vector())), chains=4, cores=4)


# fake data ARMA-GARCH
set.seed(10)
N <- 2000
mu <- 0.02
phi <- 0.22
theta <- 0.06
alpha <- c(0.01, 0.4)
beta <- c(0.03)
premium <- 1
sigma <- vector()
sigma[1] <- sqrt(alpha[1] / (1 - (alpha[2] + beta[1])))
returns <- vector()
returns[1] <- 0
sigma <- 0.1
returns <- vector()
returns[1] <- mu
for(i in 2:N) {
  sigma[i] <- sqrt(alpha[1] + alpha[2] * (returns[i-1] - mu)^2 + beta[1] * sigma[i-1]^2)
  returns[i] <- rnorm(1, mu + premium*sigma[i]^2 + phi * returns[i-1] + theta * (returns[i-1] - mu), sigma[i]);
}


# fake data ARMA-GARCH test
set.seed(10)
N <- 200
mu <- 0.02
phi <- 0.22
theta <- 0.06
alpha <- c(0.01, 0.4)
beta <- c(0.03)
premium <- 0
sigma <- vector()
sigma[1] <- sqrt(alpha[1] / (1 - (alpha[2] + beta[1])))
returns <- vector()
returns[1] <- 0
sigma <- 0.1
returns <- vector()
returns[1] <- mu
for(i in 2:N) {
  sigma[i] <- sqrt(alpha[1] + alpha[2] * (returns[i-1] - mu)^2 + beta[1] * sigma[i-1]^2)
  returns[i] <- rnorm(1, mu + premium*sigma[i]^2 + phi * returns[i-1] + theta * (returns[i-1] - mu), sigma[i]);
}
