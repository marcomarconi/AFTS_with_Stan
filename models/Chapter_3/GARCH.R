library(forecast)
library(tidyverse)
library(rstan)
library(loo)
library(xts)

# fake data GARCH
set.seed(10)
N <- 1000
mu <- 0.02
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
  returns[i] <- mu + shocks[i]
}

# 3.3
sp500 <- scan("sp500.dat.txt")
fit_sp500_garch <- stan(file = "GARCH.stan", data = list(T = length(sp500), r = sp500), chains = 4, cores = 4)
fit_sp500_ar_garch <- stan(file = "AR-GARCH.stan", data = list(T = length(sp500), r = sp500, ar_order=3), chains = 4, cores = 4)
loo_compare(loo(fit_sp500_garch), loo(fit_sp500_ar_garch)) # AR-GARCH is actually better



# fake data IGARCH
N <- 1000
mu <- 0.02
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

