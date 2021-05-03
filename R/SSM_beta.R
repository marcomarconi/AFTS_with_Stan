## @knitr load_packages
require(tidyverse)
require(rstan)
library(xts)
library(quantmod)
library(FKF)
library(MARSS)
rstan_options(auto_write = TRUE)
rstan_options(javascript = FALSE)


## @knitr load_series
getSymbols("EWC", from="2006-04-26", to="2012-04-09")
getSymbols("EWA", from="2006-04-26", to="2012-04-09")
EWC_price <- EWC$EWC.Adjusted %>% as.vector
EWA_price <- EWA$EWA.Adjusted %>% as.vector
plot(EWC_price, EWA_price)

## @knitr fit_lm_series
fit_lm <- lm(EWC_price ~ EWA_price)
fit_lm
plot(fit_lm$residuals)

## @knitr fit_stan_series
SSM_beta <- stan_model("models/EpChan/SSM_beta.stan")
trainset <- 1:1250
testset <- 1251:1500
#SSM_beta_chan <- sampling(SSM_beta, data= list(N = length(trainset), y = EWC_price[trainset],  x = EWA_price[trainset]), iter = 500, chains=1)
SSM_beta_chan <- readRDS("models/EpChan/SSM_beta_chan.rds")
  
## @knitr print_stan_series
print(SSM_beta_chan, pars=c("sigma_y", "Sigma"),  digits=5)

## @knitr kalman_stan
# Retrive the estimates for the noise deviations
sigma_y <- extract(SSM_beta_chan, pars="sigma_y")[[1]] %>% mean
Sigma <- extract(SSM_beta_chan, pars="Sigma")[[1]] %>% apply(., c(2, 3), mean)
# the initial values for the filter
a0=c(1,1);P0=matrix(c(0,0,0,0), byrow = T, ncol=2)
# Now fill all the matrices according to the package documentation https://cran.r-project.org/web/packages/FKF/FKF.pdf
N <- length(EWC_price)
y <- EWC_price
Zt <- array(data = t(cbind(1, EWA_price)), dim = c(1, 2, N))
Tt <- array(data = c(1,0,0,1), dim = c(2, 2, N))
ct <- array(0, dim = c(1, N))
dt <- array(0, dim = c(2, N))
HHt = array(data = Sigma, dim = c(2, 2, N))
GGt = array(data = matrix(sigma_y^2), dim = c(1, 1, N))
fkf.obj <- fkf(a0=a0, P0=P0, ct = ct, dt = dt, Tt = Tt, Zt = Zt, HHt = HHt, GGt = GGt, yt = rbind(y))
par(mfrow=c(1,2))
plot(index(EWC), fkf.obj$att[1,], col="blue", type="l", ylab="Offset", xlab="Date")
plot(index(EWC), fkf.obj$att[2,], col="blue", type="l", ylab="Hedge Ratio", xlab="Date")

## @knitr  fit_stan_kf_series
SSM_beta_KF <- stan_model("models/EpChan/SSM_beta_KF.stan")
trainset <- 1:1250
testset <- 1251:1500
#SSM_beta_KF_chan <- sampling(SSM_beta_KF, data= list(N = length(trainset), y = EWC_price[trainset],  x = EWA_price[trainset], m0=c(1,1), P0=matrix(0, nrow=2, ncol=2)), iter = 500, chains=1)
SSM_beta_KF_chan <- readRDS("models/EpChan/SSM_beta_KF_chan.rds")
print(SSM_beta_KF_chan, pars=c("R", "Q"),  digits=5)

## @knitr  fit_MARSS_series
trainset <- 1:1250
TT <- length(trainset)
m <- 2
y <- EWC_price[trainset] 
x <- EWA_price[trainset] 
B <- diag(m)  
U <- matrix(0, nrow = m, ncol = 1) 
Q <- matrix("c", m, m)  
diag(Q) <- c("q.alpha", "q.beta")  
Z <- array(NA, c(1, m, TT)) 
Z[1, 1, ] <- rep(1, TT)  
Z[1, 2, ] <- x 
A <- matrix(0)  
R <- matrix("r") 
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)
dlm_chan <- MARSS(y, inits = inits_list, model = mod_list)

## @knitr kalman_MARSS
# Retrive the estimates for the noise deviations
sigma_y <- sqrt(dlm_chan$par$R)
Sigma <- matrix(c(dlm_chan$par$Q[1], dlm_chan$par$Q[2], dlm_chan$par$Q[2], dlm_chan$par$Q[3]), byrow=TRUE, nrow=2, ncol=2)
# the initial values for the filter
a0=c(1,1);P0=matrix(c(0,0,0,0), byrow = T, ncol=2)
# Now fill all the matrices according to the package documentation https://cran.r-project.org/web/packages/FKF/FKF.pdf
N <- length(EWC_price)
y <- EWC_price
Zt <- array(data = t(cbind(1, EWA_price)), dim = c(1, 2, N))
Tt <- array(data = c(1,0,0,1), dim = c(2, 2, N))
ct <- array(0, dim = c(1, N))
dt <- array(0, dim = c(2, N))
HHt = array(data = Sigma, dim = c(2, 2, N))
GGt = array(data = matrix(sigma_y^2), dim = c(1, 1, N))
fkf.obj <- fkf(a0=a0, P0=P0, ct = ct, dt = dt, Tt = Tt, Zt = Zt, HHt = HHt, GGt = GGt, yt = rbind(y))
par(mfrow=c(1,2))
plot(index(EWC), fkf.obj$att[1,], col="blue", type="l", ylab="Offset", xlab="Date")
plot(index(EWC), fkf.obj$att[2,], col="blue", type="l", ylab="Hedge Ratio", xlab="Date")

## @knitr plot_trades
plot(index(EWC)[-(1:3)], fkf.obj$vt[1,-(1:3)], ylab="Forecast error", xlab="Date" )
lines(index(EWC)[-(1:3)], sqrt(fkf.obj$Ft[1,1,-(1:3)]), col="blue")
lines(index(EWC)[-(1:3)], -sqrt(fkf.obj$Ft[1,1,-(1:3)]), col="blue")

## @knitr plot_cumret
v <- fkf.obj$vt
S <- fkf.obj$Ft[1,1,]
y2 <- cbind(EWA_price, EWC_price)
shortsExit <- v < sqrt(S);
shortsEntry <- v > sqrt(S);
longsExit <- v > -sqrt(S);
longsEntry <- v < -sqrt(S);
numUnitsLong=vector(length = nrow(EWC));
numUnitsShort=vector(length = nrow(EWC));
numUnitsLong[1]=0;
numUnitsLong[longsEntry]=1;
numUnitsLong[longsExit]=0;
numUnitsShort[1]=0;
numUnitsShort[shortsEntry]=-1;
numUnitsShort[shortsExit]=0;
numUnits=numUnitsLong+numUnitsShort;
positions <- numUnits * cbind(-fkf.obj$att[1,], 1) * y2
pnl <- rowSums(lag(positions, 1) * (y2-lag(y2, 1)) / lag(y2, 1), na.rm = TRUE)
ret <- pnl / rowSums(abs(lag(positions, 1)), na.rm = TRUE);
ret[is.infinite(ret)] <- 0
ret[is.nan(ret)] <- 0
par(mfrow=c(1,2))
plot(index(EWC)[trainset], (cumprod(1+ret[trainset])-1),ylab="Forecast error", xlab="Date", type="l")
plot(index(EWC)[testset], (cumprod(1+ret[testset])-1), ylab="Forecast error", xlab="Date", type="l")



