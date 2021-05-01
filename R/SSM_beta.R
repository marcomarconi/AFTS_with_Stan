source('common.R')


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
SSS_beta <- stan_model("models/EpChan/SSM_beta.stan")
trainset <- 1:1250
testset <- 1251:1500
#SSS_beta_chan <- sampling(SSS_beta, data= list(N = length(trainset), y = EWC_price[trainset],  x = EWA_price[trainset]), iter = 500, chains=1)
SSS_beta_chan <- readRDS("models/EpChan/SSS_beta_chan.rds")
  
## @knitr print_stan_series
print(SSS_beta_chan, pars=c("sigma_y", "Sigma"),  digits=5)

## @knitr kalman
# Retrive the estimates for the noise deviations
sigma_y <- rstan::extract(SSS_beta_chan, pars="sigma_y")[[1]] %>% mean
Sigma <- rstan::extract(SSS_beta_chan, pars="Sigma")[[1]] %>% apply(., c(2, 3), mean)
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
shortsExit <- v < (S/2);
shortsEntry <- v > (S/2);
longsExit <- v > -(S/2);
longsEntry <- v < -(S/2);
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



