source('common.R')

## @knitr load_8.4
ibmsp2608 <- read.table("data/m-ibmsp2608.txt", header = T)
ibmsp2608$date <- ibmsp2608$date %>% as.character() %>%  as.Date(format="%Y%m%d")
ibmsp2608[,2:3] <- ibmsp2608[,2:3] * 100
par(mfrow=c(2,1), mar=c(2,5,1,4))
plot(ibmsp2608$date, ibmsp2608[,2], type="l", ylab="IBM", xlab="Date", cex.lab=2, cex.axis=1.5)
plot(ibmsp2608$date, ibmsp2608[,3], type="l", ylab="S&P 500", xlab="Date", cex.lab=2, cex.axis=1.5)

## @knitr fit_8.4
VAR_QR <- stan_model("models/Chapter_8/VAR_QR.stan")
fit_VAR_QR <- readRDS("models/Chapter_8/saved_RDS/fit_VAR_QR.rds")

## @knitr print_8.4
print(fit_VAR_QR, pars=c("mu", "Sigma", "phi"))

## @knitr corr_8.4
Sigma <- extract(fit_VAR_QR, pars=c("Sigma"))[[1]] 
corr <- (Sigma[,1,2] / sqrt(Sigma[,1,1] * Sigma[,2,2]))
hist(corr, xlim=c(-1,1), xlab="Concurrent Correlation", cex.lab=1.5, cex.axis=1.25, main="")

## @knitr residuals_8.4
library(portes)
r <- ibmsp2608[,2:3]
mu_hat <- extract(fit_VAR_QR, pars=c("mu"))[[1]] %>% colMeans()
phi_hat <- extract(fit_VAR_QR, pars=c("phi"))[[1]] %>% apply(., 2, function(x)colMeans(x))  
a_hat <- matrix(nrow = nrow(fit_VAR_QR), ncol = 2)
for(t in 6:nrow(fit_VAR_QR)) {
  mu_tmp <- mu_hat
  for(i in 1:5)
    mu_tmp <- mu_tmp + (matrix(phi_hat[,i], ncol=2) %>% t) %*% unlist(r[t-i,])
  a_hat[t,] <- unlist(r[t,]) - mu_tmp  
}
LjungBox(a_hat %>% na.omit, lags = 1:8)

## @knitr forecast_8.4
pars <- rstan::extract(fit_VAR_QR)
h <- as.matrix(tail(r, 6))
l <- 6
q1 <- matrix(0, nrow=l, ncol=3)
q2 <- matrix(0, nrow=l, ncol=3)
for(i in 1:l){
  tmp <- matrix(0, nrow=nrow(pars$mu), ncol=ncol(pars$mu))
  for(p in 1:5)
    tmp <- tmp + t(sapply(1:dim(pars$phi[,p,,])[1], function(x)  pars$phi[x,p,,] %*% h[nrow(h)-p+1,]))
  tmp <- sapply(1:nrow(tmp), function(x) pars$mu[x,] + tmp[x,])
  h <- rbind(h, t(rowMeans(tmp)))
  q1[i,] <- quantile(tmp[1,], probs=c(0.05, 0.5, 0.95))
  q2[i,] <- quantile(tmp[2,], probs=c(0.05, 0.5, 0.95)) 
}

## @knitr fit_8.5
ibmsp2608 <- read.table("data/m-ibmsp2608.txt", header = T)
r <- ibmsp2608[,2:3] * 100
VMA <- stan_model("models/Chapter_8/VMA.stan")
fit_VMA <- readRDS("models/Chapter_8/saved_RDS/fit_VMA.rds")

## @knitr print_8.5
print(fit_VMA, pars=c("mu", "Sigma", "Theta"))

## @knitr plot_8.6
r <- read.table("data/m-gs1n3-5301.txt", header = F)
colnames(r) <- c("1-year", "3-year", "date")
par(mfrow=c(2,1), mar=c(2,5,1,4))
plot(log(r[,1]), type="l", ylab="1-year", xlab="Date", cex.lab=2, cex.axis=1.5)
plot(log(r[,2]), type="l", ylab="3-year", xlab="Date", cex.lab=2, cex.axis=1.5)

## @knitr fit_8.6
VARMA <- stan_model("models/Chapter_8/VARMA.stan")
fit_VARMA <- readRDS("models/Chapter_8/saved_RDS/fit_VARMA.rds")

## @knitr print_8.6
print(fit_VARMA, pars=c("mu", "phi", "theta", "Sigma", "Omega", "tau"))

## @knitr residuals_8.6
err <- extract(fit_VARMA, pars=c("err"))[[1]] %>% colMeans()
colnames(err) <- c("1-year", "3-year")
plot.ts(err[-1,], main="Residuals", cex.axis=1.5, cex.lab=1.5, xlab="")
LjungBox(err[-1,], lags = 1:8)

## @knitr corr_8.6
Omega <- extract(fit_VARMA, pars=c("Omega"))[[1]] 
hist(Omega[,1,2], xlim=c(-1,1), xlab="Concurrent Correlation", cex.lab=1.5, cex.axis=1.25, main="")


## @knitr coint_plot
library(mvtnorm)
set.seed(1986)
N <- 200
Phi <- matrix(c(0.5, 0.5, 0.25, 0.5), nrow=2, byrow=TRUE)
Sigma <- matrix(c(3.58,2.5,2.5,2.19), nrow=2)
a <- rmvnorm(N, c(0,0), Sigma)
r <- matrix(nrow = N, ncol = 2)
r[1,] <- 0
for(t in 2:N)
  r[t,] <- Phi %*% r[t-1,] + a[t,]
par(mfrow=c(2,1))
matplot(r, type="l", ylab="")
matplot(apply(r, 2, cumsum), type="l", ylab="")
Acf(apply(r, 2, cumsum))


## @knitr plot_8.6.5
treasury <- read.table("data/w-tb3n6ms.txt", header = T)
r <- treasury[,1:2]
plot.ts(r, cex.axis=1.5, cex.lab=1.5, xlab="", main="")

## @knitr fit_8.6.5
ECM <- stan_model("models/Chapter_8/ECM.stan")
fit_ECM_treasury <- readRDS("models/Chapter_8/saved_RDS/fit_ECM_treasury.rds")

## @knitr print_8.6.5
print(fit_ECM_treasury)
library(tsDyn)
summary(VECM(r, 2, estim="ML"))

## @knitr residuals_8.6.5
N <- nrow(r)
r <- as.matrix(r)
pars <- extract(fit_ECM_treasury)
c0 <- colMeans(pars$c0)
alpha <- colMeans(pars$alpha)
beta <- mean(pars$beta)
phi <- apply(pars$phi, c(2,3,4), mean)
residuals <- matrix(nrow = N, ncol = 2)
residuals[1:3,] <- 0
for(t in 4:N)
  residuals[t,] <-  c0 + (alpha %*% (c(1,-beta) %*% r[t-1,])) + phi[,,1] %*% (r[t-1,] - r[t-2,]) + phi[,,2] %*% (r[t-2,] - r[t-3,]) 
plot.ts(residuals, cex.axis=1.5, cex.lab=1.5, xlab="")
LjungBox(residuals, lags = 1:8)

## @knitr cointegration_8.6.5
r_hat <- matrix(nrow = N, ncol = 2)
r_hat[1:3,] <- 0
for(t in 4:N)
  r_hat[t,] <- r[t-1,] + c0 + (alpha %*% (c(1,-beta) %*% r[t-1,])) + phi[,,1] %*% (r[t-1,] - r[t-2,]) + phi[,,2] %*% (r[t-2,] - r[t-3,]) 
plot(r_hat[,1] - r_hat[,2] , type="l", cex.axis=1.5, cex.lab=1.5, xlab="", ylab="Cointegrating residuals")

## @knitr plot_8.8.3
BHP <- read.table("data/d-bhp0206.txt", header=T)
VALE <- read.table("data/d-vale0206.txt", header=T)
Stocks <- data.frame(BHP=BHP$adjclose %>% log, VALE=VALE$adjclose %>% log)
matplot(Stocks, ylab="Log price", type="l", col=c("black", "red"), cex.axis=1.5, cex.lab=1.5, xlab=""); 
legend('topleft',legend = colnames(Stocks), col = c('black','red'), lty=1:2, cex=1.5)

## @knitr test_8.8.3
library(tseries)
fit <- lm(BHP ~ VALE, Stocks)
summary(fit)
plot(fit$residuals, type="l", cex.axis=1.5, cex.lab=1.5, xlab="")
Acf(fit$residuals)
adf.test(fit$residuals, k=2)


## @knitr fit_8.8.3
fit_ECM_pairs <- readRDS("models/Chapter_8/saved_RDS/fit_ECM_pairs.rds")


## @knitr print_8.8.3
beta <- rstan::extract(fit_ECM_pairs, pars="beta")[[1]]
alpha <- rstan::extract(fit_ECM_pairs, pars="alpha")[[1]]
c0 <- rstan::extract(fit_ECM_pairs, pars="c0")[[1]]
m <- t(apply(Stocks, 1, function(x) quantile((-x[1] + x[2] * beta ) * alpha[,1] - c0[,1] , prob=c(0.05,0.5,0.95))))
m <- m - mean(m[,2]) # I have to center the spread, it looks like I a missing something from the previous calculation
m <- as.data.frame(m)
m$Date <- BHP %>% unite(Mon, day, year, sep=".",col="Date") %>% mutate(Date=as.Date(Date, format="%m.%d.%Y")) %>% pull(Date)
theme_set(theme_bw())
ggplot(m) + geom_ribbon(aes(Date, ymin=`5%`, ymax=`95%`), alpha=0.5, fill="blue") + 
  geom_line(aes(Date, `50%`)) + geom_hline(yintercept = c(sd(m[,2]), 0, -sd(m[,2]))) +
  xlab("") + ylab("Spread") + theme(text=element_text(size=24))

