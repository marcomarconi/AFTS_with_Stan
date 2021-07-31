source('common.R')

## @knitr load_10.1
x <- read.table("data/d-hkjp0608.txt", header=T)
hk <- x$HK %>% log %>% diff %>% `*`(100)
jp <- x$JP %>% log %>% diff %>% `*`(100)
r <- cbind(hk, jp)
plot.ts(r, main="log returns", cex.axis=2, cex.lab=2, cex.main=2)

## @knitr fit_10.1
GARCH <- stan_model("models/Chapter_3/GARCH.stan")
fit_HK <- sampling(GARCH, data = list(N = length(hk), y=hk, sigma1=sd(hk), K=1), chains = 1, cores = 4, refresh=0)
fit_JP <- sampling(GARCH, data = list(N = length(jp), y=jp, sigma1=sd(jp), K=1), chains = 1, cores = 4, refresh=0)
pars1 <- rstan::extract(fit_HK, pars="sigma")[[1]]
pars2 <- rstan::extract(fit_JP, pars="sigma")[[1]]
plot.ts(cbind(hk=colMeans(pars1), jp=colMeans(pars2)), main="", cex.axis=2, cex.lab=2, )


## @knitr EWMA_10.1
EWMA <- stan_model("models/Chapter_10/EWMA.stan")
fit_EWMA <- readRDS("models/Chapter_10/saved_RDS/fit_EWMA.rds")


## @knitr EWMA_plot_10.1
print(fit_EWMA, pars=names(fit_EWMA)[1:3])
Sigma <- rstan::extract(fit_EWMA, pars="Sigma")[[1]]
Sigma <- colMeans(Sigma) %>% apply(c(1,2), mean)
colnames(Sigma) <- c("hk", "jp")
plot.ts(Sigma, main="", cex.axis=2, cex.lab=2, )

## @knitr load_10.2
r <- read.table("data/m-pfemrk6508.txt", header=T)
plot.ts(r[,-1], main="", cex.axis=2, cex.lab=2)

## @knitr fit_10.2
DVEC <- stan_model("models/Chapter_10/DVEC.stan")
fit_DVEC <- readRDS("models/Chapter_10/saved_RDS/fit_DVEC.rds")

## @knitr plot_10.2
Sigma <- rstan::extract(fit_DVEC, pars="Sigma")[[1]]
sigmas <- colMeans(Sigma) %>% apply(c(1,2), mean)
colnames(sigmas) <- c("pfe", "mrk")
plot.ts(sigmas, main="", cex.axis=2, cex.lab=2, )
corr <- apply(Sigma, c(2,3,4), mean) %>% apply(.,1,cov2cor) %>% t 
plot.ts(corr[,2])

## @knitr load_10.4
x <- read.table("data/d-hkjp0608.txt", header=T)
hk <- x$HK %>% log %>% diff %>% `*`(100)
jp <- x$JP %>% log %>% diff %>% `*`(100)
r <- cbind(hk, jp)

## @knitr fit_10.4
MGARCH <- stan_model("models/Chapter_10/MGARCH.stan")
fit_MGARCH_hkjp <- readRDS("models/Chapter_10/saved_RDS/fit_MGARCH_hkjp.rds")

## @knitr print_10.4
print(fit_MGARCH_hkjp, pars=names(fit_MGARCH_hkjp)[1:13])
pars <- rstan::extract(fit_MGARCH_hkjp)
sigmas <- pars$Xi %>% colMeans() %>% sqrt
colnames(sigmas) <- c("hk", "jp")
plot.ts(sigmas, cex.axis=2, cex.lab=2)

## @knitr residuals_10.4
library(portes)
mu <- pars$mu %>% colMeans()
a <- r - mu
residuals <- a / sigmas
LjungBox(residuals, lags = 1:12)

## @knitr load_10.5
x <- read.table("data/m-ibmsp2608.txt", header=T)
r <- x[,2:3]*100

## @knitr fit_10.5
MGARCH_chol <- stan_model("models/Chapter_10/MGARCH_cholesky.stan")
fit_MGARCH_chol <- readRDS("models/Chapter_10/saved_RDS/fit_MGARCH_chol.rds")

## @knitr print_10.5
print(fit_MGARCH_chol, names(fit_MGARCH_chol)[1:13])

## @knitr plot_10.5
library(lubridate)
rho <- rstan::extract(fit_MGARCH_chol, pars="phi")[[1]]
q <- t(apply(rho, 2, function(x) quantile(x, probs=c(0.05, 0.5, 0.95))))
q <- data.frame(Date=as_date(x[,1] %>% as.character()), q)
ggplot(q) + geom_ribbon(aes(Date, ymin=`X5.`, ymax=`X95.`),  color="lightslateblue", fill="lightblue", size=0.2) + geom_line(aes(Date, `X50.`)) + ylab("Time-varying correlation") + xlab("") + ylim(c(0,1)) + theme(text=element_text(size=32))

## @knitr load_10.6
x <- read.table("data/d-fxsk9904.txt")*100
x <- as.matrix(x); colnames(x) <- c("USEU", "USJP", "IBM", "Dell")
plot.ts(x, cex.axis=2, cex.lab=2)

## @knitr Psi_10.6
x <- read.table("data/d-fxsk9904.txt")*100
x <- as.matrix(x); colnames(x) <- c("USEU", "USJP", "IBM", "Dell")
N<-nrow(x)
K <- 4
m <- 69
Psi <- array(dim = c(N,K,K))
Psi[1,,] <- diag(rep(1, K))
a <- matrix(nrow=N, ncol=K)
a[1,] <- rep(0, K)
for(t in 2:N) {
  a[t,] <- x[t,]
  if(t > m)
    Psi[t,,] <- cor(a[(t-m):(t-1),])
  else
    Psi[t,,] <- diag(rep(1, K))
}

## @knitr fit_10.6
DCC <- stan_model("models/Chapter_10/DCC.stan")
fit_DCC <- readRDS("models/Chapter_10/saved_RDS/fit_DCC.rds")

## @knitr print_10.6
print(fit_DCC, names(fit_DCC)[1:15])

## @knitr D_10.6
library(roll)
x <- read.table("data/d-fxsk9904.txt")*100
x <- as.matrix(x); colnames(x) <- c("USEU", "USJP", "IBM", "Dell")
rolling_var <- roll_sd(x, 69)^2
D <- rstan::extract(fit_DCC, pars="D")[[1]]; 
sigmas <- apply(D, c(2,3,4), mean) 
sigmas <- apply(sigmas, 1, diag)^2 %>% t; colnames(sigmas) <- c("USEU", "USJP", "IBM", "Dell")
par(mfrow=c(4,1), mar=c(2,5,1,4))
plot(sigmas[,1], cex.axis=2, cex.lab=2, main="", xlab="", ylab=colnames(sigmas)[1], type="l")
lines(rolling_var[,1], lty=2)
plot(sigmas[,2], cex.axis=2, cex.lab=2, main="", xlab="", ylab=colnames(sigmas)[2], type="l")
lines(rolling_var[,2], lty=2)
plot(sigmas[,3], cex.axis=2, cex.lab=2, main="", xlab="", ylab=colnames(sigmas)[3], type="l")
lines(rolling_var[,3], lty=2)
plot(sigmas[,4], cex.axis=2, cex.lab=2, main="", xlab="", ylab=colnames(sigmas)[4], type="l")
lines(rolling_var[,4], lty=2)

## @knitr Rho_10.6
rolling_cor <- roll_cor(x, x, 69)
Rho_m <- rstan::extract(fit_DCC, pars="Rho_m")[[1]]; 
USEU_vs_JPUS <- Rho_m[1,,1,2]
IBM_vs_USEU <- Rho_m[1,,1,3]
IBM_vs_JPUS <- Rho_m[1,,2,3]
Dell_vs_USEU <- Rho_m[1,,1,4]
Dell_vs_JPUS <- Rho_m[1,,2,4]
Dell_vs_IBM <- Rho_m[1,,3,4]
par(mfrow=c(6,1), mar=c(2,5,1,4))
plot(USEU_vs_JPUS, cex.axis=2, cex.lab=2, main="", xlab="", ylab="USEU_vs_JPUS", ylim=c(-0.75,0.75), type="l")
lines(rolling_cor[1,2,], lty=2)
plot(IBM_vs_USEU, cex.axis=2, cex.lab=2, main="", xlab="", ylab="IBM_vs_USEU", ylim=c(-0.75,0.75), type="l")
lines(rolling_cor[1,3,], lty=2)
plot(IBM_vs_JPUS, cex.axis=2, cex.lab=2, main="", xlab="", ylab="IBM_vs_JPUS", ylim=c(-0.75,0.75), type="l")
lines(rolling_cor[2,3,], lty=2)
plot(Dell_vs_USEU, cex.axis=2, cex.lab=2, main="", xlab="", ylab="Dell_vs_USEU", ylim=c(-0.75,0.75), type="l")
lines(rolling_cor[1,4,], lty=2)
plot(Dell_vs_JPUS, cex.axis=2, cex.lab=2, main="", xlab="", ylab="Dell_vs_JPUS", ylim=c(-0.75,0.75), type="l")
lines(rolling_cor[2,4,], lty=2)
plot(Dell_vs_IBM, cex.axis=2, cex.lab=2, main="", xlab="", ylab="Dell_vs_IBM", ylim=c(-0.75,0.75), type="l")
lines(rolling_cor[3,4,], lty=2)


