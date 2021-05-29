source('common.R')

## @knitr load_4.3
d_ibmvwewsp6203 <- read.table("data/d-ibmvwewsp6203.txt", header=F )
y <- (d_ibmvwewsp6203$V2 %>% `*`(100) %>% ts) 
plot.ts(y, ylab="Log return", xlab="Time")


## @knitr fit_4.3
#TAR <- stan_model(file = "models/Chapter_4/AR-TAR-GARCH.stan")
#fit_TAR_ibmvwewsp6203 <- vb(TAR, data = list(N = length(y), y = y, K=2,threshold=0, sigma1=sd(y), mu1=mean(y)), iter=500)
#fit_ARGARCH_ibmvwewsp6203 <- stan(file = "models/Chapter_3/AR-GARCH.stan", data = list(N = length(y), y = y, K=2), chains = 1, iter=500)
#fit_TGARCH_ibmvwewsp6203 <- stan(file = "models/Chapter_3/TGARCH.stan", data = list(N = length(y), y = y), chains = 1, iter=500)
fit_TAR_ibmvwewsp6203 <- readRDS("models/Chapter_4/saved_RDS/fit_TAR_ibmvwewsp6203.rds")
fit_ARGARCH_ibmvwewsp6203 <- readRDS("models/Chapter_4/saved_RDS/fit_ARGARCH_ibmvwewsp6203.rds")
fit_TGARCH_ibmvwewsp6203 <- readRDS("models/Chapter_4/saved_RDS/fit_TGARCH_ibmvwewsp6203.rds")

## @knitr means_4.3
pars <- rstan::extract(fit_ARGARCH_ibmvwewsp6203)
mean_ARGARCH <- pars$ar0 / (1 - pars$ar[[1]] - pars$ar[[2]])
pars <- rstan::extract(fit_TGARCH_ibmvwewsp6203)
mean_TGARCH <- pars$mu
pars <- rstan::extract(fit_TAR_ibmvwewsp6203)
mean_TAR <- pars$ar0 / (1 - pars$ar[[1]] - pars$ar[[2]])
plot(density(mean_ARGARCH), xlim=c(0, 0.15), ylim=c(0,35), main="", xlab="", ylab="")
lines(density(mean_TGARCH), col="blue")
lines(density(mean_TAR), col="red")
abline(v=0.039, lty=2)
legend("topright", legend=c("AR-GARCH", "TGARCH", "AR-TAR-GARCH"), col=c("black", "blue", "red"), lty=rep(1,3))


## @knitr fit_4.4

ARCH <- stan_model("models/Chapter_3/ARCH.stan")
STAR <- stan_model("models/Chapter_4/STAR.stan")
m_3m4608<- read.table("data/m-3m4608.txt", header=T )
fit_m_3m4608_ARCH <- sampling(ARCH, data = list(N = length(m_3m4608$rtn), y = m_3m4608$rtn, Kar=2, family=0), iter=1000,chains = 4, cores = 4, refresh=0)
fit_m_3m4608_STAR <- sampling(STAR, data = list(N = length(m_3m4608$rtn), y = m_3m4608$rtn, K=2, sigma1=sd(m_3m4608$rtn),gamma2_mean=1000, gamma2_sd=500 ), iter=1000,chains = 4, cores = 4, refresh=0)
print(fit_m_3m4608_ARCH, pars=names(fit_m_3m4608_ARCH)[1:4])
print(fit_m_3m4608_STAR, pars=names(fit_m_3m4608_STAR)[1:7])

## @knitr curve_4.4
pars <- extract(fit_m_3m4608_STAR)
x <- seq(-0.01, 0.01, 0.0001); 
num <- sapply(1:2000, function(i) pars$gamma0[i] + pars$gamma1[i]*x^2)
den <- sapply(1:2000, function(i) pars$gamma2[i]*x)
res <- (num / (1 + exp(-den)))
q <- apply(res, 1, function(i) quantile(i, probs = c(0.05, 0.5, 0.95)))
plot(NA,xlim=c(-0.01, 0.01),ylim=c(0,0.2e-2),ylab="STAR effect",xlab="Previous Shock")
polygon(c(x, rev(x)), c(q[1,], rev(q[3,])), lty=2, col="gray")
lines(x, q[2,], ylim=c(0,1e-3))

## @knitr loo_4.4
loo_compare(loo(fit_m_3m4608_ARCH), loo(fit_m_3m4608_STAR))


## @knitr fit_4.5
q_gnp4791 <- scan("data/q-gnp4791.txt" )*100
MSA <- stan_model("models/Chapter_4/MSA.stan")
fit_q_gnp4791 <- sampling(MSA, data = list(N = length(q_gnp4791), y=q_gnp4791, K=2, Kar=4), chains = 4, cores = 4, iter=1000, init=0, refresh=0)
plot(fit_q_gnp4791, pars=names(fit_q_gnp4791)[1:14])

## @knitr plot_m_4.5
pars <- extract(fit_q_gnp4791)
m1 <- pars$c[,1] / (1 - (pars$phi[,1,1] + pars$phi[,1,2] + pars$phi[,1,3] + pars$phi[,1,4]))
m2 <- pars$c[,2] / (1 - (pars$phi[,2,1] + pars$phi[,2,2] + pars$phi[,2,3] + pars$phi[,2,4]))
plot(density(m1), xlim=c(-7,7), col="blue")
lines(density(m2), xlim=c(-7,7), col="red")

## @knitr plot_w_4.5
w1 <- 1/pars$p[,1] 
w2 <- 1/pars$p[,2] 
plot(density(w1), xlim=c(0,6),col="blue")
lines(density(w2),xlim=c(0,6), col="red")

## @knitr plot_s_4.5
xi <- pars$xi %>% colMeans()
matplot(xi, type="l", col=c("blue", "red"))
colMeans(xi)
