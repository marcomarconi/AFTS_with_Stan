source('common.R')

## @knitr load_5.2
ibm91_ads <- read.table("data/ibm91-ads.dat", header=F )
colnames(ibm91_ads) <- c("A", "D", "S")
plot(ibm91_ads$S*ibm91_ads$D, xlab="Time", ylab="Price Change")

## @knitr fit_5.2
ADS <- stan_model("models/Chapter_5/ADS.stan")
fit_ADS <- readRDS("models/Chapter_5/saved_RDS/fit_ADS.rds")

## @knitr print_5.2
print(fit_ADS)

## @knitr plot_A_5.2
pars <- rstan::extract(fit_ADS)
pA1_0 <- invlogit(pars$beta0 + 0*pars$beta1)
pA1_1 <- invlogit(pars$beta0 + 1*pars$beta1)
plot(density(pA1_0), xlim=c(0.2,0.5), col="red", xlab="Probability", main="")
lines(density(pA1_1), xlim=c(0.2,0.5), col="blue")
legend("topright", legend=c("P(Ai = 1 | Ai-1 = 0)", "P(Ai = 1 | Ai-1 = 1)"), col=c("red", "blue"), lty=1:1, lwd=2, cex=1.2)

## @knitr plot_D_5.2
pD1_D0_A0 <- invlogit(pars$gamma0 + 0*pars$gamma1)
pD1_D1_A1 <- invlogit(pars$gamma0 + 1*pars$gamma1)
pD1_D_1_A1 <- invlogit(pars$gamma0 + -1*pars$gamma1)
plot(density(pD1_D0_A0), xlim=c(0,1), ylim=c(0, 200), col="red", xlab="Probability", main="")
lines(density(pD1_D1_A1), xlim=c(0,1), col="blue")
lines(density(pD1_D_1_A1), xlim=c(0,1), col="green")
legend("topright", legend=c("P(Di = 1 | Di-1 = 0)", "P(Di = 1 | Di-1 = 1)", "P(Di = 1 | Di-1 = -1)"), col=c("red", "blue", "green"), lty=1:1, lwd=2, cex=1.2)

## @knitr plot_S_5.2
pS <- pgeom(0, prob=invlogit(mean(pars$phi0_u) + (1:10)*mean(pars$phi1_u))) 
par(mar=c(5,5,5,5))
plot(pS, ylab = "Probality of Price increase by 1", xlab=expression(S[t-1]), ylim=c(0,1),cex.lab=1.5, cex.axis=1.5)

## @knitr load_5.3
ibm <- scan("data/ibm1to5-dur.txt")
par(mfrow=c(1,3))
plot.ts(ibm, ylab="Adjusted Duration", main="IBM trading duration")
hist(ibm)
Acf(ibm, ylab="Adjusted Duration", ylim=c(-0.1,0.2))

## @knitr fit_5.3
ACD <- stan_model("models/Chapter_5/ACD.stan")
fit_ACD_ibm <- readRDS("models/Chapter_5/saved_RDS/fit_ACD_ibm.rds")

## @knitr print_5.3
print(fit_ACD_ibm, pars=names(fit_ACD_ibm)[1:4])

## @knitr plot_5.3
pars <- rstan::extract(fit_ACD_ibm)
residuals <- ibm / colMeans(pars$psi)
Acf(residuals)

## @knitr test_1_5.4
Tsay(residuals, p = 4)

## @knitr fit_5.4
TARACD <- stan_model("models/Chapter_5/TAR-ACD.stan")
fit_TARACD_ibm <- readRDS("models/Chapter_5/saved_RDS/fit_TARACD_ibm.rds")

## @knitr print_5.4
print(fit_TARACD_ibm, pars=names(fit_TARACD_ibm)[1:8])

## @knitr test_2_5.4
pars <- rstan::extract(fit_TARACD_ibm)
residuals_nl <- ibm / colMeans(pars$psi)
Tsay(residuals, p = 4)
Tsay(residuals_nl, p = 4)
loo_compare(loo(fit_ACD_ibm), loo(fit_TARACD_ibm))

## @knitr load_5.5
day15 <- read.csv("data/day15v.dat", header = F, sep="\t")
colnames(day15) <- c("Day", "Time", "TT", "D", "S", "N", "I", "C")
par(mfrow=c(2,2))
hist(log(day15$TT)); hist(day15$S); hist(day15$D); hist(day15$N)

## @knitr fit_5.5
PCD <- stan_model("models/Chapter_5/PCD.stan")
fit_PCD <- readRDS("models/Chapter_5/saved_RDS/fit_PCD.rds")

## @knitr print_5.5
print(fit_PCD)

