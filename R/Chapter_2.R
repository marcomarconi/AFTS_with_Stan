source('common.R')

## @knitr init_stan
AR_model <- stan_model("../models/AR.stan")

## @knitr plot_gnp

gnp <- scan(file="../data/dgnp82.txt")
gnp1 <- ts(gnp,frequency=4,start=c(1947,2))
plot(gnp1)
points(gnp1)

## @knitr fit_AR_gnp

#AR_model <- stan_model("../models/Chapter_2/AR.stan")
AR_model <- readRDS("../models/Chapter_2/AR.rds")
AR_gnp <- sampling(AR_model, data=list(N=length(gnp1), y=gnp1, K=3), chains=4, cores=4)

## @knitr print_AR_gnp

AR_gnp %>% print(pars=names(AR_gnp)[1:5],digits=4)

## @knitr print_Box.text

y_hat <-  colMeans(extract(AR_gnp, pars="y_hat")[[1]])
(gnp - y_hat)  %>%  Box.test(lag=12, type="Ljung-Box")

## @knitr AR_predictions

# refit the model omitting the last 12 observations
AR_gnp_12 <- sampling(AR_model, data=list(N=length(gnp1)-12, y=gnp1[1:(length(gnp1)-12)], K=3), chains=4, cores=4)
# extract the estimated parameters and calculate the forecasting
pars <- AR_gnp_12 %>% as.matrix()
r_hat <- matrix(nrow=nrow(pars), ncol=13)
r_hat[,1] <- gnp1[length(gnp1)-12] # set the last known observation as first value
for(i in 1:12){
  r_hat[,i+1] <- pars[,1] + rowSums(sapply(1:3, function(p) pars[,2:4][,p] * r_hat[,i] )) # see pag. 56 of the book for this formula
}
# get mean and quantiles of the predictions
r_hat_df <- data.frame(m=colMeans(r_hat), q1=apply(r_hat, 2, function(x)quantile(x, probs=0.025)), q2=apply(r_hat, 2, function(x)quantile(x, probs=0.975)))
ggplot(data.frame(gnp1=gnp1)) + geom_line(aes(x=1:length(gnp1), y=gnp1)) + geom_point(aes(x=1:length(gnp1), y=gnp1)) + geom_ribbon(data=r_hat_df, aes(x=(length(gnp1)-12):length(gnp1), ymin=q1, ymax=q2 ), alpha=0.5, col="blue", fill="blue")  + xlab("Time")


## @knitr AR_loo_compare

AR1_gnp <- sampling(AR_model, data=list(N=length(gnp1), y=gnp1, K=1), chains=4, cores=4)
loo_compare(loo(AR_gnp), loo(AR1_gnp))

## @knitr fit_AR_GAS

uic <- scan("../data/UIC.txt")
AR_GAS_model <- stan_model("../models/Chapter_2/AR-GAS.stan")
AR_GAS_model <- readRDS("../models/Chapter_2/AR-GAS.rds")
AR_GAS_uic <- sampling(AR_GAS_model, data=list(N=length(uic), y=uic), chains=4, cores=4)

## @knitr plot_AR_GAS

par(mfrow=c(2,1), mai = c(1, 1, 0.1, 0.1))
plot.ts(uic)
ar <- extract(AR_GAS_uic, pars="ar")[[1]] %>% colMeans()
plot.ts(2*( invlogit(ar)-0.5), ylab="AR")

## @knitr load_TSE
y1 <- read.table("../data/w-gs1yr.txt", header=T)
y3 <- read.table("../data/w-gs3yr.txt", header=T)
y <- ts(cbind(y1$rate, y3$rate), start = c(1962, 1), frequency = 52)
colnames(y) <- c("Treasury1", "Treasury3")
plot.ts(y, plot.type = "single", ylab="Percent")

## @knitr lm_TSE
fit <- lm(Treasury3 ~ Treasury1, data=y)
summary(fit)
plot.ts(fit$residuals)
acf(fit$residuals)

## @knitr lm_diff_TSE
y_diff <- diff(y)
fit_diff <- lm(Treasury3 ~ 0+Treasury1, data=y_diff)
summary(fit_diff)
plot.ts(fit_diff$residuals)
acf(fit_diff$residuals)

## @knitr stan_TSE
TSE <- stan_model("../models/Chapter_2/TSE.stan")
tse <- sampling(TSE, data=list(N=nrow(y_diff), y=y_diff[,2], x=y_diff[,1]), chains=1)

## @knitr plot_stan_TSE
extract(tse, pars="a")[[1]] %>% colMeans() -> a
plot.ts(a)
acf(a)
Box.test(a, lag=12, type="Ljung-Box")
