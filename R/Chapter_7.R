source('common.R')

## @knitr load_7.2
ibm <- read.table("data/d-ibm6298.txt", header=T)
ibm$date <- as.Date(ibm$date %>% as.character(), format="%Y%m%d")
plot(ibm$date,ibm$rtn)

## @knitr fit_7.2
IGARCH <- stan_model("models/Chapter_3/IGARCH.stan")
fit_ibm <- readRDS("models/Chapter_7/saved_RDS/fit_ibm.rds")

## @knitr print_7.2
print(fit_ibm, names(fit_ibm)[1:4])

## @knitr plot_7.2
pars <- extract(fit_ibm)
sd <- pars$sigma
res <- sapply(1:nrow(sd), function(i)  ibm$rtn[i] / sd[,i] )
sd_1 <- sqrt( sd[,ncol(sd)]^2 + pars$alpha1 * sd[,ncol(sd)]^2 * (res[,ncol(res)]^2 - 1))
VaR_0.05 <- (1.65 * (sd_1) * 10e+6) 
VaR_0.01 <- (2.326 * (sd_1) * 10e+6) 
ggplot(rbind(data.frame(VaR = VaR_0.01, Quantile="VaR 99%"), data.frame(VaR = VaR_0.05, Quantile="VaR 5%"))) +
  geom_density(aes(VaR/1e6, fill=Quantile)) + xlab("VaR (millions $)") 

## @knitr load_7.5.3
dates <- split(ibm$date, ceiling(seq_along(ibm$date)/21)) %>% lapply(function(x) head(x,1)) %>% unlist %>% as.Date
mins <- split(ibm$rtn, ceiling(seq_along(ibm$rtn)/21))  %>% lapply(min) %>% unlist 
maxs <- split(ibm$rtn, ceiling(seq_along(ibm$rtn)/21))  %>% lapply(max) %>% unlist 
names(maxs) <- dates
names(mins) <- dates
par(mfrow=c(2,1), mar=c(2,2,2,2))
barplot(maxs, cex.names = 1)
barplot(mins, cex.names = 1)

## @knitr fit_7.5.3
GEV <- stan_model("models/Chapter_7/GEV.stan")
fit_gev_max <- sampling(GEV, data = list(N = length(maxs), y=maxs*100), chains = 4, cores = 4)
fit_gev_min <- sampling(GEV, data = list(N = length(mins), y=-mins*100), chains = 4, cores = 4)
print(fit_gev_max)
print(fit_gev_min)

## @knitr residuals_7.5.3
pars <- extract(fit_gev_max, pars=c("mu", "sigma", "xi")) 
pars_mean <- pars %>%  lapply(mean)
w <- (1 + pars_mean$xi * (maxs*100 - pars_mean$mu) / pars_mean$sigma)^(-1/pars_mean$xi)
par(mfrow=c(1,2))
plot(w, ylab="Residuals")
qqplot(w, rexp(1000), ylab="Exponential quantiles")


## @knitr var_7.5.3
# VaR based on extreme value theory
pars <- extract(fit_gev_min, pars=c("mu", "sigma", "xi")) 
n <- 21
p <- 0.01
VaR_EVT <- (pars$mu - pars$sigma/pars$xi * (1- (-n * log(1-p))^(-pars$xi))) / 100 * 10e+6
# Var based on empirical quantile estimation
nibm=-log(ibm[,2]+1)*100
VaR_quantile <- sapply(1:1000, function(x) quantile(sample(nibm, size = length(nibm), replace = TRUE),  0.99) / 100 * 10e+6 )
# plot all together
df <- rbind(data.frame(VaR = VaR_0.01, Measure="RiskMetrics"), data.frame(VaR = VaR_EVT, Measure="EVT"), data.frame(VaR = VaR_quantile, Measure="Quantiles"))
ggplot(df) + geom_density(aes(VaR/1e6, fill=Measure), alpha=0.75) + xlab("VaR (millions $)") 



## @knitr fit_7.7.4
GPD <- stan_model("models/Chapter_7/GPD.stan")
threshold <- 0.03
y <- (-ibm$rtn)[(-ibm$rtn) > 0.03]
fit_gpd <- sampling(GPD, data = list(ymin = threshold, N = length(y), y=y), chains = 4, cores = 4)
print(fit_gpd, names(fit_gpd)[1:2], digits=5)

## @knitr var_7.7.4
pars <- extract(fit_gpd)
D <- 252
p <- 0.01
# formula (7.36) in the book. Notice that I have to use the sample mean as location.
VaR_GPD <- (mean(y)-((pars$sigma-pars$k*(threshold-mean(y)))/pars$k)*(1-(-D * log(1-p))^(-pars$k)) )  * 10e6
df <- rbind(data.frame(VaR = VaR_0.01, Measure="RiskMetrics"), 
            data.frame(VaR = VaR_EVT, Measure="EVT"), 
            data.frame(VaR = VaR_quantile, Measure="Quantiles"),
            data.frame(VaR = VaR_GPD, Measure="GPD")
            ) %>% mutate(Measure=factor(Measure, levels=c("RiskMetrics", "Quantiles", "EVT", "GPD")))
ggplot(df) + geom_density(aes(VaR/1e6, fill=Measure), alpha=0.75) + xlab("VaR (millions $)") 
