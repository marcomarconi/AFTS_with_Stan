source('common.R')

## @knitr load_9.2.1
m <- read.table("data/m-fac9003.txt", header=T)

## @knitr fit_9.2.1
SFM <- stan_model("models/Chapter_9/SFM.stan")
fit_SFM <- readRDS("models/Chapter_9/saved_RDS/fit_SFM.rds")

## @knitr plot_9.2.1
pars <- rstan::extract(fit_SFM)
cbind(beta=colMeans(pars$beta), sd=colMeans(pars$tau))
q <- t(apply(pars$beta, 2, function(x) quantile(x,probs=c(0.05,0.5,0.95))))
rownames(q) <- colnames(m[,1:(ncol(m)-1)])
ggplot(as.data.frame(q)) + geom_bar(aes( x=rownames(q), y=`50%`),  stat = "identity") + geom_errorbar(aes(rownames(q), ymin=`5%`, ymax=`95%`, width=0.1)) + xlab("")+ ylab("beta")

## @knitr corr_9.2.1
library(corrplot)
rho <- apply(pars$Omega, c(2,3), mean)
colnames(rho) <- rownames(rho) <- colnames(m[,1:(ncol(m)-1)])
corrplot(rho, tl.cex = 2,  cl.cex = 1.5)

## @knitr load_9.3.1
da <- read.table("data//m-barra-9003.txt", header=T)
rm = matrix(apply(da,2,mean),1)
rtn = da - matrix(1,168,1)%*%rm # mean-correct the returns
finance <- c(rep(1,4),rep(0,6))
technology <- c(rep(0,4),rep(1,3),rep(0,3))
other <- c(rep(0,7),rep(1,3))
ind.dum <- cbind(Finance=finance,Technology=technology,Other=other)
print(ind.dum)

## @knitr fit_9.3.1
BARRA <- stan_model("models/Chapter_9/BARRA.stan")
fit_BARRA <- readRDS("models/Chapter_9/saved_RDS/fit_BARRA.rds")

## @knitr plot_9.3.1
extract(fit_BARRA, pars="eff")[[1]] %>% colMeans() %>% plot.ts

