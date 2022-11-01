{
  require(tidyverse)
  require(rstan)
  require(mvtnorm)
  require(posterior)
  require(astsa)
  rstan_options(auto_write = TRUE)
  rstan_options(javascript = FALSE)
  theme_set(theme_classic(base_size = 24))
  
}


{
library(mvtnorm)
N <- 132   
x <- matrix(NA, nrow=N, ncol=4)
y <- matrix(NA, nrow=N, ncol=1)
alpha1 <- 1.5
alpha2 <- -0.6
beta0 <- 0.25
beta1 <- -0.12
sigma1 <- 0.025
sigma2 <- 0.1
sigmav <- 0.002
Z <- matrix( c(alpha1, alpha2, 0, 0,
               1, 0, 0, 0,
               0, 0, beta1, 0,
               0, 0, 0, 1), nrow=4, byrow=TRUE  )
A <- list()
A[[1]] <- matrix( c(1, 0, 0, 1), nrow=1)
A[[2]] <- matrix( c(1, 0, 1, 1), nrow=1)
  
Q <- diag(c(sigma1^2, 0, sigma2^2, 0))
R <- sigmav
x[1,] <- c(0.3, 0.3, 0.2, -0.1)
y[1,] <- x[1,1] + x[1,3]
s <- rep(1, N)
theta <- matrix( c(0.75, 0.25, 0.25, 0.75), nrow=2, byrow=TRUE)
for(i in 2:N) {
  s[i] <- sample(c(s[i-1], (s[i-1] %% 2) + 1), size=1, prob=theta[s[i-1],])
  x[i, ] <- Z %*% x[i-1, ] + c(0, 0, beta0, 0) + t(rmvnorm(1, c(0, 0, 0, 0), Q))
  y[i, ] <- A[[s[i]]] %*% x[i, ] +  rnorm(1, 0, R)
  
}
par(mfrow=c(2,1), mar=c(2,2,1,0.5))
plot(y, col=s, pch=16); lines(y, lwd=0.5)
matplot(x[,-2], type="l")
}

m <- cmdstan_model("models/Misc/flu.stan")
fit_ <- m$sample(data=list(N=length(flu), y=c(0, diff(as.vector(flu))), m0=rep(0, 4), P0=diag(c(1,1,1,1))), parallel_chains = 4, cores = 4,iter_warmup = 250, iter_sampling = 250)
fit <- fit1
s <- fit$draws(variables = "states") %>% merge_chains() %>% colMeans() %>% matrix(nrow=nrow(y)) 
k <- 6; y = as.matrix(flu); num = length(y); nstate = 4;
Time <- as.matrix(time(flu))
regime <- ifelse(s[,1] < s[,2], 1, 2);
par(mfrow=c(3,1), mar=c(2,3,1,1)+.1)
plot(Time, y, type="n", ylab="")
grid(lty=2); lines(Time, y, col=gray(.7))
text(Time, t(y), col=regime, labels=regime, cex=2)
text(1979,.95,"(a)")
a <- fit$draws(variables = "m") %>% merge_chains() %>% colMeans() %>% as.vector
m1 <- matrix(a[seq(1, 1056, 2)], ncol=4)
m2 <- matrix(a[seq(2, 1056, 2)], ncol=4)
p1 <- matrix(rep(s[,1], 4), ncol=4)
p2 <- matrix(rep(s[,2], 4), ncol=4)
m3 <- m1 * p1 + m2 * p2
matplot(m3)
mu <- fit$draws(variables = "mu") %>% merge_chains() %>% colMeans()
mu1 <- mu[seq(1, 264, 2)]
mu2 <- mu[seq(2, 264, 2)]
mu3 <- mu1 * s[,1] + mu2 * s[,2]
plot(mu3, type="o")
points(y)

## @knitr load_flu
require(tidyverse)
require(rstan)
require(mvtnorm)
require(cmdstanr)
require(posterior)
library(astsa)

## @knitr plot_flu
plot.default(flu, ylab="Mortality", xlab="", cex.axis=1.5, cex.lab=1.5, type="o")

## @knitr fit_flu
# m <- cmdstan_model("models/Misc/flu.stan")
# fit <- m$sample(data=list(N=length(flu), y=as.vector(flu), m0=rep(0, 4), P0=diag(c(1,1,1,1))), parallel_chains = 4, cores = 4,iter_warmup = 250, iter_sampling = 250)
fit <- readRDS("models/Misc/flu.RDS")

## @knitr print_flu
fit$print(max_rows = 12, digits = 4)

## @knitr regimes_flu
k <- 6; y <- as.matrix(flu); num <- length(y); nstate <- 4;
s <- fit$draws(variables = "regimes") %>% merge_chains() %>% colMeans() %>% matrix(nrow=nrow(y)) 
Time <- as.matrix(time(flu))
regime <- ifelse(s[,1] < s[,2], 1, 2);
plot(Time, y, type="n", ylab="", cex=2, cex.axis=2, xlab="")
grid(lty=2); lines(Time, y, col=gray(.7))
text(Time, t(y), col=regime, labels=regime, cex=2)

## @knitr states_flu
a <- fit$draws(variables = "m") %>% merge_chains() %>% colMeans()
m1 <- matrix(a[seq(1, 1056, 2)], ncol=4)
m2 <- matrix(a[seq(2, 1056, 2)], ncol=4)
p1 <- matrix(rep(s[,1], 4), ncol=4)
p2 <- matrix(rep(s[,2], 4), ncol=4)
m3 <- m1 * p1 + m2 * p2
matplot(Time, m3[,-2], type="o", cex=1.5, ylab="", cex.axis=2)

## @knitr mu_flu
mu <- fit$draws(variables = "mu") %>% merge_chains() %>% colMeans()
S <- fit$draws(variables = "S") %>% merge_chains() %>% colMeans()
mu1 <- mu[seq(1, 264, 2)]
mu2 <- mu[seq(2, 264, 2)]
mu3 <- mu1 * s[,1] + mu2 * s[,2]
S1 <- S[seq(1, 264, 2)]
S2 <- S[seq(2, 264, 2)]
S3 <- 2*sqrt(S1 * s[,1] + S2 * s[,2])
plot(Time, mu3, type="n", ylab="", xlab="", ylim=c(0,1), cex=2, cex.axis=2)
grid(lty=2);
points(Time, as.vector(y), pch=16, cex=1.5)
xx = c(Time, rev(Time))
yy = c(mu2-S3, rev(mu3+S3))
polygon(xx, yy, border=8, col=gray(.6, alpha=.3))

