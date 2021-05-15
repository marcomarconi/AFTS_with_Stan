source('common.R')

## @knitr load_intel
intel <- read.table("data/m-intc7308.txt", header=T)
intel <- xts(intel$rtn, as.Date(as.character(intel$date), format = "%Y%m%d")) 
plot(intel)

## @knitr fit_ARCH_normal
ARCH <- stan_model("models/Chapter_3/ARCH.stan")
fit_intel_3 <- sampling(ARCH, data = list(N = length(intel), y = intel %>% as.vector(), Kar = 3, family = 0), chains = 4, refresh=0)
fit_intel_1 <- sampling(ARCH, data = list(N = length(intel), y = intel %>% as.vector(), Kar = 1, family = 0), chains = 4, refresh=0)
print(fit_intel_3, pars=names(fit_intel_3)[1:5])
print(fit_intel_1, pars=names(fit_intel_1)[1:3])


## @knitr plot_ARCH_residuals
sq <- extract(fit_intel_1, pars="sigma")[[1]] %>% colMeans()
mu <- extract(fit_intel_1, pars="mu")[[1]] %>% mean()
residuals <- ((intel - mu) / sq)
plot(residuals)
Box.test(residuals, lag=10, type="Ljung-Box")

## @knitr ARCH_compare
a <- extract_log_lik(fit_intel_3, merge_chains = FALSE)[,,4:(length(intel))]
b <- extract_log_lik(fit_intel_1, merge_chains = FALSE)[,,4:(length(intel))]
r_eff_a <- relative_eff(exp(a))
r_eff_b <- relative_eff(exp(b))
loo_compare(loo(a, r_eff = r_eff_a, moment_match = TRUE), loo(b, r_eff = r_eff_b, moment_match = TRUE))


## @knitr fit_ARCH_t
fit_intel_3_t <- sampling(ARCH, data = list(N = length(intel), y = intel %>% as.vector(), Kar = 3, family = 1), chains = 4, refresh=0)
print(fit_intel_3_t, pars=names(fit_intel_3_t)[1:5])

## @knitr plot_ARCH_t_residuals
sq <- extract(fit_intel_3_t, pars="sigma")[[1]] %>% colMeans()
mu <- extract(fit_intel_3_t, pars="mu")[[1]] %>% mean()
residuals <- ((intel - mu) / sq)
plot(residuals)
Box.test(residuals, lag=12, type="Ljung-Box")

## @knitr ARCH_t_compare
a <- extract_log_lik(fit_intel_3, merge_chains = FALSE)[,,4:(length(intel))]
b <- extract_log_lik(fit_intel_3_t, merge_chains = FALSE)[,,4:(length(intel))]
r_eff_a <- relative_eff(exp(a))
r_eff_b <- relative_eff(exp(b))
loo_compare(loo(a, r_eff = r_eff_a, moment_match = TRUE), loo(b, r_eff = r_eff_b, moment_match = TRUE))


## @knitr plot_S&P500
sp500 <- scan("data/sp500.dat.txt")
plot.ts(sp500)
Acf(sp500)
Pacf(sp500)

## @knitr fit_GARCH_sp500
GARCH <- stan_model("models/Chapter_3/GARCH.stan")
ARGARCH <- stan_model("models/Chapter_3/AR-GARCH.stan")
fit_garch_sp500 <- sampling(GARCH, data = list(N = length(sp500), y = sp500, K = 1, sigma1 = sd(sp500)), chains = 4, refresh=0)
fit_argarch_sp500 <- sampling(ARGARCH, data = list(N = length(sp500), y = sp500, K = 3), chains = 4, refresh=0)
print(fit_garch_sp500, pars=names(fit_garch_sp500)[1:4], digits=5)
print(fit_argarch_sp500, pars=names(fit_argarch_sp500)[1:8], digits=5)

## @knitr loo_GARCH_sp500
a <- extract_log_lik(fit_garch_sp500, merge_chains = FALSE)[,,4:(length(sp500))]
b <- extract_log_lik(fit_argarch_sp500, merge_chains = FALSE)[,,4:(length(sp500))]
r_eff_a <- relative_eff(exp(a))
r_eff_b <- relative_eff(exp(b))
loo_compare(loo(a, r_eff = r_eff_a, moment_match = TRUE), loo(b, r_eff = r_eff_b, moment_match = TRUE))


## @knitr plot_GARCH_sigma_sp500
extract(fit_garch_sp500, pars="sigma")[[1]] %>% colMeans()  %>% plot(type="l")

## @knitr plot_IBM
ibm_sp500 <- read.table("data/m-ibmsplnsu.dat", header=T)
plot.ts(ibm_sp500$ibm)

## @knitr fit_IBM_garch
fit_argarch_ibm2 <- sampling(ARGARCH, data = list(N = length(ibm_sp500$ibm), y = ibm_sp500$ibm, K = 1), chains = 1, refresh=0)


## @knitr fit_IBM_garch_I
ARGARCH_I <- stan_model("models/Chapter_3/AR-GARCH-I.stan")
fit_argarch_I_ibm2 <- sampling(ARGARCH_I, data = list(N = length(ibm_sp500$ibm), y = ibm_sp500$ibm, u = ibm_sp500$summer, K = 1), chains = 1, refresh=0)
plot(fit_argarch_I_ibm2, pars=names(fit_argarch_I_ibm2)[1:6])

## @knitr loo_IBM
loo_compare(loo(fit_argarch_ibm2, moment_match = TRUE), loo(fit_argarch_I_ibm2, moment_match = TRUE))

## @knitr fit_SP_garch_x
GARCH <- stan_model("models/Chapter_3/GARCH.stan")
fit_garch_sp500 <- sampling(GARCH , data=list(N = length(ibm_sp500$sp), y = ibm_sp500$sp, K = 2, sigma1=sd(ibm_sp500$sp)), chains = 4, refresh=0)
GARCH_x <- stan_model("models/Chapter_3/GARCH_x.stan")
fit_garch_x_sp500 <- sampling(GARCH_x , data=list(N = length(ibm_sp500$sp), y = ibm_sp500$sp, x=ibm_sp500$ibm, K = 2), chains = 4, refresh=0)


## @knitr compare_SP_garch_x
print(fit_garch_x_sp500, pars=names(fit_garch_x_sp500)[1:6])
a <- extract_log_lik(fit_garch_sp500, merge_chains = FALSE)[,,4:(nrow(ibm_sp500))]
b <- extract_log_lik(fit_garch_x_sp500, merge_chains = FALSE)[,,4:(nrow(ibm_sp500))]
r_eff_a <- relative_eff(exp(a))
r_eff_b <- relative_eff(exp(b))
loo_compare(loo(a, r_eff = r_eff_a, moment_match = TRUE), loo(b, r_eff = r_eff_b, moment_match = TRUE))
