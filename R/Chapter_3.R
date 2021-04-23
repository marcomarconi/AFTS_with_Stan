source('common.R')

## @knitr load_intel
intel <- read.table("data/m-intc7308.txt", header=T)
intel <- xts(intel$rtn, as.Date(as.character(intel$date), format = "%Y%m%d")) 
plot(intel)

## @knitr fit_ARCH_normal
ARCH <- stan_model("models/Chapter_3/ARCH.stan")
fit_intel_3 <- sampling(ARCH, data = list(N = length(intel), y = intel %>% as.vector(), Kar = 3, family = 0), chains = 1, refresh=0)
fit_intel_1 <- sampling(ARCH, data = list(N = length(intel), y = intel %>% as.vector(), Kar = 1, family = 0), chains = 1, refresh=0)
print(fit_intel_3, pars=names(fit_intel_3)[1:5])
print(fit_intel_1, pars=names(fit_intel_1)[1:3])


## @knitr plot_ARCH_residuals
sq <- extract(fit_intel_1, pars="sigma")[[1]] %>% colMeans()
mu <- extract(fit_intel_1, pars="mu")[[1]] %>% mean()
residuals <- ((intel - mu) / sq)
plot(residuals)
Box.test(residuals, lag=10, type="Ljung-Box")

## @knitr ARCH_compare
loo_compare(loo(fit_intel_3), loo(fit_intel_1))


## @knitr fit_ARCH_t
fit_intel_3_t <- sampling(ARCH, data = list(N = length(intel), y = intel %>% as.vector(), Kar = 3, family = 1), chains = 1, refresh=0)
print(fit_intel_3_t, pars=names(fit_intel_3_t)[1:5])

## @knitr plot_ARCH_t_residuals
sq <- extract(fit_intel_3_t, pars="sigma")[[1]] %>% colMeans()
mu <- extract(fit_intel_3_t, pars="mu")[[1]] %>% mean()
residuals <- ((intel - mu) / sq)
plot(residuals)
Box.test(residuals, lag=12, type="Ljung-Box")

## @knitr ARCH_t_compare
loo_compare(loo(fit_intel_3), loo(fit_intel_3_t))
