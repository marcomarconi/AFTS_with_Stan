#' Data reader and common utility functions

## @knitr load_packages
require(tidyverse)
require(rstan)
require(ggplot2)
require(bayesplot)
require(loo)
require(LaplacesDemon)
rstan_options(auto_write = TRUE)
rstan_options(javascript = FALSE)
theme_set(theme_classic(base_size = 24))
