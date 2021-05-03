#' Data reader and common utility functions

## @knitr load_packages
require(tidyverse)
require(rstan)
require(bayesplot)
require(loo)
require(LaplacesDemon)
library(xts)
library(quantmod)
library(FKF)
library(MARSS)
rstan_options(auto_write = TRUE)
rstan_options(javascript = FALSE)
theme_set(theme_classic(base_size = 24))
