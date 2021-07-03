#' Data reader and common utility functions

## @knitr load_packages
{
  require(tidyverse)
  require(rstan)
  rstan_options(auto_write = TRUE)
  rstan_options(javascript = FALSE)
  theme_set(theme_classic(base_size = 24))
}
