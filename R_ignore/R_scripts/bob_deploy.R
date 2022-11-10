
# bob_deploy.R
#
# Author: Bob Verity
# Date: 2022-11-10
#
# Purpose:
# Test package functions.

# RStudio shortcuts:
#    cmd+shift+L     : load package from local version
#    cmd+shift+D     : document (NB, option must be activated in Build Tools)
#    cmd+shift+E     : check
#    cmd+shift+T     : test

# Useful commands:
# devtools::install()  # install package
# pkgdown::build_site() # build all pages of pkgdown website
# pkgdown::build_article('blocks')  # build single vignette
# check('.', args = '--no-examples')  # run checks without examples
# covr::report()    # interactive coverage report
# devtools::build_vignettes()
# ------------------------------------------------------------------

#library(Rcpp)
library(dplyr)
library(ggplot2)

set.seed(1)

# ------------------------------------------------------------------

# make up some data
n_ind <- 5
n_haplo <- 10
samp_time <- seq(0, 20, 2)
n_samp <- length(samp_time)
haplo_freqs <- rep(2 / n_haplo, n_haplo)
lambda <- rep(0.1, n_ind)
decay_rate <- 0.1
sens <- 0.9

# simulate cohort
dat_list <- sim_cohort(4, samp_time, haplo_freqs, lambda, decay_rate, sens)
df_data <- dat_list$df_data


# run MCMC
my_mcmc <- run_mcmc(df_data = df_data,
                    haplo_freqs = haplo_freqs,
                    lambda = lambda,
                    burnin = 1e2,
                    samples = 1e3,
                    beta = beta)







plot(my_mcmc$diagnostics$MC_accept_burnin, ylim = c(0, 1))
plot(my_mcmc$diagnostics$MC_accept_sampling, ylim = c(0, 1))

# have a look at output
head(my_mcmc$draws)
tail(my_mcmc$draws)

# plot posterior mu draws
my_mcmc$draws %>%
  ggplot() + theme_bw() +
  geom_point(aes(x = iteration, y = mu_1), size = 0.5) +
  geom_point(aes(x = iteration, y = mu_2), size = 0.5, col = "red") +
  ylim(c(0, 1)) +
  ylab("mu")

# plot posterior sigma draws
my_mcmc$draws %>%
  ggplot() + theme_bw() +
  geom_point(aes(x = iteration, y = sigma))

# plot posterior w draws
my_mcmc$draws %>%
  ggplot() + theme_bw() +
  geom_point(aes(x = iteration, y = w)) +
  ylim(c(0, 1))
