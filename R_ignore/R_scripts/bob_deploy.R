
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

set.seed(4)

# ------------------------------------------------------------------

# make up some data
n_ind <- 20
n_haplo <- 10
samp_time <- seq(0, 20, 1)
n_samp <- length(samp_time)
haplo_freqs <- rep(3 / n_haplo, n_haplo)
lambda <- rep(0.04, n_ind)
decay_rate <- 0.2
sens <- 0.9

# simulate cohort
dat_list <- sim_cohort(n_ind, samp_time, haplo_freqs, lambda, decay_rate, sens)
df_data <- dat_list$df_data

# get data.frame of true infection times
t_inf_list <- mapply(function(x) x$t_inf, dat_list$raw_list, SIMPLIFY = FALSE)
t_end_list <- mapply(function(x) x$t_end, dat_list$raw_list, SIMPLIFY = FALSE)
t_inf_df <- data.frame(ind = rep(1:n_ind, times = mapply(length, t_inf_list)),
                       t_inf = unlist(t_inf_list),
                       t_end = unlist(t_end_list)) %>%
  mutate(t_end = ifelse(t_end > max(samp_time), max(samp_time), t_end))

# plot data and overlay true infection times
df_data %>%
  group_by(ind, time) %>%
  summarise(COI = sum(positive)) %>%
  mutate(width = c(diff(samp_time), 1)[match(time, samp_time)]) %>%
  ggplot() + theme_bw() +
  geom_tile(aes(x = time + width / 2, y = ind, fill = COI, width = width)) +
  scale_fill_viridis_c() +
  geom_point(aes(x = t_inf, y = ind), col = "white", size = 0.7, data = t_inf_df) +
  geom_segment(aes(x = t_inf, xend = t_end, y = ind, yend = ind), col = "white",
               alpha = 0.5, data = t_inf_df)

# run MCMC
my_mcmc <- run_mcmc(df_data = df_data,
                    haplo_freqs = haplo_freqs,
                    lambda = lambda,
                    decay_rate,
                    sens = sens,
                    burnin = 1e2,
                    samples = 1e3)


i <- 9
z <- t(mapply(function(x) x[[i]], my_mcmc$time_inf_sampling))
plot(z[,1], type = 'l', ylim = range(samp_time))
lines(z[,2], col = 2)
lines(z[,3], col = 3)
abline(h = dat_list$raw_list[[i]]$t_inf, lty = 3)
dat_list$raw_list[[i]]$t_inf


plot(my_mcmc$diagnostics$MC_accept_burnin, ylim = c(0, 1))
plot(my_mcmc$diagnostics$MC_accept_sampling, ylim = c(0, 1))

# have a look at output
head(my_mcmc$draws)
tail(my_mcmc$draws)

