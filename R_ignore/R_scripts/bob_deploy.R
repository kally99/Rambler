
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

set.seed(3)

# ------------------------------------------------------------------

# make up some data
n_ind <- 20
n_haplo <- 10
samp_time <- seq(0, 40, 2)
n_samp <- length(samp_time)
haplo_freqs <- rep(3 / n_haplo, n_haplo)
lambda <- rep(0.05, n_ind)
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
               alpha = 0.5, data = t_inf_df) +
  xlab("Time") + ylab("Individual")

#ggsave("/Users/rverity/Desktop/COI_plot.png")

# run MCMC
burnin <- 1e2
samples <- 1e3
my_mcmc <- run_mcmc(df_data = df_data,
                    haplo_freqs = haplo_freqs,
                    lambda = lambda,
                    decay_rate,
                    sens = sens,
                    burnin = burnin,
                    samples = samples)


# trace plot
i <- 14
my_mcmc$output %>%
  filter(phase == "sampling" & ind == i) %>%
  ggplot() + theme_bw() +
  geom_point(aes(x = iteration, y = value, color = param)) +
  ylim(range(samp_time)) +
  geom_hline(yintercept = dat_list$raw_list[[i]]$t_inf, linetype = "dashed") +
  xlab("Iteration") + ylab("Infection time")

#ggsave("/Users/rverity/Desktop/trace_plot.png")

# get bandwidth over all samples together
bw <- my_mcmc$output %>%
  filter(phase == "sampling") %>%
  pull(value) %>%
  bw.nrd0()

# get kernel density over all individuals
df_density <- my_mcmc$output %>%
  filter(phase == "sampling") %>%
  group_by(param, ind) %>%
  summarise(param = param[1],
            ind = ind[1],
            x = seq(first(samp_time), last(samp_time), l = 201),
            y = length(value) / samples * density(value, from = first(samp_time), to = last(samp_time), n = 201, bw = 0.5)$y)

# plot
df_density %>%
  ggplot() + theme_bw() +
  geom_area(aes(x = x, y = y, fill = param)) +
  facet_wrap(~ind) +
  geom_segment(aes(x = t_inf, xend = t_inf, y = 0.2, yend = 0),
               arrow = arrow(length = unit(0.1, "npc")), data = t_inf_df) +
  xlab("Time") + ylab("Posterior probability")

#ggsave("/Users/rverity/Desktop/density_plot.png")


#plot(my_mcmc$diagnostics$MC_accept_burnin, ylim = c(0, 1))
#plot(my_mcmc$diagnostics$MC_accept_sampling, ylim = c(0, 1))
