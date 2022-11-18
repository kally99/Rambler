
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

set.seed(2)

# ------------------------------------------------------------------

# make up some data
n_ind <- 10
n_haplo <- 20
samp_time <- seq(0, 20, 2)
n_samp <- length(samp_time)
haplo_freqs <- rep(3 / n_haplo, n_haplo)
#lambda <- rep(0.05, n_ind)
lambda <- rep(0.3, n_ind)
#lambda <- seq(0.01, 0.5, l = n_ind)
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
burnin <- 1e3
samples <- 1e3
#beta <- seq(0, 1, 0.1)
beta <- 1
my_mcmc <- run_mcmc(df_data = df_data,
                    haplo_freqs = haplo_freqs,
                    burnin = burnin,
                    samples = samples,
                    beta = beta,
                    silent = FALSE,
                    sens_shape1 = 9000, sens_shape2 = 1000, mu = -1.75, sigma = 1.2)

#plot(my_mcmc$diagnostics$MC_accept_burnin)

z <- my_mcmc$output %>%
  #filter(ind == 3) %>%
  filter(phase == "sampling") %>%
  filter(!(param %in% c("sensitivity", "decay_rate")))
z2 <- expand_grid(iteration = 1:samples+burnin,
                  ind = 1:n_ind,
                  param = sprintf("inf_time_%s", 1:30))
z3 <- left_join(z2, z) %>%
  group_by(ind, iteration) %>%
  summarise(n = sum(!is.na(value))) %>%
  pull(n)

y <- tabulate(z3 + 1) / sum(tabulate(z3 + 1))
xvec <- seq_along(y) - 1
plot(xvec, y)
lines(xvec, dpois(xvec, lambda = diff(range(samp_time))*0.3))


# trace of decay rate
my_mcmc$output %>%
  filter(param == "decay_rate") %>%
  ggplot() + theme_bw() +
  geom_point(aes(x = iteration, y = value)) +
  geom_hline(yintercept = decay_rate, linetype = "dashed") +
  ggtitle("decay_rate")

# trace of lambda
i <- 8
my_mcmc$output %>%
  filter(param == sprintf("lambda_%s", i)) %>%
  ggplot() + theme_bw() +
  geom_point(aes(x = iteration, y = value)) +
  geom_hline(yintercept = lambda[i], linetype = "dashed") +
  ggtitle("lambda")

# trace of sensitivity
my_mcmc$output %>%
  filter(param == "sensitivity") %>%
  ggplot() + theme_bw() +
  geom_point(aes(x = iteration, y = value)) +
  geom_hline(yintercept = sens, linetype = "dashed") +
  ggtitle("sensitivity")


quantile_95 <- function(x) {
  quantile(x, probs = c(0.025, 0.5, 0.975))
}

# CrI of lambda
my_mcmc$output %>%
  filter(phase == "sampling") %>%
  filter(str_detect(param, "lambda")) %>%
  group_by(ind) %>%
  summarise(ind = rep(ind[1], 3),
            summary = c("Q2.5", "Q50", "Q97.5"),
            value = quantile_95(value)) %>%
  pivot_wider(names_from = summary) %>%
  ggplot() + theme_bw() +
  geom_pointrange(aes(x = ind, ymin = Q2.5, y = Q50, ymax = Q97.5)) +
  geom_point(x = 1:n_ind, y = lambda, color = "red")

# trace plot of infection times
i <- 2
my_mcmc$output %>%
  filter(phase == "sampling") %>%
  filter(ind == i) %>%
  filter(str_detect(param, "inf_time")) %>%
  ggplot() + theme_bw() +
  geom_point(aes(x = iteration, y = value, color = param)) +
  ylim(range(samp_time)) +
  geom_hline(yintercept = dat_list$raw_list[[i]]$t_inf, linetype = "dashed") +
  xlab("Iteration") + ylab("Infection time")

#ggsave("/Users/rverity/Desktop/trace_plot.png")

# get bandwidth over all samples together
bw <- my_mcmc$output %>%
  filter(phase == "sampling") %>%
  filter(str_detect(param, "inf_time")) %>%
  pull(value) %>%
  bw.nrd0()

# get kernel density over all individuals
df_density <- my_mcmc$output %>%
  filter(phase == "sampling") %>%
  filter(str_detect(param, "inf_time")) %>%
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
