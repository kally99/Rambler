
#------------------------------------------------
#' @title Check that Rambler package has loaded successfully
#'
#' @description Simple function to check that Rambler package has loaded 
#'   successfully. Prints "Rambler loaded successfully!" if so.
#'
#' @export

check_rambler_loaded <- function() {
  message("Rambler loaded successfully!")
}

#------------------------------------------------
# check data format, and restruture into format expected by C++
#' @noRd
restructure_data <- function(df_data) {
  
  # basic checks
  assert_dataframe(df_data)
  assert_in(c("ind", "haplo", "time", "positive"), names(df_data))
  assert_pos_int(df_data$haplo)
  assert_pos(df_data$time)
  assert_in(df_data$positive, c(0, 1))
  
  # split into nested list over individuals then haplotypes. Only store observed
  # positive status
  dat_list <- mapply(function(x) {
    mapply(function(y) {
      y$positive
    }, split(x, f = x$haplo), SIMPLIFY = FALSE)
  }, split(df_data, df_data$ind), SIMPLIFY = FALSE)
  
  # return list
  ret <- list(dat_list = dat_list,
              samp_time = unique(df_data$time))
  return(ret)
}

#------------------------------------------------
#' @title Test function to run an example MCMC
#'
#' @description This should be replaced by a more carefully thought out
#'   structure, probably involving defining a parameters data.frame outside of
#'   this function that can be loaded in.
#'
#' @param df_data data.frame of haplotypes in each individual at each time
#'   point, in long format.
#' @param haplo_freqs vector giving the probability that each haplotype is
#'   transmitted in a given infectious bite (not quite the same thing as
#'   haplotype frequencies in the population).
#' @param decay_rate_meanlog,decay_rate_sdlog mean and standard deviation (on
#'   the log scale) of the prior distribution on the decay rate.
#' @param sens_shape1,sens_shape2 shape parameters of beta prior on sensitivity.
#' @param mu_mean,mu_sd mean and standard deviation of the normal hyper-prior on
#'   mu, which is the log-mean of the log-normal prior on lambda.
#' @param sigma_shape,sigma_scale shape and scale parameters of the inverse
#'   gamma hyper-prior on sigma. The mean of this distribution is given by:
#'   \deqn{\beta / (\alpha - 1)} for \eqn{\alpha > 1}, and the variance is given
#'   by \deqn{\beta^2 / ((\alpha - 1)^2(\alpha - 2))} where \eqn{\alpha} is the
#'   shape and \eqn{\beta} is the scale.
#' @param burnin the number of burn-in iterations.
#' @param samples the number of sampling iterations.
#' @param beta vector of thermodynamic powers. Final value in the vector should
#'   always be 1.
#' @param pb_markdown whether to run progress bars in markdown mode, meaning
#'   they are only updated when they reach 100\%. This avoids large amounts of
#'   output being printed to markdown files.
#' @param silent whether to suppress all console output.
#'
#' @importFrom stats setNames
#' @importFrom utils txtProgressBar
#' @export

run_mcmc <- function(df_data,
                     haplo_freqs,
                     decay_rate_meanlog = 0.0,
                     decay_rate_sdlog = 1.0,
                     sens_shape1 = 10,
                     sens_shape2 = 1,
                     mu_mean = 0.0,
                     mu_sd = 10.0,
                     sigma_shape = 1.0,
                     sigma_scale = 1.0,
                     burnin = 1e2,
                     samples = 1e3,
                     beta = 1,
                     pb_markdown = FALSE,
                     silent = FALSE) {
  
  # check data format, and restruture into format expected by C++
  data_processed <- restructure_data(df_data)
  n_ind <- length(data_processed$dat_list)
  n_haplo <- length(data_processed$dat_list[[1]])
  n_samp <- length(data_processed$samp_time)
  
  # check inputs
  assert_vector_bounded(haplo_freqs)
  assert_single_numeric(decay_rate_meanlog)
  assert_single_pos(decay_rate_sdlog)
  assert_single_pos(sens_shape1)
  assert_single_pos(sens_shape2)
  assert_single_numeric(mu_mean)
  assert_single_pos(mu_sd)
  assert_single_pos(sigma_shape)
  assert_single_pos(sigma_scale)
  assert_single_pos_int(burnin, zero_allowed = FALSE)
  assert_single_pos_int(samples, zero_allowed = FALSE)
  assert_vector_bounded(beta)
  assert_eq(beta[length(beta)], 1)
  assert_single_logical(pb_markdown)
  assert_single_logical(silent)
  
  # make a list of model parameters
  args_params <- list(haplo_freqs = haplo_freqs,
                      decay_rate_meanlog = decay_rate_meanlog,
                      decay_rate_sdlog = decay_rate_sdlog,
                      sens_shape1 = sens_shape1,
                      sens_shape2 = sens_shape2,
                      mu_mean = mu_mean,
                      mu_sd = mu_sd,
                      sigma_shape = sigma_shape,
                      sigma_scale = sigma_scale)
  
  # make a list of MCMC parameters
  args_MCMC <- list(burnin = burnin,
                    samples = samples,
                    beta = beta,
                    pb_markdown = pb_markdown,
                    silent = silent)
  
  # make progress bars
  pb_burnin <- txtProgressBar(min = 0, max = burnin, initial = NA, style = 3)
  pb_samples <- txtProgressBar(min = 0, max = samples, initial = NA, style = 3)
  args_progress <- list(pb_burnin = pb_burnin,
                        pb_samples = pb_samples)
  
  # R functions to pass to C++
  args_functions <- list(update_progress = update_progress)
  
  # ---------- run MCMC ----------
  
  # run efficient C++ function
  output_raw <- run_mcmc_cpp(data_processed, args_params, args_MCMC,
                             args_progress, args_functions)
  
  #return(output_raw)
  
  # ---------- process output ----------
  
  # get data.frame of infection times for all individuals. Note that the number
  # of infection events can change from one iteration to the next, hence this is
  # in long form
  time_inf_list <- c(output_raw$time_inf_burnin,
                     output_raw$time_inf_sampling)
  df_time_inf <- mapply(function(i) {
    x <- time_inf_list[[i]]
    lx <- mapply(length, x)
    data.frame(phase = ifelse(i <= burnin, "burnin", "sampling"),
               iteration = i,
               ind = rep(1:n_ind, times = lx),
               param = sprintf("inf_time_%s", unlist(mapply(seq_len, lx))),
               value = unlist(x))
  }, seq_along(time_inf_list), SIMPLIFY = FALSE) %>%
    bind_rows()
  
  # get data.frame of lambda
  lambda_list <- c(output_raw$lambda_burnin,
                   output_raw$lambda_sampling)
  df_lambda <- mapply(function(i) {
    x <- lambda_list[[i]]
    data.frame(phase = ifelse(i <= burnin, "burnin", "sampling"),
               iteration = i,
               ind = 1:n_ind,
               param = sprintf("lambda_%s", 1:n_ind),
               value = unlist(x))
  }, seq_along(lambda_list), SIMPLIFY = FALSE) %>%
    bind_rows()
  
  # get data.frame of decay rate
  df_decay_rate <- rbind(data.frame(phase = "burnin",
                                    iteration = 1:burnin,
                                    ind = NA,
                                    param = "decay_rate",
                                    value = output_raw$decay_rate_burnin),
                         data.frame(phase = "sampling",
                                    iteration = 1:samples + burnin,
                                    ind = NA,
                                    param = "decay_rate",
                                    value = output_raw$decay_rate_sampling))
  
  # get data.frame of sensitivity
  df_sens <- rbind(data.frame(phase = "burnin",
                              iteration = 1:burnin,
                              ind = NA,
                              param = "sensitivity",
                              value = output_raw$sens_burnin),
                   data.frame(phase = "sampling",
                              iteration = 1:samples + burnin,
                              ind = NA,
                              param = "sensitivity",
                              value = output_raw$sens_sampling))
  
  # get MCMC diagnostics
  diagnostics = list(MC_accept_burnin = output_raw$MC_accept_burnin / burnin,
                     MC_accept_sampling = output_raw$MC_accept_sampling / samples)
  
  
  # return list
  ret <- list(output = rbind(df_time_inf,
                             df_lambda,
                             df_decay_rate,
                             df_sens),
              diagnostics = diagnostics)
  return(ret)
}

#------------------------------------------------
# update progress bar
# pb_list = list of progress bar objects
# name = name of this progress bar
# i = new value of bar
# max_i = max value of bar (close when reach this value)
# close = whether to close when reach end
#' @importFrom utils setTxtProgressBar
#' @noRd
update_progress <- function(pb_list, name, i, max_i, close = TRUE) {
  setTxtProgressBar(pb_list[[name]], i)
  if (i == max_i & close) {
    close(pb_list[[name]])
  }
}
