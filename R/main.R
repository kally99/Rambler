
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
#' @param lambda vector of FOI in each individual.
#' @param decay_rate rate at which each haplotype clears.
#' @param sens sensitivity of sequencing (assumed the same for all haplotypes).
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
                     lambda,
                     decay_rate,
                     sens,
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
  assert_vector_pos(lambda)
  assert_single_pos(decay_rate)
  assert_single_bounded(sens)
  assert_single_pos_int(burnin, zero_allowed = FALSE)
  assert_single_pos_int(samples, zero_allowed = FALSE)
  assert_vector_bounded(beta)
  assert_eq(beta[length(beta)], 1)
  assert_single_logical(pb_markdown)
  assert_single_logical(silent)
  
  # make a list of model parameters
  args_params <- list(haplo_freqs = haplo_freqs,
                      lambda = lambda,
                      decay_rate = decay_rate,
                      sens = sens)
  
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
  
  # get dataframe of infection times for all individuals. Note that the number
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
  
  # return
  ret <- list(output = df_time_inf,
              diagnostics = list(MC_accept_burnin = output_raw$MC_accept_burnin / burnin,
                                 MC_accept_sampling = output_raw$MC_accept_sampling / samples))
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
