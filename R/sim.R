
#------------------------------------------------
#' @title Simulate a single individual
#'
#' @description Simulate a single individual from the generative model. Returns outpus as
# list, including observed data as well as hidden true states. Note that the
# method for simulating data is deliberately different from the way this is
# represented in inference - this is deliberate as it will allow us to check
# that our different representations of the process are consistent
#'
#' @export

sim_ind <- function(samp_time, haplo_freqs, lambda, decay_rate, sens) {
  
  # get dimensions
  n_haplo <- length(haplo_freqs)
  n_samp <- length(samp_time)
  
  # define a matrix for holding the true infection state of every haplotype
  # (rows) at every timepoint (columns)
  state_true <- matrix(0, n_haplo, n_samp)
  
  # consider whether this individual initialises positive for each haplotype. The
  # probability of initialising positive is assumed to be given by the relative
  # rates of acquiring haplotypes vs. losing them.
  for (j in 1:n_haplo) {
    p_init_infected <- lambda*haplo_freqs[j] / (lambda*haplo_freqs[j] + decay_rate)
    init_infected <- rbinom(1, 1, p_init_infected)
    if (init_infected) {
      t_decay <- rexp(1, rate = decay_rate)
      state_true[j, (samp_time < t_decay)] <- 1
    }
  }
  
  # draw the number of new infections that occur during the observation period and
  # the timings of these infections
  samp_period <- diff(range(samp_time))
  n_inf <- rpois(1, lambda*samp_period)
  t_inf <- sort(runif(n_inf, 0, samp_period))
  
  # loop through all new infections (if any present)
  for (k in seq_len(n_inf)) {
    
    # draw which haplotypes are introduced. For each haplotype, draw time at which
    # it decays and update state matrix accordingly
    which_haplos <- which(rbinom(n_haplo, 1, haplo_freqs) == 1)
    for (j in seq_along(which_haplos)) {
      t_decay <- t_inf[k] + rexp(1, rate = decay_rate)
      state_true[j, (samp_time > t_inf[k]) & (samp_time < t_decay)] <- 1
    }
  }
  
  # draw the observed state from the true state by taking into account sensitivity
  state_obs <- state_true * matrix(rbinom(n_haplo * n_samp, 1, sens), n_haplo)
  
  # return as list
  ret <- list(state_true = state_true,
              state_obs = state_obs,
              t_inf = t_inf)
  return(ret)
}



#------------------------------------------------
#' @title Simulate a complete cohort
#'
#' @description Simulate a cohort of individuals from the generative model. Runs
#'   \code{sim_ind()} multiple times, and reformats observed output into
#'   data.frame. Also keeps hold of the true states of all individuals for
#'   reference.
#'
#' @export

sim_cohort <- function(n, samp_time, haplo_freqs, lambda, decay_rate, sens) {
  
  # simulate each individual
  raw_list <- list()
  for (i in 1:n) {
    raw_list[[i]] <- sim_ind(samp_time, haplo_freqs, lambda[i], decay_rate, sens)
  }
  
  # get observed data into data.frame
  df_data <- mapply(function(i) {
    data.frame(ind = i,
               haplo = seq_along(haplo_freqs),
               time = rep(samp_time, each = length(haplo_freqs)),
               positive = as.vector(raw_list[[i]]$state_obs))
  }, seq_along(raw_list), SIMPLIFY = FALSE) %>%
    bind_rows()
  
  # return as list
  ret <- list(df_data = df_data,
              raw_list = raw_list)
  return(ret)
}
