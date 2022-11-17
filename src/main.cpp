
#include <chrono>
#include "main.h"
#include "MCMC.h"

using namespace std;

//------------------------------------------------
// run basic example mcmc
Rcpp::List run_mcmc_cpp(Rcpp::List args_data, Rcpp::List args_params,
                        Rcpp::List args_MCMC, Rcpp::List args_progress,
                        Rcpp::List args_functions) {
  
  // start timer
  chrono::high_resolution_clock::time_point t0 = chrono::high_resolution_clock::now();
  
  // extract R functions
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // create MCMC object and load arguments
  MCMC mcmc(args_data, args_params, args_MCMC, args_progress);
  
  // run burn-in and sampling phases of MCMC
  mcmc.run_mcmc_burnin(update_progress);
  mcmc.run_mcmc_sampling(update_progress);
  
  // end timer
  double t_diff = chrono_timer(t0, "chain completed in ", "\n", !mcmc.silent);
  
  // return outputs in list
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("time_inf_burnin") = mcmc.time_inf_burnin,
                                      Rcpp::Named("time_inf_sampling") = mcmc.time_inf_sampling,
                                      Rcpp::Named("lambda_burnin") = mcmc.lambda_burnin,
                                      Rcpp::Named("lambda_sampling") = mcmc.lambda_sampling,
                                      Rcpp::Named("decay_rate_burnin") = mcmc.decay_rate_burnin,
                                      Rcpp::Named("decay_rate_sampling") = mcmc.decay_rate_sampling,
                                      Rcpp::Named("sens_burnin") = mcmc.sens_burnin,
                                      Rcpp::Named("sens_sampling") = mcmc.sens_sampling,
                                      Rcpp::Named("theta_burnin") = mcmc.theta_burnin,
                                      Rcpp::Named("theta_sampling") = mcmc.theta_sampling,
                                      Rcpp::Named("MC_accept_burnin") = mcmc.MC_accept_burnin,
                                      Rcpp::Named("MC_accept_sampling") = mcmc.MC_accept_sampling,
                                      Rcpp::Named("t_diff") = t_diff);
  return ret;
}
