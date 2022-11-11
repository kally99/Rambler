
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
                                      Rcpp::Named("t_diff") = t_diff);
  return ret;
}
