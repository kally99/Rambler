
#include "System.h"
#include "misc_v16.h"

using namespace std;

//------------------------------------------------
// load parameters and functions
void System::load(Rcpp::List args_data, Rcpp::List args_params) {
  
  // data
  data_array = rcpp_to_array_bool(args_data["dat_list"]);
  n_ind = data_array.size();
  n_haplo = data_array[0].size();
  samp_time = rcpp_to_vector_double(args_data["samp_time"]);
  n_samp = samp_time.size();
  samp_time_start = samp_time[0];
  samp_time_end = samp_time[n_samp - 1];
  
  // other parameters
  haplo_freqs = rcpp_to_vector_double(args_params["haplo_freqs"]);
  theta = rcpp_to_double(args_params["theta"]);
  decay_rate_meanlog = rcpp_to_double(args_params["decay_rate_meanlog"]);
  decay_rate_sdlog = rcpp_to_double(args_params["decay_rate_sdlog"]);
  sens_shape1 = rcpp_to_double(args_params["sens_shape1"]);
  sens_shape2 = rcpp_to_double(args_params["sens_shape2"]);
  mu_mean = rcpp_to_double(args_params["mu_mean"]);
  mu_sd = rcpp_to_double(args_params["mu_sd"]);
  sigma_shape = rcpp_to_double(args_params["sigma_shape"]);
  sigma_scale = rcpp_to_double(args_params["sigma_scale"]);
  
  // misc
  MH_stepsize = 1.0;
  MH_target_acceptance = 0.44;
  
}
