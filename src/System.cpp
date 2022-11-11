
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
  lambda = rcpp_to_vector_double(args_params["lambda"]);
  decay_rate = args_params["decay_rate"];
  sens = args_params["sens"];
  
}
