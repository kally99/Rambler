
#pragma once

#include <Rcpp.h>
#include <vector>

//------------------------------------------------
// class holding all common objects that every particle needs access to
class System {
  
public:
  // PUBLIC OBJECTS
  
  // data
  std::vector<std::vector<std::vector<bool>>> data_array;
  int n_ind;
  int n_haplo;
  std::vector<double> samp_time;
  int n_samp;
  double samp_time_start;
  double samp_time_end;
  
  // other parameters
  std::vector<double> haplo_freqs;
  double decay_rate_meanlog;
  double decay_rate_sdlog;
  double sens_shape1;
  double sens_shape2;
  double mu;
  double sigma;
  
  // misc
  double MH_stepsize;
  double MH_target_acceptance;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  System() {};
  
  // public methods
  void load(Rcpp::List args_data, Rcpp::List args_params);
  
};
