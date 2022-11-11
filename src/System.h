
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
  std::vector<double> lambda;
  double decay_rate;
  double sens;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  System() {};
  
  // public methods
  void load(Rcpp::List args_data, Rcpp::List args_params);
  
};
