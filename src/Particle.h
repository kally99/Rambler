
#pragma once

#include "System.h"
#include "misc_v16.h"
#include "probability_v18.h"

#include <Rcpp.h>

//------------------------------------------------
// class defining MCMC particle
class Particle {
  
public:
  // PUBLIC OBJECTS
  
  // pointer to system object
  System * s_ptr;
  
  // model parameters
  std::vector<std::vector<double>> time_inf;
  
  // proposal objects
  std::vector<double> time_inf_prop;
  
  // likelihood and prior
  std::vector<double> loglike_ind;
  std::vector<double> logprior_ind;
  double loglike;
  double logprior;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Particle() {};
  
  // initialisers
  void init(System &s);
  
  // likelihood and prior
  double get_loglike_ind(int ind, const std::vector<double> &time_inf);
  double get_logprior_ind(int ind, const std::vector<double> &time_inf);
  
  // update functions
  void update(double beta);
  void MH_time_inf(double beta);
  void split_merge_time_inf(double beta);
  
};
