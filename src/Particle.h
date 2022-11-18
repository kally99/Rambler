
#pragma once

#include "System.h"
#include "misc_v16.h"
#include "probability_v19.h"

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
  std::vector<double> lambda;
  double decay_rate;
  double sens;
  
  // proposal objects
  std::vector<double> time_inf_prop;
  std::vector<double> loglike_ind_prop;
  
  // likelihood and prior
  std::vector<double> loglike_ind;
  double loglike;
  double logprior;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Particle() {};
  
  // initialisers
  void init(System &s);
  
  // likelihood and prior
  double get_loglike_ind(int ind, double lambda, double decay_rate,
                         double sens, const std::vector<double> &time_inf);
  double get_logprior_time_inf(const std::vector<double> &time_inf, double lambda);
  double get_logprior_lambda(double lambda);
  double get_logprior_decay_rate(double decay_rate);
  double get_logprior_sens(double sens);
  
  // update functions
  void update(double beta, double &time_inf_bw, std::vector<double> &lambda_bw,
              double &decay_rate_bw, double &sens_bw,
              int rep, bool Robbins_Monro);
  void MH_time_inf(double beta, double &time_inf_bw, int rep, bool Robbins_Monro);
  void split_merge_time_inf(double beta);
  void MH_lambda(double beta, std::vector<double> &lambda_bw, int rep, bool Robbins_Monro);
  void MH_decay_rate(double beta, double &decay_rate_bw, int rep, bool Robbins_Monro);
  void MH_sens(double beta, double &sens_bw, int rep, bool Robbins_Monro);
  
};
