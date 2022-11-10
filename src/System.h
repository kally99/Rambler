
#pragma once

#include <Rcpp.h>
#include <vector>

//------------------------------------------------
// class holding all common objects that every particle needs access to
class System {
  
public:
  // PUBLIC OBJECTS
  
  // data
  //std::vector<int> a;
  //std::vector<int> r;
  //int n_loci;
  
  // other parameters
  /*
  std::vector<double> p;
  double c;
  */
  
  // PUBLIC FUNCTIONS
  
  // constructors
  System() {};
  
  // public methods
  void load(Rcpp::List args_data, Rcpp::List args_params);
  
};
