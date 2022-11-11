
#include "Particle.h"

using namespace std;

//------------------------------------------------
// initialise/reset particle
void Particle::init(System &s) {
  
  // pointer to system object
  s_ptr = &s;
  
  // initialise model parameters
  time_inf = vector<vector<double>>(s_ptr->n_ind, vector<double>(1));
  loglike_ind = vector<double>(s_ptr->n_ind);
  logprior_ind = vector<double>(s_ptr->n_ind);
  for (int i = 0; i < s_ptr->n_ind; ++i) {
    
    // draw infection times uniformly over sampling period and place in
    // ascending order
    for (int j = 0; j < time_inf[i].size(); ++j) {
      time_inf[i][j] = runif1(s_ptr->samp_time_start, s_ptr->samp_time_end);
    }
    sort(time_inf[i].begin(), time_inf[i].end());
    
    // calculate loglikelihood and logprior for this individual
    loglike_ind[i] = get_loglike_ind(i, time_inf[i]);
    logprior_ind[i] = get_logprior_ind(i, time_inf[i]);
  }
  loglike = sum(loglike_ind);
  logprior = sum(logprior_ind);
  
  // initialise proposal objects
  time_inf_prop = vector<double>(time_inf[0].size());
}

//------------------------------------------------
// calculate loglikelihood for a single individual
double Particle::get_loglike_ind(int ind, const vector<double> &time_inf) {
  
  // the final likelihood for this individual will be the product over all
  // haplotypes, and all time series of observations for this haplotype. This
  // will be the product of a large number of terms, and so has a high chance of
  // underflow if not done correctly. At the same time, it would not be
  // efficient to work in log space throughout the calculation, as we will be
  // summing terms inside the likelihood as well as multiplying them. For this
  // reason, stick with a middle ground in which terms are multiplied, but at
  // each step we check if we are getting dangerously small (i.e. <1e-200) and
  // if so we take out a scaling factor and renormalise. loglike_running will
  // keep track of the total scaling factor (logged) throughout, and eventually
  // will come to hold the final loglikelihood
  double loglike_running = 0.0;
  
  // loop over haplotypes
  for (int j = 0; j < s_ptr->n_haplo; ++j) {
    
    // calculate the probability of initialising in the positive state, based on
    // relative rates of infection and clearance
    double prob_init_positive = s_ptr->lambda[ind]*s_ptr->haplo_freqs[j] / (s_ptr->lambda[ind]*s_ptr->haplo_freqs[j] + s_ptr->decay_rate);
    
    // forward_positive will store the joint probability of 1) being in the
    // positive state at the current time, and 2) the probability of all data up
    // to and including the current time. forward_negative does the same for the
    // negative state. These will be updated as we walk through the HMM, but for
    // now should hold the probabilities at initialisation
    double forward_positive = (s_ptr->data_array[ind][j][0] == 1) ? prob_init_positive * s_ptr->sens : prob_init_positive * (1 - s_ptr->sens);
    double forward_negative = (s_ptr->data_array[ind][j][0] == 1) ? 0 : (1 - prob_init_positive);
    
    // as we walk through the HMM we will focus on sequential time intervals.
    // The focal interval will always run from t0 to t1
    double t0 = s_ptr->samp_time[0];
    double t1 = t0;
    
    // we want to loop through the time intervals in our data, but also we need
    // to consider known infection times. If there is an infection within an
    // observation period then our new focal interval should run up to this
    // infection time. Hence, we need a second loop within our observation loop
    // that runs through infection times. inf_index gives the index of this
    // second loop. Note that this is not a simple nested loop, as the second
    // loop progresses synchronously with the first, and is updated as needed
    int n_inf = time_inf.size();
    int inf_index = 0;
    
    // loop through all sampling times. Start at observation k=1 (the second
    // observation) meaning we are at the right hand side of the focal interval
    for (int k = 1; k < s_ptr->n_samp; ++k) {
      
      // deal with any possible infections within this observation interval. Use
      // a while loop as we don't know how many infections there may be within
      // this interval.
      if (inf_index < n_inf) {
        while (time_inf[inf_index] < s_ptr->samp_time[k]) {
          
          // calculate all transition probabilities to the positive (1) and
          // negative (0) states, acccounting for the possibility that the
          // haplotype is introduced within this infection
          t1 = time_inf[inf_index];
          double trans_11 = s_ptr->haplo_freqs[j] + (1 - s_ptr->haplo_freqs[j])*exp(-s_ptr->decay_rate*(t1 - t0));
          double trans_10 = 1 - trans_11;
          double trans_01 = s_ptr->haplo_freqs[j];
          double trans_00 = 1 - trans_01;
          
          // calculate new forward positive and negative probabilities. Note
          // that there is no data within this calculation (i.e. no emmission
          // probability) because this is an intermediate time between
          // observations.
          // Simply updating forward probabilities could lead to underflow if
          // the proposed transitions are very unlikely. Therefore, calculate
          // proposed values, and only make the update if new values will sum to
          // a reasonable value. Otherwise, deal with underflow
          double forward_positive_prop = forward_positive*trans_11 + forward_negative*trans_01;
          double forward_negative_prop = forward_positive*trans_10 + forward_negative*trans_00;
          if ((forward_positive_prop + forward_negative_prop) < 1e-200) {
            
            // deal with underflow by extracting a scaling factor into the
            // running loglikelihood, and renormalising
            double forward_sum = forward_positive + forward_negative;
            loglike_running += log(forward_sum);
            forward_positive /= forward_sum;
            forward_negative /= forward_sum;
            forward_positive_prop = forward_positive*trans_11 + forward_negative*trans_01;
            forward_negative_prop = forward_positive*trans_10 + forward_negative*trans_00;
          }
          forward_positive = forward_positive_prop;
          forward_negative = forward_negative_prop;
          
          // move time forward and the infection index forward one step
          t0 = t1;
          inf_index++;
          if (inf_index == n_inf) {
            break;  // break while loop if we reach the end of infections
          }
        }
      }
      
      // at this point we have dealt with any possible infections within the
      // interval, so we need to move forward to the next observation. Calculate
      // transition probabilities from the positive state (note that if in
      // negative state we can only remain there)
      t1 = s_ptr->samp_time[k];
      double trans_11 = exp(-s_ptr->decay_rate*(t1 - t0));
      double trans_10 = 1 - trans_11;
      
      // update forward positive and negative probabilities, taking account the
      // observed data. As above, do this in a way that avoids underflow
      double forward_positive_prop = (s_ptr->data_array[ind][j][k] == 1) ? forward_positive*trans_11*s_ptr->sens : forward_positive*trans_11*(1 - s_ptr->sens);
      double forward_negative_prop = (s_ptr->data_array[ind][j][k] == 1) ? 0 : forward_positive*trans_10 + forward_negative;
      if ((forward_positive_prop + forward_negative_prop) < 1e-200) {
        double forward_sum = forward_positive + forward_negative;
        loglike_running += log(forward_sum);
        forward_positive /= forward_sum;
        forward_negative /= forward_sum;
        forward_positive_prop = (s_ptr->data_array[ind][j][k] == 1) ? forward_positive*trans_11*s_ptr->sens : forward_positive*trans_11*(1 - s_ptr->sens);
        forward_negative_prop = (s_ptr->data_array[ind][j][k] == 1) ? 0 : forward_positive*trans_10 + forward_negative;
      }
      forward_positive = forward_positive_prop;
      forward_negative = forward_negative_prop;
      
      // step forward time
      t0 = t1;
    } // end loop over time steps
    
    // sum over the final stage of the forward algorithm to get the final
    // likelihood
    loglike_running += log(forward_positive + forward_negative);
    
  } // end loop over haplotypes
  
  return loglike_running;
}

//------------------------------------------------
// calculate logprior for a single individual
double Particle::get_logprior_ind(int ind, const vector<double> &time_inf) {
  
  // exponential waiting time between infections
  double ret = 0.0;
  double t0 = s_ptr->samp_time_start;
  for (int i = 0; i < time_inf.size(); ++i) {
    ret += log(s_ptr->lambda[ind]) - s_ptr->lambda[ind]*(time_inf[i] - t0);
    t0 = time_inf[i];
  }
  
  // chance of seeing no more infections until end of sampling period
  ret += -s_ptr->lambda[ind]*(s_ptr->samp_time_end - t0);
  
  return ret;
}

//------------------------------------------------
// propose new parameter values and accept/reject
void Particle::update(double beta) {
  
  // distinct update steps for each free parameter
  MH_time_inf(beta);
  split_merge_time_inf(beta);
  
}

//------------------------------------------------
// Metropolis-Hastings on infection times for all individuals
void Particle::MH_time_inf(double beta) {
  
  // update all individuals
  for (int i = 0; i < s_ptr->n_ind; ++i) {
    
    // skip if no infections
    int n_inf = time_inf[i].size();
    if (n_inf == 0) {
      continue;
    }
    
    // copy over infection times
    time_inf_prop = time_inf[i];
    
    // loop through infections
    for (int j = 0; j < n_inf; ++j) {
      
      // get start and end times of proposal interval. This runs from the
      // previous to the next infection, but is truncated at the start and end
      // times of sampling
      double t0 = (j == 0) ? s_ptr->samp_time_start : time_inf_prop[j - 1];
      double t1 = (j == (n_inf - 1)) ? s_ptr->samp_time_end : time_inf_prop[j + 1];
      
      // propose new infection time by drawing from reflected normal within
      // interval
      time_inf_prop[j] = rnorm1_interval(time_inf_prop[j], 1.0, t0, t1);
      
      double loglike_prop = get_loglike_ind(i, time_inf_prop);
      double logprior_prop = get_logprior_ind(i, time_inf_prop);
      
      // calculate Metropolis-Hastings ratio
      double MH = beta*(loglike_prop - loglike_ind[i]) + (logprior_prop - logprior_ind[i]);
      
      double acceptance_linear = 1.0;
      bool accept_move = true;
      if (MH < 0) {
        acceptance_linear = exp(MH);
        accept_move = (runif_0_1() < acceptance_linear);
      }
      
      // accept or reject move
      if (accept_move) {
        time_inf[i][j] = time_inf_prop[j];
        loglike += loglike_prop - loglike_ind[i];
        logprior += logprior_prop - logprior_ind[i];
        loglike_ind[i] = loglike_prop;
        logprior_ind[i] = logprior_prop;
      } else {
        time_inf_prop[j] = time_inf[i][j];
      }
    }
    
  } // end loop through individuals
  
}

//------------------------------------------------
// propose new infection times by dropping an existing infection OR splitting an
// existing interval to add a new infection
void Particle::split_merge_time_inf(double beta) {
  
  // update all individuals
  for (int i = 0; i < s_ptr->n_ind; ++i) {
    
    // copy over infection times
    time_inf_prop = time_inf[i];
    
    // get number of infections
    int n_inf = time_inf_prop.size();
    
    // draw whether split or merge step
    if (rbernoulli1(0.5)) { // split
      
      // choose an interval at random with equal probability. Note that
      // interval_index = 0 denotes the interval from start of sampling to the
      // first infection. Therefore, the start time for subsequent intervals is
      // time_inf_prop[interval_index - 1], i.e. we have to account for the fact
      // we are effectively 1-indexed. Calculate t0 and t1 as the start and end
      // times of the chosen interval
      int interval_index = sample2(0, n_inf);
      double t0 = (interval_index == 0) ? s_ptr->samp_time_start : time_inf_prop[interval_index - 1];
      double t1 = (interval_index == n_inf) ? s_ptr->samp_time_end : time_inf_prop[interval_index];
      
      // propose a new infection time and slot this in at the correct place in
      // the vector
      double prop_time = runif1(t0, t1);
      if (n_inf == 0) {
        time_inf_prop.push_back(prop_time);
      } else {
        time_inf_prop.insert(time_inf_prop.begin() + interval_index, prop_time);
      }
      
      // calculate new likelihood and prior
      double loglike_prop = get_loglike_ind(i, time_inf_prop);
      double logprior_prop = get_logprior_ind(i, time_inf_prop);
      
      // calculate forward and backward proposal probabilities
      double logprop_forward = -log((n_inf + 1) * (t1 - t0));
      double logprop_backward = -log(n_inf + 1);
      
      // calculate Metropolis-Hastings ratio
      double MH = beta*(loglike_prop - loglike_ind[i]) + (logprior_prop - logprior_ind[i]) +
        (logprop_backward - logprop_forward);
      
      double acceptance_linear = 1.0;
      bool accept_move = true;
      if (MH < 0) {
        acceptance_linear = exp(MH);
        accept_move = (runif_0_1() < acceptance_linear);
      }
      
      // accept or reject move
      if (accept_move) {
        time_inf[i] = time_inf_prop;
        loglike += loglike_prop - loglike_ind[i];
        logprior += logprior_prop - logprior_ind[i];
        loglike_ind[i] = loglike_prop;
        logprior_ind[i] = logprior_prop;
      }
      
    } else {  // merge
      
      // nothing to merge if no infections
      if (n_inf == 0) {
        continue;
      }
      
      // choose an infection time at random with equal probability and drop this
      // time. Note that interval_index now represents the index of the
      // infection we are dropping, and so is 0-indexed
      int interval_index = sample2(1, n_inf) - 1;
      double t0 = (interval_index == 0) ? s_ptr->samp_time_start : time_inf_prop[interval_index - 1];
      double t1 = (interval_index == (n_inf - 1)) ? s_ptr->samp_time_end : time_inf_prop[interval_index + 1];
      
      // drop infection and calculate new likelihood and prior
      time_inf_prop.erase(time_inf_prop.begin() + interval_index);
      double loglike_prop = get_loglike_ind(i, time_inf_prop);
      double logprior_prop = get_logprior_ind(i, time_inf_prop);
      
      // calculate forward and backward proposal probabilities
      double logprop_forward = -log(n_inf);
      double logprop_backward = -log(n_inf * (t1 - t0));
      
      // calculate Metropolis-Hastings ratio
      double MH = beta*(loglike_prop - loglike_ind[i]) + (logprior_prop - logprior_ind[i]) +
        (logprop_backward - logprop_forward);
      
      double acceptance_linear = 1.0;
      bool accept_move = true;
      if (MH < 0) {
        acceptance_linear = exp(MH);
        accept_move = (runif_0_1() < acceptance_linear);
      }
      
      // accept or reject move
      if (accept_move) {
        time_inf[i] = time_inf_prop;
        loglike += loglike_prop - loglike_ind[i];
        logprior += logprior_prop - logprior_ind[i];
        loglike_ind[i] = loglike_prop;
        logprior_ind[i] = logprior_prop;
      }
      
    }
  } // end loop over individuals
    
    
    /*
    // copy over infection times
    time_inf_prop = time_inf[i];
    
    // loop through infections
    for (int j = 0; j < n_inf; ++j) {
      
      // get start and end times of proposal interval. This runs from the
      // previous to the next infection, but is truncated at the start and end
      // times of sampling
      double t0 = (j == 0) ? s_ptr->samp_time_start : time_inf_prop[j - 1];
      double t1 = (j == (n_inf - 1)) ? s_ptr->samp_time_end : time_inf_prop[j + 1];
      
      // propose new infection time by drawing from reflected normal within
      // interval
      time_inf_prop[j] = rnorm1_interval(time_inf_prop[j], 1.0, t0, t1);
      
      double loglike_prop = get_loglike_ind(i, time_inf_prop);
      double logprior_prop = get_logprior_ind(i, time_inf_prop);
      
      // calculate Metropolis-Hastings ratio
      double MH = beta*(loglike_prop - loglike_ind[i]) + (logprior_prop - logprior_ind[i]);
      
      double acceptance_linear = 1.0;
      bool accept_move = true;
      if (MH < 0) {
        acceptance_linear = exp(MH);
        accept_move = (runif_0_1() < acceptance_linear);
      }
      
      // accept or reject move
      if (accept_move) {
        time_inf[i][j] = time_inf_prop[j];
        loglike += loglike_prop - loglike_ind[i];
        logprior += logprior_prop - logprior_ind[i];
        loglike_ind[i] = loglike_prop;
        logprior_ind[i] = logprior_prop;
      } else {
        time_inf_prop[j] = time_inf[i][j];
      }
    }
    
  } // end loop through individuals
  */
}

