#include<iostream>
#include <vector>
#include<Rcpp.h>
#include <math.h>
#include<cmath>
#include <algorithm>

using namespace std;

// Function to sample an integer from a sequence from
// 0:n
//[[Rcpp::export]]
int sample_int(int n){
  return rand() % n;
}

// // Sample an uniform value from a double vector
// //[[Rcpp::export]]
// double sample_double(Rcpp::NumericVector vec, int n_min_size){
//
//   // Getting the range of the vector
//   std::sort(vec.begin(),vec.end() );
//   int vec_size = vec.size();
//   Rcpp::NumericVector valid_values = vec[Rcpp::Range((n_min_size-1),(vec_size-n_min_size-1))];
//   int sample_index = rand() % valid_values.size();
//
//   return valid_values(sample_index);
// }

// Sample an uniform value from a double vector
// [[Rcpp::export]]
double sample_double(Rcpp::NumericVector vec){

  // Getting the range of the vector
  int sample_index = floor(R::runif(0,vec.size()));
  // cout << "Sample Index:  " << vec.size() << endl;

  return vec(sample_index);
}

// [[Rcpp::export]]
double sample_rule(Rcpp::NumericVector x_cut_var,
                   Rcpp:: NumericVector x_obs_var,
                   int node_min_size){

  // Getting the range of the vector
  // Rcpp::NumericVector x_clone = clone(x_obs_var);
  // std::sort(x_clone.begin(),x_clone.end() );
  double min_v = min(x_obs_var);
  double max_v = max(x_obs_var);
  // cout << "Min value " << min_v << endl;
  // cout << "Max value " << max_v << endl;

  // Vector with only valid x_cut_values
  Rcpp::NumericVector x_cut_f;

  // Filtering only the valid cut values
  for(int i=0;i<x_cut_var.size();i++){
    // cout << x_cut_var(i) << endl ;
    if(((x_cut_var(i))>=(min_v)) & ((x_cut_var(i))<=max_v)){
      x_cut_f.push_back(x_cut_var(i));
    }
  }

  if(x_cut_f.size() == 0){
    return 0.0003;
  }

  return x_cut_f(sample_int(x_cut_f.size()));
}


