#include<iostream>
#include <vector>
#include <math.h>
#include<cmath>
#include<RcppArmadillo.h>
#include <algorithm>

using namespace std;
using namespace Rcpp;

// Function to sample an integer from a sequence from
// 0:n
int sample_int(int n){
  return rand() % n;
}

// Sample an uniform value from a double vector
//[[Rcpp::export]]
double sample_double(Rcpp::NumericVector vec, int n_min_size){

  // Getting the range of the vector
  std::sort(vec.begin(),vec.end() );
  int vec_size = vec.size();
  Rcpp::NumericVector valid_values = vec[Range((n_min_size-1),(vec_size-n_min_size-1))];
  int sample_index = rand() % valid_values.size();

  return valid_values(sample_index);
}

// Sample an uniform value from a double vector
// // [[Rcpp::export]]
// double sample_double(arma::vec vec, int n_min_size){
//
//   // Getting the range of the vector
//   int sample_index = floor(R::runif(0,vec.size()));
//   return vec(sample_index);
// }

