#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
bool bolean_test(double x, double y){
  return (x>5) && (y<10);
}

// Assing a vector due for
//[[Rcpp::export]]
NumericVector new_vec(NumericVector old_vec){

    NumericVector new_vector;
    int vec_size = old_vec.size();
    for(int i=0;i<vec_size;i++){
      new_vector.push_back(old_vec[i]);
    }

    for(int j = 0; j<vec_size;j++){
      new_vector.erase(j);
    }
    return new_vector;
}

//[[Rcpp::export]]
NumericVector filling_vec(NumericVector x){

// NEED TO SELECT ONLY SPLIT RULE FOR THAT TREE
  NumericVector x_curr(x.size());

  for(int i=0;i<x.size();i++){
    x_curr(i) = x[i];
  }

  return x_curr;
}
