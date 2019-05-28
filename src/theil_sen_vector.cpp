#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector c_theil_sen_vector(NumericVector x,
                                 NumericVector t) {
  int i,j,k,n = x.size();
  NumericVector w( n * (n-1)/2);
  k = 0;
  for(i=0; i < n-1; i++) {
    for(j = i+1; j < n; j++) {
      w(k) = (x(j)-x(i)) /(t(j) - t(i));
      k++;
    }
  }
  return w;
}

