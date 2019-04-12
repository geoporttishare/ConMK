#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List c_mann_kendall_test(NumericVector x) {
  int i, j, n = x.size();
  //IntegerVector ties(n); // ADD THIS
  int S=0;
  double s2, d;
  for(i = 0; i < n - 1; i++) { //if(!R_isnancpp(x(i)))
    for(j = i + 1; j < n; j++) {
      d = x(j) - x(i);
      S += d > 0 ? 1 : d < 0 ? -1 : 0;
    }
  }
  s2 = n * (n-1.0)*(2.0*n + 5.0) / 18.0;
  return List::create(Named("S")=(double)S, Named("s2")=s2);
}

// Extend the test to include Kendall's rho-coefficient estimation

// [[Rcpp::export]]
List c_mann_kendall_test_and_beta(NumericVector x, NumericVector t) {
  int i, j, n = x.size();
  // IntegerVector ties(n); //ADD THIS TODO @HERE
  NumericVector X(n*(n-1)/2);
  int S=0;
  double s2, d;
  int k = 0;
  for(i = 0; i < n - 1; i++)
    for(j = i + 1; j < n; j++){
      d = x(j) - x(i);
      S += d > 0 ? 1 : d < 0 ?-1 : 0; // for trend test
      X(k) = d / (t(j)-t(i)); // for estimation, use R's median on these.
      k++;
    }
  s2 = n * (n-1.0)*(2.0*n + 5.0) / 18.0;

  return List::create(Named("S")=(double)S, Named("s2")=s2, Named("X") = X);
}
