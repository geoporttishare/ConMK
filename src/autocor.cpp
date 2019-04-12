#include <Rcpp.h>
using namespace Rcpp;


// lag one autocorrelation
// [[Rcpp::export]]
NumericVector c_autocorrelation_1(NumericVector x) {
  int n = x.size(), i;
  double Sx = 0, Sxx = 0, Sxy = 0, Sxpy = 0;

  for(i = 0; i < n-1; i++) {
    Sx += x(i);
    Sxx += x(i)*x(i);
    Sxy += x(i)*x(i+1);
    Sxpy += x(i)+x(i+1);
  }
  Sx += x(i);
  Sxx += x(i)*x(i);
  //
  double nn = (double)n;
  double mu = Sx/nn;
  double s2 = (Sxx - 2.0*mu*Sx + mu*mu*nn)/(nn-1.0);
  double rho1 = (Sxy - mu*Sxpy + mu*mu*(nn-1.0) ) / ( (nn-1.0) * s2 );
  return NumericVector::create(mu, s2, rho1);
}

