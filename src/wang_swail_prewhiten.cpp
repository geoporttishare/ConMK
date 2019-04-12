#include <Rcpp.h>
#include <R.h>
using namespace Rcpp;


double autocor_part(NumericVector x, double beta) {
  int n = x.size(), i;
  double Sx = 0, Sxx = 0, Sxy = 0, Sxpy = 0;
  double a, b;
  for(i = 0; i < n-1; i++) {
    a = x(i) - (i+1.0)*beta;
    b = x(i+1) - (i+2.0)*beta;
    Sx += a;
    Sxx += a*a;
    Sxy += a*b;
    Sxpy += a+b;
  }
  a = x(i) - (i+1.0)*beta;
  Sx += a;
  Sxx += a*a;
  //
  double nn = (double)n;
  double npairs = nn * (nn-1.0)/2.0;
  double mu = Sx/nn;
  double s2 = (Sxx - 2.0*mu*Sx + mu*mu*nn)/(nn-1.0);
  double rho1 = (Sxy - mu*Sxpy + mu*mu*(nn-1.0) ) / ( (nn-1.0) * s2 );
  return rho1;
}

double beta_part(NumericVector x) {
  int i, j, n = x.size();
  // IntegerVector ties(n); //ADD THIS TODO @HERE
  NumericVector X(n*(n-1)/2);
  int S=0;
  double s2, d;
  int k = 0;
  for(i = 0; i < n - 1; i++)
    for(j = i + 1; j < n; j++){
      d = x(j) - x(i);
      X(k) = d / (j-i);
      k++;
    }
  return median(X);
}

double abs(double a){
  if( a < 0) return -a;
  return a;
}


// Wang&Swail prewhitening assuming AR(1) noise. Equidistant data, t ignored below.


// [[Rcpp::export]]
List c_wang_swail_prewithen_1d(NumericVector x, NumericVector t,
                            double eps, double rho_th, int itmax) {
  int i,j, it = 0, n = x.size();

  NumericVector W(n-1);

  double rho0 = 2;
  double rho = 0;
  double beta = 0, beta0 = 0;

  while( (abs(rho0 - rho) > eps  || abs(beta0 - beta) > eps) && it < itmax) {
    beta0 = beta;
    rho0 = autocor_part(x, beta0);
    for(i=0; i < n-1; i++) W(i) = (x(i+1)-rho0*x(i))/(1-rho0);
    beta0 = beta_part(W);
    rho = autocor_part(x, beta0);
    for(i=0; i < n-1; i++) W(i) = (x(i+1)-rho*x(i))/(1-rho);
    beta = beta_part(W);
    it++;
  }
  //if(it == itmax) Rprintf("itmax reached.\n");
  return List::create(Named("W")=W, Named("rho")=rho, Named("beta")=beta, Named("iter")=it);
}


