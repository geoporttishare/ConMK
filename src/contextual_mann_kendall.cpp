#include <Rcpp.h>

#include "neighbourhoods.h"

using namespace Rcpp;

// We assume the values are centered timewise.
// compute local mean and cell variance. Covariances elsewhere.
NumericMatrix local_stats(NumericMatrix x, int nrow,
                           NumericVector MKS, int neighbourhood = 2){
  int i, j, k, t, N = x.nrow();
  int ncol = N/nrow;
  int ntime = x.ncol();
  NumericMatrix out(N, 3); // mean_S, sd, nneigh
  NumericVector nei;

  // neighbourhood
  IntegerVector (*forward_neighbour_cells)(int, int, int);
  if(neighbourhood == 2) forward_neighbour_cells = &forward_neighbour_cells_queen_col_row;
  if(neighbourhood == 0) forward_neighbour_cells = &neighbour_cells_none;
  //double ss = 0;
  //
  for(i = 0; i < N; i++)
  if( !R_isnancpp( x(i,0) ) )  {
    out(i,0) += MKS(i);
    for(t = 0; t < ntime; t++)  out(i,1) += x(i, t) * x(i, t);
    out(i,1) = sqrt( out(i,1) / (ntime-1.0) );
    nei = forward_neighbour_cells(i, nrow, ncol);
    if( nei.size() ) {
        for(k = 0; k < nei.size(); k++) {
          j = nei(k);
          if(!R_isnancpp( x(j,0) )) {
            out(i,0) += MKS(j);
            out(j,0) += MKS(i);
            out(i,2)++;
            out(j,2)++;
          }
        }
    }
    out(i, 0) /=  1.0 + out(i,2); // mean S
  }
  else for(k=0; k < out.ncol(); k++) out(i,k) = R_NaReal;

  return out; // Smean, var, nneig
}

/* Mann-Kendall S-statistic */
// [[Rcpp::export]]
NumericVector mann_kendall_S(NumericMatrix x){
  double d;
  int i,j,k;
  double S;
  int ntime = x.ncol();
  NumericVector out(x.nrow());
  for(i = 0; i < x.nrow(); i++) if(!R_isnancpp(x(i,0))){
    S = 0;
    for(j = 0; j < ntime; j++)
      for(k = j+1; k < ntime; k++){
        d = x(i,k) - x(i,j);
        S += d > 0 ? 1 : d < 0 ? -1 : 0;
      }
    out(i) = S;
  }else out(i) = R_NaReal;
  return out;
}


// we assume x[,i] is given in col-row order, raster does this

// [[Rcpp::export]]
List c_contextual_mann_kendall(NumericMatrix x, int nrow,
                               int neighbourhood = 2)  // 2 queen, 0 none
  {
  int N = x.nrow();
  int ncol = N/nrow;
  int ntime = x.ncol();
  // stastic at each cell
  NumericVector MKS = mann_kendall_S(x);
  // compute covs for the covariance formula
  NumericMatrix covs = local_stats(x, nrow, MKS, neighbourhood);

  // we have all components now. Let's compute the the test statistic
  // and variance
  NumericMatrix S_and_s2(N,2);
  NumericVector nei;
  int i,j1,j2,k1,k2,t,m;
  double ss, vs, Sm;
  double s2 = ntime * (ntime-1.0)*(2.0*ntime + 5.0) / 18.0;

  // neighbourhood
  IntegerVector (*neighbour_cells)(int, int, int);
  if(neighbourhood == 2) neighbour_cells = &neighbour_cells_queen_col_row;
  if(neighbourhood == 0) neighbour_cells = &neighbour_cells_none;

  for(i = 0; i < N; i++) { //covs is  Smean, sd, nneigh
    Sm = covs(i, 0);
    if(!R_isnancpp(Sm) ){
      m = 1;
      vs = s2;//covs(i,1) * covs(i,1);
      nei = neighbour_cells(i, nrow, ncol);
      for(k1 = 0; k1 < nei.size(); k1++) {
        j1 = nei(k1);
        if(!R_isnancpp(covs(j1,1))){
          m++;
          vs += s2;//covs(j1,1) * covs(j1,1);
          // correlation with focal
          ss = 0;
          for(t=0; t < ntime; t++ ) ss += x(j1,t)*x(i,t);
          ss /= ntime * covs(j1,1) * covs(i,1);
          vs += 2 * ss * s2;
          // need also correlation between other neighs
          for(k2 = k1+1; k2 < nei.size(); k2++) {
            j2 = nei(k2);
            if(!R_isnancpp(covs(j2,1))){
              ss = 0;
              for(t=0; t < ntime; t++ ) ss += x(j1,t)*x(j2,t);
              ss /= ntime * covs(j1,1) * covs(j2,1);
              vs += 2 * ss * s2;
            }
          }
        }
      }
      vs /= (double) m * m;
      S_and_s2(i, 0) = Sm;
      S_and_s2(i, 1) = vs;
    }else{
      S_and_s2(i,0) = R_NaReal;
      S_and_s2(i,1) = R_NaReal;
    }
  }

  return List::create(covs, MKS, S_and_s2, N, ncol, ntime);
}

