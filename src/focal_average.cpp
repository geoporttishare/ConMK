#include <Rcpp.h>
#include "neighbourhoods.h"
using namespace Rcpp;

// Local mean of each cell in a raster layer of a stack
// [[Rcpp::export]]
NumericMatrix c_focal_average(NumericMatrix x,
                              int nrow,
                              int neighbourhood = 2){
  int i, j, k, t, N = x.nrow();
  int ncol = N/nrow;
  int ntime = x.ncol();
  NumericMatrix out(N, ntime); // mean
  IntegerVector nn(N);
  NumericVector nei;
  // neighbourhood
  IntegerVector (*forward_neighbour_cells)(int, int, int);
  if(neighbourhood == 0) forward_neighbour_cells = &neighbour_cells_none;
  if(neighbourhood == 1) forward_neighbour_cells = &forward_neighbour_cells_rook_col_row;
  if(neighbourhood == 2) forward_neighbour_cells = &forward_neighbour_cells_queen_col_row;
  //
  for(t = 0; t < ntime; t++) {
    for(i = 0; i < N; i++)
      if( !R_isnancpp( x(i,t) ) )  {
        out(i,t) += x(i,t);
        nei = forward_neighbour_cells(i, nrow, ncol);
        if( nei.size() ) {
          for(k = 0; k < nei.size(); k++) {
            j = nei(k);
            if(!R_isnancpp( x(j,t) )) {
              out(i,t) += x(j,t);
              out(j,t) += x(i,t);
              nn(i)++;
              nn(j)++;
            }
          }
        }
        out(i,t) /=  1.0 + nn(i);
      }
      else out(i,t) = R_NaReal;
  }
  return out; // mean per time
}
