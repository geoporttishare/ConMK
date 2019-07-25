#include <Rcpp.h>
#include "neighbourhoods.h"
using namespace Rcpp;

// We assume the values are centered timewise.

// Upper tri only
// [[Rcpp::export]]
List c_stack_covariance_sparse_UT(NumericMatrix x, //  x is matrix (nrow*ncol) x ntime
                                  int nrow, // this is raster nrow
                                  int neighbourhood = 2  // 2 queen, 1 rook, 0 none
)
{
  int N = x.nrow();
  int ncol = N/nrow; // raster nrow, not x.
  int ntime = x.ncol();
  // neighbourhood
  IntegerVector (*neighbour_cells)(int, int, int);
  IntegerVector nn(4+4+1); //cell count contributions to upper-tri cov-matrix
    //                          in ul corner, ur corner, br corner, bl corner, left edge, upper edge, right edge, bottom edge, interior
  if(neighbourhood == 0){ neighbour_cells = &neighbour_cells_none;          nn(0)=1; nn(1)=1; nn(2)=1; nn(3)=1; nn(4)=1; nn(5)=1; nn(6)=1; nn(7)=1; nn(8)=1;}
  if(neighbourhood == 1){ neighbour_cells = &forward_neighbour_cells_rook_col_row;  nn(0)=3; nn(1)=2; nn(2)=1; nn(3)=2; nn(4)=3; nn(5)=3; nn(6)=2; nn(7)=2; nn(8)=3;}
  if(neighbourhood == 2){ neighbour_cells = &forward_neighbour_cells_queen_col_row; nn(0)=4; nn(1)=2; nn(2)=1; nn(3)=3; nn(4)=5; nn(5)=4; nn(6)=2; nn(7)=3; nn(8)=5;}

  // need to compute the total number of neighbours:
  int M = nn(0) + nn(1) + nn(2) + nn(3) + (nrow-2) * (nn(4)+nn(6)) + (ncol-2)*(nn(5) + nn(7)) + (nrow-2)*(ncol-2)*nn(8);
  NumericVector Cov(M);
  IntegerVector ivec(M), jvec(M);

  IntegerVector nei;
  int i, k1, j, t, it = 0;
  double ss;

  for(i = 0; i < N; i++) { //covs is  Smean, sd, nneigh
    if(!R_isnancpp(x(i,0))) { // make sure not dealing with a missing value cell
      // variance
      ss = 0;
      for(t = 0; t < ntime; t++) ss+= x(i,t)*x(i,t);
      Cov(it) = ss/ntime;
      ivec(it) = i;
      jvec(it) = i;
      // covariances with neighbours.
      it++;
      nei = neighbour_cells(i, nrow, ncol);
      for(k1 = 0; k1 < nei.size(); k1++) { // do slowly. Optimise later
        j = nei(k1);
        if(!R_isnancpp(x(j,0))){
          ivec(it) = i;
          jvec(it) = j;
          ss = 0;
          for(t=0; t < ntime; t++ ) ss += x(j,t)*x(i,t);
          ss /= ntime;
          Cov(it) = ss;
          it++;
        }
      }
    }
  }
  // TODO: turn this into a sparse matrix
  return List::create(Named("i")=ivec, Named("j")=jvec, Named("x")=Cov);
}




// [[Rcpp::export]]
List c_stack_covariance_sparse(NumericMatrix x, //  x is matrix (nrow*ncol) x ntime
                               int nrow, // this is raster nrow
                               int neighbourhood = 2  // 2 queen, 1 rook, 0 none
)
{
  int N = x.nrow();
  int ncol = N/nrow; // raster nrow, not x.
  int ntime = x.ncol();
  // neighbourhood
  IntegerVector (*neighbour_cells)(int, int, int);
  IntegerVector nn(3); //cell in corner, edge, interior, number of contributions to upper-tri cov-matrix.
  if(neighbourhood == 0){ neighbour_cells = &neighbour_cells_none;          nn(0)=1; nn(1)=1; nn(2)=1; }
  if(neighbourhood == 1){ neighbour_cells = &neighbour_cells_rook_col_row;  nn(0)=3; nn(1)=4; nn(2)=5;}
  if(neighbourhood == 2){ neighbour_cells = &neighbour_cells_queen_col_row; nn(0)=1; nn(1)=6; nn(2)=9;}

  // need to compute the total number of neighbours:
  int M = 4 * nn(0) + 2*(nrow-2 + ncol-2)*nn(1) + (nrow-2)*(ncol-2)*nn(2);
  NumericVector Cov(M);
  IntegerVector ivec(M), jvec(M);

  IntegerVector nei;
  int i, k1, j, t, it = 0;
  double ss;

  for(i = 0; i < N; i++) { //covs is  Smean, sd, nneigh
    if(!R_isnancpp(x(i,0))) { // make sure not dealing with a missing value cell
      // variance
      ss = 0;
      for(t = 0; t < ntime; t++) ss+= x(i,t)*x(i,t);
      Cov(it) = ss/ntime;
      ivec(it) = i;
      jvec(it) = i;
      // covariances with neighbours.
      it++;
      nei = neighbour_cells(i, nrow, ncol);
      for(k1 = 0; k1 < nei.size(); k1++) { // do slowly. Optimise later
        j = nei(k1);
        if(!R_isnancpp(x(j,0))){
          ivec(it) = i;
          jvec(it) = j;
          ss = 0;
          for(t=0; t < ntime; t++ ) ss += x(j,t)*x(i,t);
          ss /= ntime;
          Cov(it) = ss;
          it++;
        }
      }
    }
  }
  // TODO: turn this into a sparse matrix
  return List::create(Named("i")=ivec, Named("j")=jvec, Named("x")=Cov);
}




// Full matrix calculation. Massively memoryhungry
// [[Rcpp::export]]
List c_stack_covariance(NumericMatrix x, //  x is matrix (nrow*ncol) x ntime
                        int nrow, // this is raster nrow
                        int neighbourhood = 2  // 2 queen, 1 rook, 0 none
    )
{
  int N = x.nrow();
  int ncol = N/nrow; // raster nrow, not x.
  int ntime = x.ncol();
  NumericMatrix Cov(N, N);
  // neighbourhood
  IntegerVector (*neighbour_cells)(int, int, int);
  if(neighbourhood == 0) neighbour_cells = &neighbour_cells_none;
  if(neighbourhood == 1) neighbour_cells = &neighbour_cells_rook_col_row;
  if(neighbourhood == 2) neighbour_cells = &neighbour_cells_queen_col_row;
  IntegerVector nei;
  int i, k1, j, t;
  double ss;

  for(i = 0; i < N; i++) { //covs is  Smean, sd, nneigh
    if(!R_isnancpp(x(i,0))) { // make sure not dealing with a missing value cell
      // variance
      ss = 0;
      for(t = 0; t < ntime; t++) ss+= x(i,t)*x(i,t);
      Cov(i,i) = ss/ntime;
      // covariances with neighbours
      nei = neighbour_cells(i, nrow, ncol);
      for(k1 = 0; k1 < nei.size(); k1++) { // do slowly. Optimise later
        j = nei(k1);
        if(!R_isnancpp(x(j,0))){
          ss = 0;
          for(t=0; t < ntime; t++ ) ss += x(j,t)*x(i,t);
          ss /= ntime;
          Cov(i,j) = ss;
          Cov(j,i) = ss;
        }
      }
    }
  }
  // TODO: turn this into a sparse matrix
  return List::create(Named("Cov")=Cov);
}

