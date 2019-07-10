#include "neighbourhoods.h"

/*  indices of neighbouring cells of i. 0-base   */

// values assumed to be
//   row-col
// order in the matrix
//[[Rcpp::export]]
IntegerVector neighbour_cells_queen_row_col(int i, int nrow, int ncol) {
  int ic = (int) i/nrow;
  int ir = i-ic*nrow;
  IntegerVector n;

  if(ir == 0) { // at first row
    if(ic == 0) { // at first col
      n.push_back(1); // below
      n.push_back(nrow);// right
      n.push_back(nrow+1);// b-r
    }
    else if(ic == ncol-1){ //at last col
      n.push_back(i+1);// below
      n.push_back(i-nrow);// left
      n.push_back(i+1-nrow); // l-b
    }
    else{ // in between on the first row
      n.push_back(i+nrow);// right
      n.push_back(i+nrow+1); // r-b
      n.push_back(i+1);// below
      n.push_back(i+1-nrow); // l-b
      n.push_back(i-nrow);// left
    }
  }
  else if(ir == nrow-1){ // at the last row
    if(ic == 0) { // at first col. i==ir
      n.push_back(i-1); // above
      n.push_back(i+nrow);// right
      n.push_back(i+nrow-1);// a-r
    }
    else if(ic == ncol-1){ //at last col
      n.push_back(i-nrow);// left
      n.push_back(i-1);// above
      n.push_back(i-1-nrow); //a-l
    }
    else{ // in between on the last row
      n.push_back(i-nrow);// left
      n.push_back(i-1-nrow); //a-l
      n.push_back(i-1);// above
      n.push_back(i-1+nrow);// a-r
      n.push_back(i+nrow);// right
    }
  }
  else if(ic == 0){ // left edge, but not corner
    n.push_back(i-1); // above
    n.push_back(i+1); // below
    n.push_back(i-1+nrow);// a-r
    n.push_back(i  +nrow);// right
    n.push_back(i+1+nrow);// b-r
  }
  else if(ic == ncol-1){ // right edge, but not corner
    n.push_back(i-1); // above
    n.push_back(i+1); // below
    n.push_back(i-nrow);// left
    n.push_back(i-1-nrow); //a-l
    n.push_back(i+1-nrow); // l-b
  }
  else{ // not on any edge
    n.push_back(i-nrow);// left
    n.push_back(i-nrow-1); //a-l
    n.push_back(i-1); // above
    n.push_back(i+nrow-1);// a-r
    n.push_back(i+nrow);// right
    n.push_back(i+nrow+1);// b-r
    n.push_back(i+1); // below
    n.push_back(i-nrow+1); // l-b

  }
  return n;
}


// indices of neighbouring cells of i. 0-base
// values assumed to be
//   col-row
// order in the matrix
//[[Rcpp::export]]
IntegerVector neighbour_cells_queen_col_row(int i, int nrow, int ncol) {
  int ir = (int) i/ncol;
  int ic = i-ir*ncol;
  IntegerVector n;

  if(ir == 0) { // at first row
    if(ic == 0) { // at first col
      n.push_back(ncol); // below
      n.push_back(1);// right
      n.push_back(ncol+1);// b-r
    }
    else if(ic == ncol-1){ //at last col
      n.push_back(ncol+i);// below
      n.push_back(i-1);// left
      n.push_back(ncol+i-1); // l-b
    }
    else{ // in between on the first row
      n.push_back(i+1);// right
      n.push_back(i+ncol+1); // r-b
      n.push_back(i+ncol);// below
      n.push_back(i+ncol-1); // l-b
      n.push_back(i-1);// left
    }
  }
  else if(ir == nrow-1){ // at the last row
    if(ic == 0) { // at first col
      n.push_back(i-ncol); // above
      n.push_back(i+1);// right
      n.push_back(i-ncol+1);// a-r
    }
    else if(ic == ncol-1){ //at last col
      n.push_back(i-1);// left
      n.push_back(i-ncol);// above
      n.push_back(i-ncol-1); //a-l
    }
    else{ // in between on the last row
      n.push_back(i-1);// left
      n.push_back(i-ncol-1); //a-l
      n.push_back(i-ncol);// above
      n.push_back(i-ncol+1);// a-r
      n.push_back(i+1);// right
    }
  }
  else if(ic == 0){ // left edge, but not corner
    n.push_back(i-ncol); // above
    n.push_back(i+ncol); // below
    n.push_back(i-ncol+1);// a-r
    n.push_back(i+1);// right
    n.push_back(i+ncol+1);// b-r
  }
  else if(ic == ncol-1){ // right edge, but not corner
    n.push_back(i-ncol); // above
    n.push_back(i+ncol); // below
    n.push_back(i-1);// left
    n.push_back(i-ncol-1); //a-l
    n.push_back(i+ncol-1); // l-b
  }
  else{ // not on any edge
    n.push_back(i-1);// left
    n.push_back(i-ncol-1); //a-l
    n.push_back(i-ncol); // above
    n.push_back(i-ncol+1);// a-r
    n.push_back(i+1);// right
    n.push_back(i+ncol+1);// b-r
    n.push_back(i+ncol); // below
    n.push_back(i+ncol-1); // l-b

  }
  return n;
}






/*  indices of forward neighbouring cells of i. 0-base. These are
 *  retrieving neighbours with j > i.
 */

// values assumed to be
//   row-col
// order in the matrix
//[[Rcpp::export]]
IntegerVector forward_neighbour_cells_queen_row_col(int i, int nrow, int ncol) {
  int ic = (int) i/nrow;
  int ir = i-ic*nrow;
  IntegerVector n;

  if(ir == 0) { // at first row
    if(ic == 0) { // at first col
      n.push_back(1); // below
      n.push_back(nrow);// right
      n.push_back(nrow+1);// b-r
    }
    else if(ic == ncol-1){ //at last col
      n.push_back(i+1);// below
    }
    else{ // in between on the first row
      n.push_back(i+nrow);// right
      n.push_back(i+nrow+1); // r-b
      n.push_back(i+1);// below
    }
  }
  else if(ir == nrow-1){ // at the last row
    if(ic == 0) { // at first col. i==ir
      n.push_back(i+nrow);// right
      n.push_back(i+nrow-1);// a-r
    }
    else if(ic == ncol-1){ //at last col
    }
    else{ // in between on the last row
      n.push_back(i-1+nrow);// a-r
      n.push_back(i+nrow);// right
    }
  }
  else if(ic == 0){ // left edge, but not corner
    n.push_back(i+1); // below
    n.push_back(i-1+nrow);// a-r
    n.push_back(i  +nrow);// right
    n.push_back(i+1+nrow);// b-r
  }
  else if(ic == ncol-1){ // right edge, but not corner
    n.push_back(i+1); // below
  }
  else{ // not on any edge
    n.push_back(i+nrow-1);// a-r
    n.push_back(i+nrow);// right
    n.push_back(i+nrow+1);// b-r
    n.push_back(i+1); // below
  }
  return n;
}


// indices of neighbouring cells of i. 0-base
// values assumed to be
//   col-row
// order in the matrix
//[[Rcpp::export]]
IntegerVector forward_neighbour_cells_queen_col_row(int i, int nrow, int ncol) {
  int ir = (int) i/ncol;
  int ic = i-ir*ncol;
  IntegerVector n;

  if(ir == 0) { // at first row
    if(ic == 0) { // at first col
      n.push_back(ncol); // below
      n.push_back(1);// right
      n.push_back(ncol+1);// b-r
    }
    else if(ic == ncol-1){ //at last col
      n.push_back(ncol+i);// below
      n.push_back(ncol+i-1); // l-b
    }
    else{ // in between on the first row
      n.push_back(i+1);// right
      n.push_back(i+ncol+1); // r-b
      n.push_back(i+ncol);// below
      n.push_back(i+ncol-1); // l-b
    }
  }
  else if(ir == nrow-1){ // at the last row
    if(ic == 0) { // at first col
      n.push_back(i+1);// right
    }
    else if(ic == ncol-1){ //at last col
    }
    else{ // in between on the last row
      n.push_back(i+1);// right
    }
  }
  else if(ic == 0){ // left edge, but not corner
    n.push_back(i+ncol); // below
    n.push_back(i+1);// right
    n.push_back(i+ncol+1);// b-r
  }
  else if(ic == ncol-1){ // right edge, but not corner
    n.push_back(i+ncol); // below
    n.push_back(i+ncol-1); // l-b
  }
  else{ // not on any edge
    n.push_back(i+1);// right
    n.push_back(i+ncol+1);// b-r
    n.push_back(i+ncol); // below
    n.push_back(i+ncol-1); // l-b

  }
  return n;
}




IntegerVector neighbour_cells_none(int i, int nrow, int ncol) {
  IntegerVector n(0);
  return n;
}
