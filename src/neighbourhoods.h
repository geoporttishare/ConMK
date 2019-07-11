#include <Rcpp.h>
using namespace Rcpp;

IntegerVector neighbour_cells_none(int i, int nrow, int ncol);

IntegerVector neighbour_cells_queen_row_col(int i, int nrow, int ncol);
IntegerVector neighbour_cells_queen_col_row(int i, int nrow, int ncol);

IntegerVector neighbour_cells_rook_row_col(int i, int nrow, int ncol);
IntegerVector neighbour_cells_rook_col_row(int i, int nrow, int ncol);

IntegerVector forward_neighbour_cells_queen_row_col(int i, int nrow, int ncol);
IntegerVector forward_neighbour_cells_queen_col_row(int i, int nrow, int ncol);

IntegerVector forward_neighbour_cells_rook_row_col(int i, int nrow, int ncol);
IntegerVector forward_neighbour_cells_rook_col_row(int i, int nrow, int ncol);
