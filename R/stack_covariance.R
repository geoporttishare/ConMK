#' Covariance Matrix of Timeseries RasterStack
#'
#' Compute the covariance of two timeseries in neighbouring cells.
#'
#' @param x Raster stack.
#' @param neighbourhood 0: none 1:rook 2:queen (default is 2)
#' @param sparse Logical: Use spareMatrix? (default TRUE, highly recommended)
#' @param ver leave at 2. (internal development only.)
#' @import Matrix
#'
#' @details Compute the covariances between timeserieses in rasterStack cells, limited to neighbouring cells.
#'
#' @export

stack_covariance <- function(x, neighbourhood = 2, sparse = TRUE, ver = 2) {
  if(!neighbourhood %in% c(0,1,2)) stop("Only '0', '1' and '2' neighbourhoods implemented.")
  # center first
  if( !canProcessInMemory(x) ) stop("'canProcessInMemory' return FALSE")
  t0 <- Sys.time()
  x <- x - mean(x)
  X <- values(stack(x))
  nr <- nrow(x)
  nc <- ncol(x)
  N <- nc * nr
  #
  if(!sparse) c_stack_covariance(X, nr, neighbourhood)[[1]]
  else{
    if(ver == 1) {
    ijx <- c_stack_covariance_sparse(X, nr, neighbourhood)
    sparseMatrix(i = ijx$i, j = ijx$j, x = ijx$x, symmetric = FALSE, index1 = FALSE, dims = c(N, N))
    }else{
      ijx <- c_stack_covariance_sparse_UT(X, nr, neighbourhood)
      sparseMatrix(i = ijx$i, j = ijx$j, x = ijx$x, symmetric = TRUE, index1 = FALSE, dims = c(N, N))
    }
  }
}
