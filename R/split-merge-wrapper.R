#' Split-Calc-Merge Wrapper
#'
#' Takes a raster stack, splits it to pieces, apply fun for each sub-stack, and then merges together.
#'
#' @param x Raster stack
#' @param nx passed on to splitRaster
#' @param ny passed on to splitRaster
#' @param buffer passed on to splitRaster
#' @param fun The function to call per piece. Should return either a raster or stack
#' @param ... passed to fun in argument list after the piece
#' @param dbg Print messages?
#' @param neighbourhood 2:queen, 1:rook, 0:none (classical M-K test)
#'
#' @details The result should be identical to just calling 'fun(x, ...)'.
#'
#' @note If at stitching a message about different buffers appears (from 'mergeRaster()'), the (nx,ny) pair
#' is problematic for 'mergeRaster' and the resulting raster will have artifacts in the splitted edges. Different (nx,ny) is the only way to get around this at the moment. (Should be fixed by devs!)
#'
#' @references
#' Neeti, N. and Eastman J.R. (2011) A Contextual Mann-Kendall Approach for the Assesment of Trend Significance in Image Time Series, \emph{Transactions in GIS}
#' @useDynLib ConMK
#' @import raster progress SpaDES.tools
#' @export
split_calc_wrapper <- function(x, nx, ny, buffer = c(1,1), fun, ..., dbg = FALSE){
  #
  # Something wrong with the splitRaster-mergeRaster of SpaDES.


  # create a list of extents. copied from SpaDES.tools::splitRaster
  wlist <- vector("list", length = nx * ny)
  ext <- extent(x)
  n <- 1
  for (i in seq_len(nx) - 1) {
    for (j in seq_len(ny) - 1) {
      x0 <- ext@xmin +       i * ((ext@xmax - ext@xmin)/nx) - buffer[1] * xres(x)
      x1 <- ext@xmin + (i + 1) * ((ext@xmax - ext@xmin)/nx) + buffer[1] * xres(x)
      y0 <- ext@ymin +       j * ((ext@ymax - ext@ymin)/ny) - buffer[2] * yres(x)
      y1 <- ext@ymin + (j + 1) * ((ext@ymax - ext@ymin)/ny) + buffer[2] * yres(x)
      wlist[[n]] <- extent(x0, x1, y0, y1)
      n <- n + 1L
    }
  }
  #wlist <- splitRaster(raster(x), nx = nx, ny = ny, buffer = buffer)

  # Crop the original stack to each subwindow and calculate.
  if(dbg) message("number of pieces: ", length(wlist))
  estlist <- list()
  T0 <- Sys.time()
  tl <- NULL

  if(dbg) pb <- progress_bar$new(format = "[:bar]:current/:total eta :eta", total = length(wlist))
  if(dbg) message("calculating:")
  #
  # TODO Parallisation here.
  if(dbg) pb$tick(0)
  for(wi in seq_along(wlist)){
    w <- wlist[[wi]]
    #substack <- crop(x = x, y = w)
    estlist[[wi]] <- fun( crop(x = x, y = w, snap = "out") ,  ...   )
    if(dbg) pb$tick()
  }
  t1 <- Sys.time()-T0
  # Stitch: just rasters,
  if(dbg) message("stitching:", appendLF=FALSE)
  #browser()
  t2 <- Sys.time()

  if(!is(estlist[[1]], "RasterStack")) {
    out <- mergeRaster( estlist )

  }else{
    vnams <- names(estlist[[1]])
    # needs optimisation here
    outl <- lapply(vnams, function(v)
      #suppressMessages(
        mergeRaster( lapply(estlist, function(z) z[[v]]) ))
      #) # suppress the message about using mosaic
    out <- stack(outl)
  }
  #
  t2 <- Sys.time() - t2
  #browser()
  if(dbg){ message(". done!"); pb$terminate() }
  #
  list(result = out,
       extents = wlist,
       split_info = list(nx=nx, ny=ny, buffer=buffer),
       timing = list(crop_n_calc = t1, stitch=t2))
}
