#' Split-Calc-Merge Wrapper
#'
#' Takes a raster stack, splits it into pieces, applies a function on each sub-stack, and then merges the results together.
#'
#' @param x Raster stack
#' @param nx Number of x-dimension divisions
#' @param ny Number of y-dimension divisions
#' @param buffer Extent each sub-raster* by this many cells
#' @param fun The function to call per piece. Should return either a raster or stack of identical dimensions as input for it.
#' @param ... passed to fun in argument list after the piece
#' @param dbg Print messages?
#' @param byrowcol Split by row-col indices instead of spatial extents? Recommended to keep this TRUE
#'
#' @details
#' This function automatically takes advantage of raster-package's multicore setup.
#' See \link{beginCluster} for details how to use it.
#'
#' The 'fun' should take a rasterStack as input and also should return a rasterStack
#' (or raster).
#'
#' @note If byrowcol=FALSE: a stitching message about different buffers might appear (from 'mergeRaster()'), the (nx,ny) pair
#' is problematic for 'mergeRaster' and the resulting raster will have artifacts in the splitted edges. Different (nx,ny) is the only way to get around this at the moment.
#'
#' @return The output is a list. The first element of this list, called 'result', should be identical to the output of simply calling 'fun(x, ...)'. The remaining elements contain various information related to the calculation, e.g. time spent.
#'
#' @references
#' Neeti, N. and Eastman J.R. (2011) A Contextual Mann-Kendall Approach for the Assesment of Trend Significance in Image Time Series, \emph{Transactions in GIS}
#' @useDynLib ConMK
#' @import raster progress
#' @export
split_calc_wrapper <- function(x, nx, ny, buffer = c(1,1), fun, ...,
                               dbg = FALSE, byrowcol = TRUE){
  #
  #################################################################################
  # details of the splits
  wlist <- vector("list", length = nx * ny)
  rclist <- vector("list", length = nx * ny)
  if(byrowcol){
    # Crop by row-column based extents
    wlist <- vector("list", length = nx * ny)
    n <- 1
    rc <- dim(x)[-3]
    dx <- ceiling(rc[2]/nx)
    dy <- ceiling(rc[1]/ny)
    for (i in seq_len(nx) - 1) {
      for (j in seq_len(ny) - 1) {
        x0 <- 1 +       i * dx
        x1 <- min(rc[2], 1 + (i + 1) * dx )
        y0 <- 1 +       j * dy
        y1 <- min(rc[1], 1 + (j + 1) * dy )
        rclist[[n]] <- c(y0, y1, x0, x1) # original picture rc without buffer
        wlist[[n]] <- extent(x, max(1,     y0 - buffer[2]),
                                min(rc[1], y1 + buffer[2]),
                                max(1,     x0 - buffer[1]),
                                min(rc[2], x1 + buffer[1]))
        n <- n + 1L
      }
    }
  }
  else{
    #################################################################################
    # OR create a list of extents using spatial intervals.
    # copied from SpaDES.tools::splitRaster.
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
  }
  #
  if(dbg) message("number of pieces: ", length(wlist))
  #
  ##################################################################################
  #
  # Crop the original stack to each subwindow and calculate.
  #
  estlist <- list()
  T0 <- Sys.time()

  if(dbg) pb <- progress_bar$new(format = "[:bar]:current/:total eta :eta", total = length(wlist))
  if(dbg) message("calculating:")
  #

  if(dbg) pb$tick(0)
  #
  # This will compute one piece
  clfun <- function(i) {
    w <- wlist[[i]]
    fun( crop(x = x, y = w) ,  ...   )
  }
  #
  npieces <- length(wlist)
  doparallel <- raster:::.doCluster() & npieces > 1
  ##################################################################
  # TODO Parallisation here.
  # Use low-level snow functions in case we need to do intermediate
  # write-to-disk.
  #
  if(doparallel) {
    cl <- getCluster() # latch to the cluster
    on.exit(returnCluster()) # remember to detach
    #
    nnodes <- length(cl)
    for(i in 1:min(npieces, nnodes) ) snow::sendCall( cl[[i]], clfun, i, tag=i )
    #browser()
    for(i in 1:npieces) {
      d <- snow::recvOneData(cl)
      if(!d$value$success) stop("cluster error")
      b <- d$value$tag
      estlist[[b]] <- d$value$value
      ni <- nnodes + i
      if(ni <= npieces){
        snow::sendCall(cl[[d$node]], clfun, ni, tag = ni)
      }
      if(dbg) pb$tick()
    }
  }
  ##################################################################
  # No parallelisation
  else{
    for(wi in seq_along(wlist)){
      estlist[[wi]] <- clfun(wi)
      if(dbg) pb$tick()
    }
  }
  ##################################################################
  #
  t1 <- Sys.time()-T0
  ##
  # Stitch: just rasters,
  if(dbg) message("stitching:", appendLF=FALSE)
  #browser()
  t2 <- Sys.time()
  r0 <- raster(x)
  #browser()
  merger <- if(byrowcol) mergeRaster_cell else mergeRaster
  if(!is(estlist[[1]], "RasterStack")) {
    out <- merger( x = estlist, rclist=rclist, r0=r0, buffer = buffer )

  }else{
    vnams <- names(estlist[[1]])
    # needs optimisation here
    outl <- lapply(vnams, function(v)
      #suppressMessages(
      merger( x = lapply(estlist, function(z) z[[v]]), rclist=rclist, r0=r0, buffer = buffer )
      )
      #) # suppress the message about using mosaic
    out <- stack(outl)
    names(out) <- vnams
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


#' Merge RasterLayers
#'
#' Copied from SpaDES.tools which has large requirements.
#'
#' @param x list of rasters
#' @param fun function to apply on the overlapping regions
#'


mergeRaster <- function (x, ..., fun = NULL)
{
  xminExtent <- sort( unique(  sapply(x, xmin) ))  #%>% unique() %>% sort()
  xmaxExtent <- sort( unique(  sapply(x, xmax) )) #%>% unique() %>% sort()
  yminExtent <- sort( unique(  sapply(x, ymin) )) #%>% unique() %>% sort()
  ymaxExtent <- sort( unique(  sapply(x, ymax) )) #%>% unique() %>% sort()
  xBuffer <- unique((xmaxExtent[-length(xmaxExtent)] - xminExtent[-1])/2)
  yBuffer <- unique((ymaxExtent[-length(ymaxExtent)] - yminExtent[-1])/2)
  if (any(length(xBuffer) > 1, length(yBuffer) > 1)) {
    message(paste0("The tiles present different buffers (likely due to resampling).",
                   " mergeRaster() will use raster::mosaic()."))
    for (i in seq_along(x)) {
      if (i == 1)
        next
      rTemplate <- x[[1]]
      raster::extent(x[[i]]) <- raster::alignExtent(extent = raster::extent(x[[i]]),
                                                    object = rTemplate, snap = "near")
    }
    rasMosaicArgs <- x
    if (!is.null(fun)) {
      rasMosaicArgs$fun <- fun
    }
    else {
      rasMosaicArgs$fun <- mean
    }
    y <- do.call(what = raster::mosaic, args = rasMosaicArgs)
  }
  else {
    for (i in seq_along(x)) {
      r <- x[[i]]
      if (xmin(r) != min(xminExtent)) {
        xminCut <- xmin(r) + xBuffer
      }
      else {
        xminCut <- xmin(r)
      }
      if (xmax(r) != max(xmaxExtent)) {
        xmaxCut <- xmax(r) - xBuffer
      }
      else {
        xmaxCut <- xmax(r)
      }
      if (ymin(r) != min(yminExtent)) {
        yminCut <- ymin(r) + yBuffer
      }
      else {
        yminCut <- ymin(r)
      }
      if (ymax(r) != max(ymaxExtent)) {
        ymaxCut <- ymax(r) - yBuffer
      }
      else {
        ymaxCut <- ymax(r)
      }
      x[[i]] <- crop(r, extent(xminCut, xmaxCut, yminCut,
                               ymaxCut))
    }
    y <- do.call(raster::merge, x)
  }
  names(y) <- gsub("_tile[0-9].*$", "", names(x[[1]]))
  return(y)
}

#' Merge Rasters Split by row-col
#'
#' @param rlist list of rasters, split by row-col
#' @param rclist row-col data of the split
#' @param r0 raster with original details
#'
#' @details No checks.
#'
mergeRaster_cell <- function(x, rclist, r0, buffer) {
  z <- dim(r0)
  V <- matrix( NA, nrow = z[1], ncol = z[2])
  for(i in seq_along(rclist)) {
    rc <- rclist[[i]]
    ri <- x[[i]]
    w <- dim(ri)
    rci <- c(1, w[1], 1, w[2])
    if(rc[1]>1)    rci[1] <- rci[1] + buffer[2]
    if(rc[2]<z[1]) rci[2] <- rci[2] - buffer[2]
    if(rc[3]>1)    rci[3] <- rci[3] + buffer[1]
    if(rc[4]<z[2]) rci[4] <- rci[4] - buffer[1]

    #browser()
    r0[rc[1]:rc[2], rc[3]:rc[4]] <- ri[ rci[1]:rci[2], rci[3]:rci[4] ]
  }
  r0
}


