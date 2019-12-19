<img src="https://github.com/geoportti/Logos/blob/master/geoportti_logo_300px.png">

# R-package: Contextual Mann-Kendall trend test for rasterStacks

Compute a Mann-Kendall trend-test in each cell of a given `rasterStack` where layers are considered as observation times.

## Main features

The package is designed to work alongside the package `raster`. 

For a given `rasterStack`-object (or `rasterBrick`) with >3 layers representing observation times, the package  computes quickly the the nearest-neighbour smoothed aka "contextual" Mann-Kendall test (`contextual_mann_kendall`; basic cell-wise Mann-Kendall test included as a special case).

Additional "wrappers" for useful operations on the `rasterStack` or result `raster` include

* `snow` parallelised computation for large `rasterStacks` (`split_calc_wrapper`)
* p-value adjustment for multiple testing
* Pettitt's Test For a Change-Point
* Cox-Stuart trend test
* Wan-Swail prewhitening assuming AR1-noise

Depends on the packages `raster` , `sp` and ` progress`. 

## Installation

As usual,

```
library(devtools)
install_github("antiphon/ConMK", build_vignettes = TRUE)
```

to include vignettes.

## Usage

The packages comes with a set of synthetic `rasterStack` objects:
```
> data("test_stacks2")
> x <- test_stacks2$trend
> r <- contextual_mann_kendall(x)
> r
class      : RasterStack 
dimensions : 40, 60, 2400, 3  (nrow, ncol, ncell, nlayers)
resolution : 0.01666667, 0.01666667  (x, y)
extent     : 0, 1, 0, 0.6666667  (xmin, xmax, ymin, ymax)
crs        : +init=epsg:3067 +proj=utm +zone=35 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs 
names      :             S,            s2,             p 
min values : -3.750000e+01,  3.302724e+01,  3.255194e-06 
max values :      108.4444,      545.6284,        1.0000 

> plot(r)
> plot( p.adjust_raster(r$p) < 0.05 )
```

See the vignette for further examples,

```
vignette(package="ConMK")
```

## Usage and Citing
When using GeoCubes, the following citing should be mentioned:
"We made use of geospatial data/instructions/computing resources provided by the Open Geospatial Information Infrastructure for Research (oGIIR, urn:nbn:fi:research-infras-2016072513) funded by the Academy of Finland."
