
# groebner

<!-- badges: start -->
<!-- badges: end -->

The groebner package provides a pure R implementation of groebner bases for finding the roots of 
polynomial systems

## Installation

The current version of groebner can be installed from GitHub using the remotes package. 
```r
# install.packages("remotes")
remotes::install_github("SWotherspoon/groebner")
```


## Example

To find the intersection of the unit circle with a hyperbola

``` r
library(groebner)
## Define system
ps <- parse_polys(c("x^2+y^2-1",
                    "x*y-1/4"),
                  c("x","y"))
## Find roots
rs <- solve_polys(ps)
rs
## Check
rowSums(abs(eval_polys(ps,rs)))
```

