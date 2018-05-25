README
================

### About

The **ClustTools** package provides users with the evaluation tools cluster correlation, cluster index and SSE/SST (for r-squared computation). The `SSE` and `SST` functions are useful for r-squared computation, for the two popular hierarchical clustering functions `stats::hclust` and `factoextra::hcut` as these two functions don't contain SSE/SST computation. Moreover, the package offer the Hopkins statistic useful for cluster tendency computation.

### Practicalities

The package does not have any reverse dependencies outside the system library of R.

### Installation

To install latest version in R use following commands:

``` r
#install.packages('devtools')
#devtools::install_github('jhmadsen/ClustTools')
```

Work is currently carried out to make it available in the CRAN repository
