# locus – large-scale variational inference for combined predictor and outcome selection

[![Travis-CI Build Status](https://travis-ci.org/hruffieux/locus.svg?branch=master)](https://travis-ci.org/hruffieux/locus)
 
## Overview

**locus** is an R package providing efficient variational algorithms for
simultaneous variable selection of predictors and associated outcomes based
on multivariate regression models. Dependence across outcomes linked to the 
same predictors is captured through the model hierarchical structure 
(Hélène Ruffieux, Anthony C. Davison, Jörg Hager, Irina Irincheeva, 2016, 
arXiv:1609.03400). 

## Installation

To install, run the following commands in R:

``` r
install.packages("devtools")
devtools::install_github("hruffieux/locus")
```

## License

This project uses the GPL v2 license, see [LICENSE](LICENSE).


## Issues

To report an issue, please use the [locus issue tracker](https://github.com/hruffieux/locus/issues) at github.com.
