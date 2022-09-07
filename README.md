# SpillR

How to use our package:

```r
# install and load pkg
devtools::install_github("marcoguazzini/SpillR")
library("SpillR")

# parameters of true expressions
shape <- 9.0
rate <- 1.0

# load data from CATALYST pkg
sm <- SpillR::load_spillover()

# generate counts with spillover
counts <- SpillR::prepare_data(rate, shape, sm)

# add spillover
counts_compensated <- SpillR::compensate(counts)

# plot comparison
SpillR::plot_density(counts[,1], counts_compensated[,1])
```
 
