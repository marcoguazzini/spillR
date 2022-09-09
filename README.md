# spillR

How to use our package:

```r
# install and load pkg
devtools::install_github("marcoguazzini/spillR")
library("spillR")

# parameters of true expressions
shape <- 9.0
rate <- 1.0

# load data from CATALYST pkg
sm <- spillR::load_spillover()

# generate counts with spillover
counts <- spillR::prepare_data(rate, shape, sm)

# add spillover
counts_compensated <- spillR::compensate(counts)

# plot comparison
spillR::plot(counts[,1], counts_compensated[,1])
```
 
