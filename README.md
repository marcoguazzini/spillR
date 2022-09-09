# spillR

How to use our package:

```r
# install and load pkg
devtools::install_github("marcoguazzini/spillR")
library("spillR")

# parameters of true expressions
shape <- 9.0
rate <- 1.0
n_cells <- 1000

# channels for comparison 
channel_names <- c("Nd143Di","Nd148Di")

# load data from CATALYST pkg
sm <- spillR::load_spillover()

# generate counts with spillover
counts <- spillR::prepare_data(sm, shape, rate, n_cells)

# add spillover
counts_compensated <- spillR::compensate(counts, sm, n_cells)

# plot comparison
spillR::plot(counts, counts_compensated, channel_names)
```
 
