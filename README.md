# spillR

 # Algorithm for real data.
 
How to use our package on real data:
 ```{r}
 # install and load pkg
devtools::install_github("marcoguazzini/spillR")
library("spillR")

 # load counts data from CATALYST pkg (this is an example of counts data, one can use an available dataset of mass cytometry data)
 counts <- spillR::extract_data()
 
 # load counts from the spillover experiments
 counts_spill <- spillR:: extract_spill_distr()
 
 #load spillover matrix
 sm <- spillR::load_spillover()
 ```

```{r}

```







# Simulation with method presented on EuroBioc2022

How to use our package:

```r

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
counts_compensated <- spillR::compensate(counts, sm)

# plot comparison
spillR::plot(counts, counts_compensated, channel_names)
```
 

 
 
 
 
 
 
 
 
```r 
# DAG representing the spillover
# set this for better image resolution fig.width=10, fig.height=7, out.width="100%"
spillR::spill_dag(spillover_matrix)
```





