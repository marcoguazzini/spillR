# spillR

 # Algorithm for real data.
 
How to use our package on real data:
 ```r
 # install and load pkg
devtools::install_github("marcoguazzini/spillR")
library("spillR")

# Libraries
library(spillR)
library(flowCore)
library(ggplot2)
library(tibble)
library(dplyr)
library(RColorBrewer)
library(magrittr)
library(tidyr)
library(CATALYST)
```

```{r}
 # load counts data from CATALYST pkg (this is an example of counts data, one can use an 
 available dataset of mass cytometry data)
 counts_real <- spillR::load_real_data()
 
 # load counts from the spillover experiments
 beads_exp <- spillR:: load_beads_data()
 counts_bead <- beads_exp[[1]]
 sm <- beads_exp[[2]]
 
 ```

```r
# Compensation
target_marker <- rownames(sm)
fit<-
    lapply(target_marker, function(marker)
        EM_mixture(
          marker,
          counts_nb,
          counts_bead,
          n_iter = 20,
          sm
        )
      )
  counts_comp <- fit[[1]]
```

```r
# visualization of the results
# choosing a marker
marker <- "Nd145Di"
plot_compensation(fit[[counts_comp[,marker], 
                              marker 
                              )

```


 
 
 
 
 
 
```r 
# DAG representing the spillover
# set this for better image resolution fig.width=10, fig.height=7, out.width="100%"
spillR::spill_dag(spillR::load_spillover())
```





