# spillR

https://academic.oup.com/bioinformatics/article/40/6/btae337/7689679

## Description

Channel interference in mass cytometry can cause spillover and may result in miscounting of protein markers. We develop a nonparametric finite mixture model, and use the mixture components to estimate the probability of spillover. We implement our method using expectation-maximization to fit the mixture model.

## Installation

Our R package is available on [Bioconductor](https://bioconductor.org/packages/spillR):

```	r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("spillR")
```
