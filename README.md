# 2019 DRC Ebola analysis  
> Sean Browning  

## What?  

Random visualizations and analysis work on the current DRC Ebola outbreak.


![](output/viz.gif)

## Dependencies  

All data pulled from HDX using [rhdx](https://github.com/dickoa/rhdx)  
- [DRC case data](https://data.humdata.org/dataset/ebola-cases-and-deaths-drc-north-kivu)
- [DRC shape file](https://data.humdata.org/dataset/democratic-republic-of-congo-health-boundaries) 
- [DRC Zone de Sante denominator data](https://data.humdata.org/dataset/dr-congo-health-0)  

All packages for vis and modelling are handled via [packrat](https://github.com/rstudio/packrat)  

- Running R in project root folder should initialize packrat install  
- run `packrat::restore()` to install all packages to the local directory  

## Updating case counts  

Updating the CSVs and shape file is fully automated through [R/update_ebola_counts.R](https://github.com/beansrowning/drc_ebola_2019/blob/master/R/update_ebola_counts.R)  

```shell
Rscript R/update_ebola_counts.R
```