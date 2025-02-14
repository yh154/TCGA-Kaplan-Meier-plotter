# TCGA-Kaplan-Meier-plotter
A Shiny app for plotting TCGA Kaplan-Meier curve of a given gene 

``` r
shiny::runApp("TCGA-KM.R")
```
Packages involved: 
``` r
suppressPackageStartupMessages({
require(shiny)
require(TCGAbiolinks)
require(survival)
require(survminer)
require(dplyr)
require(DT)
})
```
TCGA projects:
``` r
TCGA-ACC, TCGA-BLCA, TCGA-BRCA, TCGA-CESC, TCGA-CHOL, TCGA-COAD, TCGA-DLBC, TCGA-ESCA, TCGA-GBM, TCGA-HNSC, TCGA-KICH, TCGA-KIRC,
TCGA-KIRP, TCGA-LAML, TCGA-LGG, TCGA-LIHC, TCGA-LUAD, TCGA-LUSC, TCGA-MESO, TCGA-OV, TCGA-PAAD, TCGA-PCPG, TCGA-PRAD, TCGA-READ,
TCGA-SARC, TCGA-SKCM, TCGA-STAD, TCGA-TGCT, TCGA-THCA, TCGA-THYM, TCGA-UCEC, TCGA-UCS, TCGA-UVM

![alt text]([https://raw.githubusercontent.com/yh154/TCGA-Kaplan-Meier-plotter/master/layout.png](https://raw.githubusercontent.com/yh154/TCGA-Kaplan-Meier-plotter/main/layout.png)
