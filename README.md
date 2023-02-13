# Selective-Sweeps
##Estimation of selective sweeps using different methods in R

A selective sweep is a rapid evolution process in which a beneficial mutation spreads throughout a population. This process can have a major impact on the genetic diversity of a species and can be used to study evolution and adaptation.

This repository contains scripts and tools for performing selective sweep analysis using R. The goal is to provide an easy-to-use and efficient way to detect selective sweeps in genomic data.

##Installation
You'll need to install SHAPEIT software (https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#citations). The output of this program **.PHASED** files are going to be the input for the following workflow

The scripts in this repository require the following R packages:
* dplyr
* rehh
* ggplot2

For quick installation
```
#Install packages + Loading Wrapping
{list.of.packages <- c("dplyr", "rehh", "ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)}
```

##Usage
To use the scripts, you need to provide the following inputs:

* A **.PHASED** file which is the output from **SHAPEIT**
* One file per each chromosome
* A **.map** file obtained from **plink** (https://www.cog-genomics.org/plink/)
