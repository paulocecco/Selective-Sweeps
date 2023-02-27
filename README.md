# Selective-Sweeps
## Estimation of selective sweeps using different methods in R

A selective sweep is a rapid evolution process in which a beneficial mutation spreads throughout a population. This process can have a major impact on the genetic diversity of a species and can be used to study evolution and adaptation.

This repository contains scripts and tools for performing selective sweep analysis using R. The goal is to provide an easy-to-use and efficient way to detect selective sweeps in genomic data.

All data in this repository has been published in https://onlinelibrary.wiley.com/doi/abs/10.1111/jbg.12733

# Index
- [Installation](#installation)
- [Usage](#usage)
- [Methods](#methods)
- [Guide](#guide)
- [Contributing](#contributing)
- [License](#license)

## Installation
You'll need to install SHAPEIT software (https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#citations). The output of this program **.PHASED** files are going to be the input for the following workflow. Please also read OutFLANK repository (https://github.com/whitlock/OutFLANK) since it's the one we are going to use.

The scripts in this repository require the following R packages:
* dplyr
* rehh
* ggplot2
* qvalue
* vcfR
* OutFLANK
* bigsnpr
* qqman
* stringr

For quick installation
```{r}
#Install packages + Loading Wrapping
{list.of.packages <- c("dplyr", "rehh", "ggplot2", "qvalue", "OutFLANK", "vcfR", "bigsnpr", "qqman", "stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)}
```

## Usage
To use the scripts, you need to provide the following inputs:

* A **vcf** file containing variant information for multiple individuals
* A **.PHASED** file which is the output from **SHAPEIT**
* One file per each chromosome
* A **.map** file obtained from **plink** (https://www.cog-genomics.org/plink/)

Here's a quick example of the workflow
![Example](https://user-images.githubusercontent.com/43005715/218570525-d0a224c9-4a77-4ea6-8027-22d14b61c884.png)


## Methods
In this script we are going to use two different approches for selective sweeps estimations.
1. Population differentiation: **FST**
2. Linkage of desequilibrum: estimation of extent haplotype homozygous **Rsb** and **XPEHH**

For the first one we are going to use the **Outflank** package while **rehh** is going to be used for the second one.

## Step-by-step Guide
- [Plink](#plink)
- [ShapeIT](#shapeit)
- [OUTFLANK](#outflank)
- [rEHH](#rehh)

### Plink
plink it's used for 4 things
1. Creating **.bed** **.vcf** extention files needed for ShapeIT and OUTFLANK softwares.
2. Split population by family tag using **.fam** file.
3. Filtering data by genotype (**--geno**) and minimum allele frequency (**--maf**), required for rEHH.
4. Chromosome splitting, needed for ShapeIT.

#### 1. Extention
Using **.bed** file extention is recommended since **.fam** file is created 
```
#From .ped to .bed
#--cow flag is used since data is from cattle data
plink2 --file InputFile --Out BedFile --make-bed --cow 
```
Afterwards **.vcf** file is created
```
#Using .bed to create .vcf
plink2 --bfile BedFile --out VCFFile --recode vcf --cow
```

#### 2. Filtering by genotype and MAF

#### 3. Chromosome splitting

### ShapeIT

### OUTFLANK

### rEHH

## Contributing
We welcome contributions to this repository. If you find any bugs or have suggestions for improvements, please open an issue or submit a pull request.

## License
The code in this repository is available under the MIT License.
