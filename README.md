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

## Guide
- [Plink](#plink)
- [ShapeIT](#shapeit)
- [OUTFLANK](#outflank)
- [rEHH](#rehh)
- [Manhattan Plot](#manhattan-plot)

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
```
#95% of genotype and a minimum allele frequency of 1% are kept
plink2 --bfile BedFile --out FilterFile --geno 0.05 --maf 0.01 --make-bed --cow 
```

#### 3. Chromosome splitting
```
#Create a file for each Chromosome
for i in $(seq 1 29); do plink2 --bfile FilterFile --out FinalFile-BTA$i --chr $i --cow --make-bed; done
```

### ShapeIT
SHAPEIT will require a bed file extention as input. To run the following script you'll need to create a directory where you have all the bed files. In this same directory it's recomended to have the script saved. SHAPEIT will create a new directory which is going to contain all the **PHASED** output files.

See [runshapeIT.sh](https://github.com/paulocecco/Selective-Sweeps/tree/main/SHAPEIT)

### OUTFLANK
OutFLANK is used to estimate FST with efective size population correction. In contrast to the classical REYNOLDS, WEIR, and COCKERHAM FST estimation which has the general assumption that all population has the same effective size, this could lead to miss estimations such as outliers or even negative values of the estimator. OUTFLANK uses population size correction to estimate FST in a more acurate and precised way.

To run the script you must have a **vcf** file obtained from plink. The output it's going to be a csv file with the SNPs obtained and 4 columns named: CHR (Chromosome), POSITION, SNP (marker name), FST (weir Cockerham estimation), FSTNoCorr (OUTFLANK correction) and pvaluesRightTail (significant SNPs p < 0.01)
For script see 

### rEHH
Extended Homozygous Haplotype (EHH) is an estimation of the area under the curve of the haplotype extention. In this repository only Rsb and XP-EHH are listed. The differences between them are subtle. While Rsb is used to reveal recent haplotype formation between populations, the cross population (XP) can reveal more ancient selective sweeps. You can see all documentation in the paper posted by Gautier, et. al 2017 (https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12634) and his repository at https://cran.r-project.org/web/packages/rehh/index.html for more information.

For this script you'll need **.PHASED.sample** and **.PHASED.haps** files from *SHAPEIT* as input. A csv file with the top 1% of the markers, is going to be created. The script is divided into two different parts. First, the [haplotype_format](https://github.com/paulocecco/Selective-Sweeps/blob/main/rEHH/haplotype_Format.R) script that it's going to create an **hap.txt** and a **scanhh.txt** files. This files are going to be the input for the calculations listed in [rEHH calculations](https://github.com/paulocecco/Selective-Sweeps/blob/main/rEHH/rEHH_calc.R).

### Manhattan Plot


## Contributing
We welcome contributions to this repository. If you find any bugs or have suggestions for improvements, please open an issue or submit a pull request.

## License
The code in this repository is available under the MIT License.
