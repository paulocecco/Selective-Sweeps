#Input files: .PHASED from SHAPEIT
#Output files: -hap.txt and -map.txt

#Install packages + Loading Wrapping
{list.of.packages <- c("dplyr", "rehh", "ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)}
#source("../rehh_patches.R")

#############################################################
#
# Rsb CALCULATIONS
#
#############################################################

#Name the populations after the name files
popname1 <- "HoLo"
popname2 <- "LoHo"

{for (i in seq(1,29)){
  haps1 <- read.table(paste0(popname1,"-chr",i,".PHASED.haps"), stringsAsFactors = FALSE)
  sample1 <- read.table(paste0(popname1,"-chr",i,".PHASED.sample"), stringsAsFactors = FALSE)
  haps2 <- read.table(paste0(popname2,"-chr",i,".PHASED.haps"), stringsAsFactors = FALSE)
  sample2 <- read.table(paste0(popname2,"-chr",i,".PHASED.sample"), stringsAsFactors = FALSE)
  
############################################
# 1) POPULATION I CALCULATIONS
############################################

#1.1) Count Frequencies 1
  
  # Count frequencies of haps1 (0,1) by row and convert to matrix
  #View(haps1[,-(1:5)])
  count_frq <- apply(haps1[,-(1:5)], 1, function(x) { table(factor(x, levels=0:1)) })
  # Column 1 = 0
  # Column 2 = 1
  tr_freq <- t(as.data.frame(count_frq))
  
  # Sum frequencies for each allele by row (SNP)
  al_a <- tr_freq[,c(T,F)] / (tr_freq[,c(T,F)] + tr_freq[,c(F,T)])
  al_b <- tr_freq[,c(F,T)] / (tr_freq[,c(T,F)] + tr_freq[,c(F,T)])
  #View(data.frame(al_a,al_b))
  # Assign 1s to major frequency for each row
  freqsMatrix <- t(apply(data.frame(al_a, al_b), 1, function(row) 1 * (row == max(row))))
  # Check number of snps (rows)
  #dim(freqsMatrix)
  
  # all 1s become 2s and 0s become 1s
  freqsMatrix[freqsMatrix == 1] = 1
  freqsMatrix[freqsMatrix == 0] = 2
  #dim(freqsMatrix)
  #View(freqsMatrix)
  
#1.2) Write MAP 1
  
  # Major Allele in Control population will be "ancestral"
  maps1 <- haps1[c(2,1,3)];
  finalMap1 <- as.matrix(cbind(maps1, freqsMatrix))
  # haps1[which(FR < 0.5),]
  #dim(finalMap)
  
  mapFileName1 <- paste0("chr-", i,"-",popname1, "-map.txt")
  
  # write to file
  write.table(
    finalMap1, 
    file = mapFileName1, 
    quote = FALSE, 
    col.names = FALSE, 
    row.names = FALSE)
  
  
#1.3) HAP Formatting 1
  
  # Remove header row from samples file
  sample1 <- sample1[-(1:2),]
  
  # Convert haplotype data frame into matrix
  haps1 <- haps1[,-(1:5)] %>% as.matrix
  
  # Assign column/row names as positions/individual IDs
  row.names(haps1) <- haps1[,3]
  colnames(haps1) <- rep(sample1$V2, each = 2)
  
  # Add 1 to each cell
  haps1 <- haps1 + 1
  
  # Transpose
  haps1 <- t(haps1)
  
  # Annotate rows
  haplo.names <- rep(sample1$V2, each = 2)
  first.haplos <- seq(from = 1, to = length(haplo.names), by = 2)
  second.haplos <- seq(from = 2, to = length(haplo.names), by = 2)
  
  haplo.reps <- c()
  haplo.reps[first.haplos] <- paste0(haplo.names[first.haplos],"_1")
  haplo.reps[second.haplos] <- paste0(haplo.names[second.haplos],"_2")
  
  haps1 <- data.frame(haps1) %>%
    mutate(ind = haplo.reps) %>%
    select(ind, everything())
  # write to file (necessary for rehh for some reason)
  
  tableHap1 <- paste0("chr-", i , "-",popname1,"-hap.txt") 
  
  write.table(haps1, file = tableHap1, quote = FALSE, col.names = FALSE, row.names = FALSE)

#1.4 Scan_hh() funtion 1
  
  # Haplotype are recoded (if recode.allele option is activated) according to the ancestral and derived allele 
  # definition available in the map file (fourth and fifth columns) as :
  # 0=missing data, 
  # 1=ancestral allele, 
  # 2=derived allele. 
  
  hap1 <- data2haplohh(
    hap_file = tableHap1, 
    map_file = mapFileName1, 
    haplotype.in.columns = F , 
    # recode.allele = TRUE, 
    # min_maf = 0.1,
    # min_perc_geno.hap = 95,
    # min_perc_geno.snp = 50,
    chr.name = i)
  ehhScan1 <- scan_hh(hap1)
  scanFileName1 <- paste0("chr-", i ,"-",popname1,"-scanhh.txt") 
  write.table(
    ehhScan1,
    file = scanFileName1)
  
#1.5 iHS 1

  #ihsFileName1 <- paste0("chr-", i , "-rehh_ihs_pop1.txt")  
  #ihs1 <- ihh2ihs(ehhScan1, freqbin=0.5)
  #write.table(ihs1, file = ihsFileName1)
  #write.csv(ihs1, file = "rehh_table_ihs_pop-1.csv")
  
  #ihsFileName1Png <- paste0("chr-", i , "-rehh_ihs_pop-1.png") 
  #ihsTitleName1 <- paste0("iHS Population 1 chr " , i)
  #png(filename = ihsFileName1Png)
  #ihsplot_mod(ihs1, plot.pval = TRUE, ylim.scan=2, main = ihsTitleName1)
  #dev.off()

############################################
# 2) POPULATION II CALCULATIONS
############################################
  
#2.1) Count Frequencies 2

  # Count frequencies of haps2 (0,1) by row and convert to matrix
  #View(haps2[,-(1:5)])
  count_frq <- apply(haps2[,-(1:5)], 1, function(x) { table(factor(x, levels=0:1)) })
  # Column 1 = 0
  # Column 2 = 1
  tr_freq <- t(as.data.frame(count_frq))
  
  # Sum frequencies for each allele by row (SNP)
  al_a <- tr_freq[,c(T,F)] / (tr_freq[,c(T,F)] + tr_freq[,c(F,T)])
  al_b <- tr_freq[,c(F,T)] / (tr_freq[,c(T,F)] + tr_freq[,c(F,T)])
  #View(data.frame(al_a,al_b))
  # Assign 1s to major frequency for each row
  freqsMatrix <- t(apply(data.frame(al_a, al_b), 1, function(row) 1 * (row == max(row))))
  # Check number of snps (rows)
  #dim(freqsMatrix)
  
  # all 1s become 2s and 0s become 1s
  freqsMatrix[freqsMatrix == 1] = 1
  freqsMatrix[freqsMatrix == 0] = 2
  #dim(freqsMatrix)
  #View(freqsMatrix)
  
#2.2 Write MAP 2
  
  # Major Allele in Control population will be "ancestral"
  maps2 <- haps2[c(2,1,3)];
  finalMap2 <- as.matrix(cbind(maps2, freqsMatrix))
  # haps2[which(FR < 0.5),]
  #dim(finalMap)
  
  
  mapFileName2 <- paste0("chr-", i ,"-",popname2,"-map.txt")
  
  # write to file
  write.table(
    finalMap2, 
    file = mapFileName2, 
    quote = FALSE, 
    col.names = FALSE, 
    row.names = FALSE)
  
  
#2.3 HAP Formatting 2
  
  # Remove header row from samples file
  sample2 <- sample2[-(1:2),]
  
  # Convert haplotype data frame into matrix
  haps2 <- haps2[,-(1:5)] %>% as.matrix
  
  # Assign column/row names as positions/individual IDs
  row.names(haps2) <- haps2[,3]
  colnames(haps2) <- rep(sample2$ID_1, each = 2)
  
  # Add 1 to each cell
  haps2 <- haps2 + 1
  
  # Transpose
  haps2 <- t(haps2)
  
  # Annotate rows
  haplo.names <- rep(sample2$V2, each = 2)
  first.haplos <- seq(from = 1, to = length(haplo.names), by = 2)
  second.haplos <- seq(from = 2, to = length(haplo.names), by = 2)
  
  haplo.reps <- c()
  haplo.reps[first.haplos] <- paste0(haplo.names[first.haplos],"_1")
  haplo.reps[second.haplos] <- paste0(haplo.names[second.haplos],"_2")
  
  haps2 <- data.frame(haps2) %>%
    mutate(ind = haplo.reps) %>%
    select(ind, everything())
  # write to file (necessary for rehh for some reason)
  
  tableHap2 <- paste0("chr-", i, "-",popname2,"-hap.txt") 
  
  write.table(haps2, file = tableHap2, quote = FALSE, col.names = FALSE, row.names = FALSE)
  

#2.4) Scan_hh() funtion 2
  
  # Haplotype are recoded (if recode.allele option is activated) according to the ancestral and derived allele 
  # definition available in the map file (fourth and fifth columns) as :
  # 0=missing data, 
  # 1=ancestral allele, 
  # 2=derived allele. 
  
  hap2 <- data2haplohh(
    hap_file = tableHap2, 
    map_file = mapFileName1, 
    haplotype.in.columns = F , 
    # recode.allele = TRUE, 
    # min_maf = 0.1,
    # min_perc_geno.hap = 95,
    # min_perc_geno.snp = 50,
    chr.name = i)
  ehhScan2 <- scan_hh(hap2)
  scanFileName2 <- paste0("chr-", i , "-",popname2,"-scanhh.txt") 
  write.table(
    ehhScan2,
    file = scanFileName2)
}
