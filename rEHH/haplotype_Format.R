#Input files: .PHASED from SHAPEIT
#Output files: -hap.txt and -map.txt
#Recomendation: The input files may contain the population name at the beggining, followed by the chromosome number, and then the SHAPEIT extention. 
#Consider pop1 then the file name would be: pop1-chr1.PHASED.sample

#Install packages + Loading Wrapping
{list.of.packages <- c("dplyr", "rehh", "ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)}

# Haplotype formating
#############################################################
popnames <- c("popname1","popname2","popname3")
for (popname in popnames) {
  for (i in seq(1,29)){
    haps <- read.table(paste0(popname,"-chr",i,".PHASED.haps"), stringsAsFactors = FALSE)
    sample <- read.table(paste0(popname,"-chr",i,".PHASED.sample"), stringsAsFactors = FALSE)
    #1.1) Count Frequencies 1
    
    # Count frequencies of haps (0,1) by row and convert to matrix
    count_frq <- apply(haps[,-(1:5)], 1, function(x) { table(factor(x, levels=0:1)) })
    tr_freq <- t(as.data.frame(count_frq))
    
    # Sum frequencies for each allele by row (SNP)
    al_a <- tr_freq[,c(T,F)] / (tr_freq[,c(T,F)] + tr_freq[,c(F,T)])
    al_b <- tr_freq[,c(F,T)] / (tr_freq[,c(T,F)] + tr_freq[,c(F,T)])

    freqsMatrix <- t(apply(data.frame(al_a, al_b), 1, function(row) 1 * (row == max(row))))
    freqsMatrix[freqsMatrix == 1] = 1
    freqsMatrix[freqsMatrix == 0] = 2
    
    #1.2) Write MAP 1
    
    # Major Allele in Control population will be "ancestral"
    maps <- haps[c(2,1,3)];
    finalMap <- as.matrix(cbind(maps, freqsMatrix))
    
    mapFileName <- paste0("chr-", i,"-",popname, "-map.txt")
    
    # write to file
    write.table(
      finalMap, 
      file = mapFileName, 
      quote = FALSE, 
      col.names = FALSE, 
      row.names = FALSE)
    
    
    #1.3) HAP Formatting 1
    
    # Remove header row from samples file
    sample <- sample[-(1:2),]
    
    # Convert haplotype data frame into matrix
    haps <- haps[,-(1:5)] %>% as.matrix
    
    # Assign column/row names as positions/individual IDs
    row.names(haps) <- haps[,3]
    colnames(haps) <- rep(sample$V2, each = 2)
    
    # Add 1 to each cell
    haps <- haps + 1
    
    # Transpose
    haps <- t(haps)
    
    # Annotate rows
    haplo.names <- rep(sample$V2, each = 2)
    first.haplos <- seq(from = 1, to = length(haplo.names), by = 2)
    second.haplos <- seq(from = 2, to = length(haplo.names), by = 2)
    
    haplo.reps <- c()
    haplo.reps[first.haplos] <- paste0(haplo.names[first.haplos],"_1")
    haplo.reps[second.haplos] <- paste0(haplo.names[second.haplos],"_2")
    
    haps <- data.frame(haps) %>%
      mutate(ind = haplo.reps) %>%
      select(ind, everything())

    tableHap <- paste0("chr-", i , "-",popname,"-hap.txt") 
    write.table(haps, file = tableHap, quote = FALSE, col.names = FALSE, row.names = FALSE)
    
    #1.4 Scan_hh() funtion 
    
    hap <- data2haplohh(
      hap_file = tableHap, 
      map_file = mapFileName, 
      haplotype.in.columns = F , 
      chr.name = i)
    ehhScan1 <- scan_hh(hap)
    scanFileName <- paste0("chr-", i ,"-",popname,"-scanhh.txt") 
    write.table(
      ehhScan1,
      file = scanFileName)
  }
}
