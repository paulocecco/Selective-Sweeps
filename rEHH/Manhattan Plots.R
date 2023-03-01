list.of.packages <- c("dplyr", "rehh", "ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

#############################################################
#
# Scan HH
#
#############################################################

for(i in 1:29) {
  hapFile1 <- paste0("chr-", i ,"-rehh_hap_pop-3.txt")
  hapFile2 <- paste0("chr-", i ,"-rehh_hap_pop-2.txt")
  mapFile <- paste0("chr-", i ,"-rehh_map_pop-3.txt")
  data1 <- data2haplohh(
    hap_file = hapFile1, 
    mapFile, 
    chr.name=i,
    min_perc_geno.hap = 50,
    allele_coding = "12")
  
  data2 <- data2haplohh(
    hap_file = hapFile2, 
    mapFile, 
    chr.name=i, 
    allele_coding = "12",
    min_perc_geno.hap = 50)
  
  ehhScan1 <- scan_hh(data1)
  ehhScan2 <- scan_hh(data2)  
  if(i==1) {
    wg.res1 <- ehhScan1
    wg.res2 <- ehhScan2 }  
  else { 
    wg.res1 <- rbind(wg.res1, ehhScan1)
    wg.res2 <- rbind(wg.res2, ehhScan2)
  }
}


#############################################################
#
# FORMATTING
#
#############################################################

# Set rowname to first column name
wg.res1rm <- tibble::rownames_to_column(wg.res1,"SNPs")
wg.res2rm <- tibble::rownames_to_column(wg.res2,"SNPs")
# Difference between res1 and res2
wg.res1filtered <- wg.res1rm[wg.res1rm$SNPs %in% wg.res2rm$SNPs,]
# Clean row names
rownames(wg.res1filtered) <- c()
wg.res1rm <- tibble::column_to_rownames(wg.res1filtered,"SNPs")

#############################################################
#
# RSB
#
#############################################################

#Calculation
wg.rsb <- ines2rsb(wg.res1rm, wg.res2)
wg.rsb <- wg.rsb[which(wg.rsb$LOGPVALUE>0) ,]
wg.rsb$pvalue <- 10^-wg.rsb$LOGPVALUE

#Collect and save top 1%
RSB_oneperc <- round(length(wg.rsb$LOGPVALUE)/100)
RSB_newdata <- wg.rsb[order(-wg.rsb$LOGPVALUE),]
RSB_newdata <- RSB_newdata[1:RSB_oneperc,]

#Collect p <0.01
RSB_sigp <- na.omitwg.rsb[wg.xpehh$pvalue<0.01,])

#############################################################
#
# XPEHH
#
#############################################################

wg.xpehh <- ies2xpehh(wg.res1rm, wg.res2)
wg.xpehh$pvalue <- 10^-wg.xpehh$LOGPVALUE

#Collect and save top 1%
xpehh_oneperc <- round(length(wg.xpehh$LOGPVALUE)/100)
xpehh_newdata <- wg.xpehh[order(-wg.xpehh$LOGPVALUE),]
xpehh_newdata <- xpehh_newdata[1:xpehh_oneperc,]

#Collect p <0.01
xpehh_sigp <- na.omit(wg.xpehh[wg.xpehh$pvalue<0.01,])

#############################################################
#
# Output
#
#############################################################

popnames <- "popname1v2"
#For top 1%
write.csv(RSB_newdata, paste0("RSB_",popnames,"_OnePerc.csv"), quote = F)
write.csv(xpehh_newdata, paste0("XPEHH_",popnames,"_OnePerc.csv"), quote = F)
#For p < 0.01
write.csv(RSB_sigp, paste0("RSB_",popnames,"_signif.csv"), quote = F)
write.csv(xpehh_sigp, paste0("XPEHH_",popnames,"_signif.csv"), quote = F)
