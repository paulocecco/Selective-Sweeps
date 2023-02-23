#############################################################
#
# FST
#
#############################################################
setwd("C:/Users/pecun/OneDrive/Desktop/Brangus Paper/data/Raw/FST")

fst <- read.table("ZeBr.fst", header=T)
fst <- fst[,c(1,2,3,5)]
fst <- fst[which(fst$FST>0),]
colnames(fst) <- c("CHR", "SNP","POSITION", "LOGPVALUE")

#Collect and save top 1%
fst_oneperc <- round(length(fst$LOGPVALUE)/100)
fst_newdata <- fst[order(-fst$LOGPVALUE),]
fst_newdata <- fst_newdata[1:fst_oneperc,]

fst_top_val <- fst_newdata[fst_oneperc,]
fst_top_val <- fst_top_val$LOGPVALUE
fst_y_topval <- round(max(fst_newdata$LOGPVALUE))
fst_y_minval <- round(min(na.omit(fst$LOGPVALUE)))
fst_thresh_val <- min(fst_newdata$LOGPVALUE)


#############################################################
#
# SCAN_HH
#
#############################################################
setwd("C:/Users/pecun/OneDrive/Desktop/Brangus Paper/data/Raw/LD")

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

#Collect and save top 1%
RSB_oneperc <- round(length(wg.rsb$LOGPVALUE)/100)
RSB_newdata <- wg.rsb[order(-wg.rsb$LOGPVALUE),]
RSB_newdata <- RSB_newdata[1:RSB_oneperc,]

RSB_top_val <- RSB_newdata[RSB_oneperc,]
RSB_top_val <- RSB_top_val$LOGPVALUE
RSB_y_topval <- round(max(RSB_newdata$LOGPVALUE))
RSB_y_minval <- round(min(na.omit(wg.rsb$LOGPVALUE)))
RSB_thresh_val <- round(min(RSB_newdata$LOGPVALUE),digit=1)



#############################################################
#
# XPEHH
#
#############################################################

wg.xpehh <- ies2xpehh(wg.res1rm, wg.res2)

#Collect and save top 1%
xpehh_oneperc <- round(length(wg.xpehh$LOGPVALUE)/100)
xpehh_newdata <- wg.xpehh[order(-wg.xpehh$LOGPVALUE),]
xpehh_newdata <- xpehh_newdata[1:xpehh_oneperc,]

XPEHH_top_val <- xpehh_newdata[xpehh_oneperc,]
XPEHH_top_val <- XPEHH_top_val$LOGPVALUE
XPEHH_y_topval <- round(max(xpehh_newdata$LOGPVALUE))
XPEHH_y_minval <- round(min(na.omit(wg.xpehh$LOGPVALUE)))
XPEHH_thresh_val <- round(min(xpehh_newdata$LOGPVALUE), digits = 1)

#############################################################
#
# GRAPH
#
#############################################################

#Zebuine #00BFC4
#Angus #F8766D
#Brangus #7CAE00

palette(c("#F8766D","#7CAE00"))

#Manhattanplot FST####
#png(filename = paste0("FST-", comparison, ".png"), width = 1240, height = 720)
maplot(fst, pval = T, threshold = fst_thresh_val,
       cex=0.7)
#dev.off()
#svg(filename = paste0("FST-", comparison, ".svg"), width = 1240, height = 720)
manplot(fst, threshold = fst_thresh_val)
#dev.off()


#Manhattanplot RSB####
#png(filename = paste0("Rsb-", comparison, ".png"), width = 1240, height = 720)
manhattanplot(wg.rsb, pval = T, ylim = c(RSB_y_minval, RSB_y_topval), 
              threshold = RSB_thresh_val, color = palette, cex = 0.7)
#dev.off()
#####
#svg(filename = paste0("Rsb-", comparison, ".svg"), width = 1240, height = 720)
#manhattanplot(wg.rsb, pval = T, ylim = c(RSB_y_minval, RSB_y_topval),threshold = RSB_thresh_val)
#dev.off()

#Manhattanplot XPEHH#####
#png(filename = paste0("XPEHH-", comparison, ".png"), width = 1240, height = 720)
manhattanplot(wg.xpehh, pval = T, ylim = c(XPEHH_y_minval, XPEHH_y_topval), 
              threshold = XPEHH_thresh_val,color = palette,
              cex=0.7)

#dev.off()
#####
#svg(filename = paste0("XPEHH-", comparison, ".svg"), width = 1240, height = 720)
#manhattanplot(wg.xpehh, pval = T, ylim = c(XPEHH_y_minval, XPEHH_y_topval), threshold = XPEHH_thresh_val)              
#dev.off()




#SAVING FILES####
write.csv(fst_newdata, paste0("FST " ,comparison, " complete.csv"), quote = F)
write.csv(wg.rsb,paste0("Rsb " ,comparison, " complete.csv"), quote = F, row.names = T)
write.csv(wg.xpehh,paste0("XPEHH " ,comparison, " complete.csv"), quote = F, row.names = T)
#write.csv(RSB_newdata, "Rsb One_Percent LBC-Ho.csv", quote = F)
#write.csv(xpehh_newdata, "XPEHH One_Percent LBC-Ho.csv", quote = F)
