{list.of.packages <- c("dplyr","gtools","qvalue", "vcfR", "OutFLANK", "bigsnpr","qqman","stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)}

##############################
## Upload Compressed VCF file
##############################

vcf.filename <- "HiLo"
obj.vcfR <- read.vcfR(paste0(vcf.filename,".vcf"), verbose=FALSE) 
#Extract Genotype, Positions, Chromosome and Locinames
geno <- extract.gt(obj.vcfR) # Character matrix containing the genotypes
position <- getPOS(obj.vcfR) # Positions in bp
chromosome <- as.integer(as.factor(getCHROM(obj.vcfR))) # Chromosome information
locinames <- obj.vcfR@fix[,"POS"]
locinames <- paste(obj.vcfR@fix[,"CHROM"], obj.vcfR@fix[,"POS"], obj.vcfR@fix[,"ID"],sep="_")
#Build-up the genotype matrix
G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))
G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2
#Check for NA. TOTAL NECESARY THAT NA THAT HAS TO BE CONVERTED INTO 9
if (sum(is.na(G)) != 0) {
  print(paste0("A total of: ", sum(is.na(G)), " NAs found in data set. Correction in progress."))
  G[is.na(G)] <- 9
  print(paste0("A total of: ", sum(is.na(G)), " NAs found in data set. Correction done")) 
} else { print ("No NAs in data set. Correction passed.")}
dim(G) == dim(geno) #Checkpoint

##############################
## Prunning using bigsnpr 
##############################
#Prun the data cause SNPs are not independent as FST estimations assumes due to LD. This will be used to estimate a more realistic mean of FST.
#A new matrix has to be made cause bigsnpr uses missing value as 3 while outFLANK recognizes as 9
G2 <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))
G2[geno %in% c("0/0", "0|0")] <- 0
G2[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G2[geno %in% c("1/1", "1|1")] <- 2
#if (sum(is.na(G2)) != 0) {
#  print(paste0("A total of: ", sum(is.na(G2)), " NAs found in data set. Correction in progress."))
#  G2[is.na(G)] <- 9
#  print(paste0("A total of: ", sum(is.na(G)), " NAs found in data set. Correction done")) 
#} else { print ("No NAs in data set. Correction passed.")}
#SNP trimming
G3 <- add_code256(big_copy(t(G2),type="raw"),code=bigsnpr:::CODE_012)
newpc <- snp_autoSVD(G=G3, infos.chr = chromosome, infos.pos = position,
                     size=100, roll.size = 40)
which_pruned <- attr(newpc, which="subset") # Indexes of remaining SNPS after pruning
length(which_pruned)
#plot(newpc)

##############################
### The individual data
##############################
ind <- read.table(paste0(vcf.filename,".fam"), 
                  header=F, 
                  col.names = c("pop", "ind", "father","mother","sex", "phenotype"))
pop <- as.character(ind$pop)

##############################
### FST Estimations
##############################
my_fst <- MakeDiploidFSTMat(SNPmat = t(G), locusNames = locinames, popNames = pop)
#sum(is.na(my_fst$FST))
#nrow(my_fst[my_fst$FST<0,])
#{NAerror <- paste0("There're a total of ", sum(is.na(my_fst$FST))," NAs values.")
#Negerror <-  paste0("There're a total of ", nrow(my_fst[my_fst$FST<0,])," negative values. ")
#TOTerror <- paste0("Total of errors for FST: ", sum(is.na(my_fst$FST)) + nrow(my_fst[my_fst$FST<0,]))
#FITdf <- paste0("Fitted dataframe length: ",
#                length(my_fst$LocusName) - (sum(is.na(my_fst$FST)) + nrow(my_fst[my_fst$FST<0,])),
#                                           " SNPs")
#cat(NAerror,"\n",Negerror,"\n",TOTerror,"\n",FITdf)}

#Heterozygosity vs. FST: Low He SNPs have high FST values.
#plot(my_fst$He, my_fst$FSTNoCorr)
#Lets plot the corrected vs. uncorrected FSTs: Here you find for loci that deviate for linear trend, as all SNPs were all genotyped in the same number of individuals, no deviation will be observed
#plot(my_fst$FST, my_fst$FSTNoCorr)
#abline(0,1)
#hist(my_fst$FSTNoCorr, breaks=seq(0,1, by=0.005))
#hist(my_fst$FSTNoCorr[my_fst$He>0.4], breaks=seq(0,1, by=0.05))

##############################
### OutFLANK
##############################

#str(out_trim)
#Without  prunned filter
#OutFLANKResultsPlotter(out_trim, withOutliers = TRUE,
#                       NoCorr = T, Hmin = 0.1, binwidth = 0.001, Zoom =
#                         T, RightZoomFraction = 0.01, titletext = NULL)
#hist(out_trim$results$pvaluesRightTail)
#With prunn filter
prunn <- OutFLANK(my_fst[which_pruned,], NumberOfSamples=2, 
                  qthreshold = 0.05, Hmin = 0.1)
#Prunned
prunn.chi1 <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = prunn$FSTNoCorrbar, 
                                dfInferred = prunn$dfInferred, 
                               qthreshold = 0.05, Hmin=0.1)

sum(prunn.chi1$OutlierFlag, na.rm = TRUE)
lociname.split.Prunn <- str_split_fixed(prunn.chi1$LocusName, "_", 3)
colnames(lociname.split.Prunn) <- c("CHR","POSITION","SNP")
prunn.chi1.f <- cbind(lociname.split.Prunn,prunn.chi1)
prunn.chi1.f <- prunn.chi1.f[,-4]
prunn.chi1.f$CHR <- as.integer(prunn.chi1.f$CHR)
prunn.chi1.f$POSITION <- as.integer(prunn.chi1.f$POSITION)
#Full Table
df.export.Prunn <- prunn.chi1.f[,c(1,2,3,5,8,12,13)]
df.export.Prunn <- na.omit(df.export.Prunn)
df.export.Prunn <- df.export.Prunn[order(-df.export.Prunn$FSTNoCorr, df.export.Prunn$pvaluesRightTail),]
df.export.Prunn.1p <- df.export.Prunn[df.export.Prunn$pvaluesRightTail<0.01,]
#Export into csv
write.csv(df.export.Prunn.1p,paste0(vcf.filename,"-FST.csv"),row.names = F)

##############################
### Manhattan Plot
##############################

#Genome Lines
MeanLine <- mean(df.export.Prunn$FSTNoCorr)
GenLine <- 0.01

if (vcf.filename == "Angus"){
  pallete <- c("#F8766D", "#7CAE00")
}  else if (vcf.filename == "Brahman"){
    pallete <- c("#00BFC4", "#7CAE00")
}  else if (vcf.filename == "HiLo"){
  pallete <- c("#d84d00", "#00d8ac")
}  else if (vcf.filename == "HoLo"){
  pallete <- c("#0028d8", "#00d8ac")
}  else if (vcf.filename == "HiHo"){
  pallete <- c("#d84d00", "#0028d8")
}  else {
  pallete <- c("#6c6960", "#403a3a")
}        

png(paste0(vcf.filename,"-Rtail.png"),width = 1366,height = 745)
manhattan(df.export.Prunn, chr="CHR", bp="POSITION", snp="POSITION",  p="pvaluesRightTail",
            suggestiveline = -log10(MeanLine),
            genomewideline = -log10(GenLine),
            col = pallete)
dev.off()
print(paste0("File name: ",vcf.filename," done Succesfully!"))
