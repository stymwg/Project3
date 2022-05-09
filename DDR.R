# This script first plots diversity vs differentiation residuals for population pairs tsv files.
# The script then identifies the upper 99th percentile of DDR outliers, and writes them to an excel file.
# Tsv files created using Gatk program installed into new conda environment, with command gatk VariantsToTable -V input.vcf -F CHROM -F POS -F TYPE -GF AD -GF DP -O output.table
# Script by: Levi Yant, altered by Matthew Gaskins
# Date: 22nd April 2022

rm(list=ls())
install.packages("readr")
library("readr")
install.packages("ggplot2")
library("ggplot2")
setwd("C:/Users/matth/OneDrive/Documents/Bioinformatics")
d1=read_table2("BRUC_GATK.tsv") #read in GATK-generated table from population 1
d2=read_table2("GAR_GATK.tsv") #read in GATK-generated table from population 2



wsize=25 #number of SNPs in window
dn=10 #number of alleles from population 1 (remember to consider ploidy - diploid in B.fruticulosa)
tn=10 #number of alles from population 2

d1$AC[is.na(d1$AC)]=0 #sets NA values equal to zero
d2$AC[is.na(d2$AC)]=0

d1$df=d1$AC/dn #calculates allele frequency for population 1
d1$tf=d2$AC/tn #calculates allele frequency for population 2
d1$diff=abs(d1$df-d1$tf) #calculates allele frequency difference between population 1 and 2
d1$tpi=(2*d1$tf*(1-d1$tf))*(tn/(tn-1)) #calculates nucleotide diversity
d1$tpi[d1$tpi<0] <- 0
 #diversityvsallelicdiff

w=matrix(data=0,ceiling(nrow(d1)/wsize),5)
wstart=1
wend=wsize

#plots diversity vs differentiation
for(i in 1:nrow(w)) {
	td=d1[wstart:wend,]
	w[i,1]=min(td$POS)
	w[i,2]=max(td$POS)
	w[i,3]=(w[i,2]-w[i,1])/2
	w[i,4]=mean(td$diff,na.rm=T)
	w[i,5]=mean(td$tpi,na.rm=T)
	wstart=wend+1
	wend=wend+wsize
	plot(w[,4], w[,5])
}

#Plots DDR density distribution
w=matrix(data=0,ceiling(nrow(d1)/wsize),6)
wstart=1
wend=wsize
for(i in 1:nrow(w)) {
  td=d1[wstart:wend,]
  w[i,1]=min(td$POS)
  w[i,2]=max(td$POS)
  w[i,3]=(w[i,2]-w[i,1])/2
  w[i,4]=mean(td$diff,na.rm=T)
  w[i,5]=mean(td$tpi,na.rm=T)
  w[i,6]=w[i,5]-w[i,4]
  wstart=wend+1
  wend=wend+wsize
}

#this code is in 25 identified SNP windows - NOT 25 nucleotide windows

hist(w[,6])
d <- density(w[,6]) #change y lab to frequency
plot(d)

#do 1% cutoff to find outliers
install.packages("dplyr")
library("dplyr")
cutoff <- round(nrow(w)*0.01) #generates cutoff value
w1 <- as.data.frame(w)
arranged_residuals <- desc(w1[,6]) #arranges from most negative DDR value to most positive
blended <- cbind(w1, arranged_residuals)
new_df <- blended %>% arrange(arranged_residuals)
final <- as.data.frame(new_df)
DDRoutliers_BRUC <- final[1:cutoff,] #extracts top 1% of DDR outliers

install.packages("writexl")
library("writexl")
write_xlsx(DDRoutliers_BRUC,"DDR_BRUCvsGAR_outliers.xlsx") #writes DDR outliers to an excel file


