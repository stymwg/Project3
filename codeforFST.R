#This code first plots FST distributions between population pairs.
#The code then extracts the upper 99th percentile of FST outliers and writes them to an excel file.


#To obtain FST values, install vcftools into a new conda environment
#Then run command "vcftools --vcf FILE --weir-fst-pop pop1.txt --weir-fst-pop pop2.txt --out pop1_pop2.FST
#pop 1.txt is a single text file containing the name of every sample from population 1, each on an individual line
#pop 2.txt is a single text file containing the name of every sample from population 2, each on an individual line
#Created by Matthew Gaskins - Date - 22nd April 2022.

install.packages("readr")
library("readr")
setwd("C:/Users/matth/OneDrive/Documents/Bioinformatics")
FST <- read_table2("PAUvsLLAN.windowed.weir.fst") #insert your pairwise FST file here
library("dplyr")
new_df <- subset(FST, WEIGHTED_FST<1 & WEIGHTED_FST>0) #subsets all FST values between 0 and 1
new_df2 <- subset(new_df, N_VARIANTS>25) #extracts all windows with more than 25 SNPs
hist(new_df2$WEIGHTED_FST) #plots histogram
d <- density(new_df2$WEIGHTED_FST)
plot(d) #plots frequency distribution


cutoff <- round(nrow(new_df2)*0.01) #cutoff for 1% outliers - can be changed depending on desired outlier number

arrangedFST <- arrange(new_df2,desc(WEIGHTED_FST)) #arranges WEIGHTED FST from largest to smallest
FSToutliers_PAULLAN <- arrangedFST[1:cutoff,] #extracts 1% outliers


install.packages("writexl")
library("writexl")
write_xlsx(FSToutliers_PAULLAN,"FST_PAULLAN_outliers_vcftools.xlsx") #writes outliers to an excel file

plot(FST$N_VARIANTS, FST$WEIGHTED_FST)
corr <- cor.test(x=FST$WEIGHTED_FST, y=FST$N_VARIANTS, method = 'spearman') #calculates Spearman's Rank Correlation Coefficient
abline(lm(FST$WEIGHTED_FST ~ FST$N_VARIANTS), col = "red", lwd = 3) #plots correlation line

