# This script plots AFD values within chosen windows.
# To adapt to your samples, simply use find and replace 'BRUC' for population 1 and 'GAR' for population 2
# Tsv files created using Gatk program installed into new conda environment, with command gatk VariantsToTable -V input.vcf -F CHROM -F POS -F TYPE -GF AD -GF DP -O output.table
# Script by Matthew Gaskins
# Date: 25th April 2022

setwd("C:/Users/matth/OneDrive/Documents/Bioinformatics")
library(readr)
rm(list=ls())
# to start, set up a new environment by calling some packages:
library(ggplot2)
library(tidyverse)

# create a dataframe with the txt file (or the vcf file!)
pop1 <- read_table2("BRUC_GATK.tsv")
pop2 <- read_table2("GAR_GATK.tsv")

results1 <- pop1 %>% filter(CHROM == "ctg1") #set contig identity to desired
results2 <- pop2 %>% filter(CHROM == "ctg1")
BRUC <- c("CHROM", "POS", "BRUC1.AD", "BRUC1.DP", "BRUC2.AD", "BRUC2.DP", "BRUC3.DP", "BRUC3.AD", "BRUC4.AD", "BRUC4.DP", "BRUC5.DP", "BRUC5.AD") #extracts all columns with this heading - CHANGE FOR DIFFERENT SAMPLE NAMES
GAR <- c("CHROM", "POS", "GAR1.AD", "GAR1.DP", "GAR2.AD", "GAR2.DP", "GAR3.DP", "GAR3.AD", "GAR4.AD", "GAR4.DP", "GAR5.DP", "GAR5.AD") #extracts all columns with this heading - CHANGE FOR DIFFERENT SAMPLE NAMES
BRUC_table <- results1[,BRUC]
GAR_table <- results2[,GAR]

BRUC_table$corrected_BRUC1.AD <- as.numeric(sub(".*,", "", BRUC_table$BRUC1.AD)) #each line of this region of code creates a new column containing the values after the comma from the AD column for each individual BRUC
BRUC_table$corrected_BRUC2.AD <- as.numeric(sub(".*,", "", BRUC_table$BRUC2.AD))
BRUC_table$corrected_BRUC3.AD <- as.numeric(sub(".*,", "", BRUC_table$BRUC3.AD))
BRUC_table$corrected_BRUC4.AD <- as.numeric(sub(".*,", "", BRUC_table$BRUC4.AD))
BRUC_table$corrected_BRUC5.AD <- as.numeric(sub(".*,", "", BRUC_table$BRUC5.AD))

GAR_table$corrected_GAR1.AD <- as.numeric(sub(".*,", "", GAR_table$GAR1.AD)) #each line of this region of code creates a new column containing the values after the comma from the AD column for each individual BRUC
GAR_table$corrected_GAR2.AD <- as.numeric(sub(".*,", "", GAR_table$GAR2.AD))
GAR_table$corrected_GAR3.AD <- as.numeric(sub(".*,", "", GAR_table$GAR3.AD))
GAR_table$corrected_GAR4.AD <- as.numeric(sub(".*,", "", GAR_table$GAR4.AD))
GAR_table$corrected_GAR5.AD <- as.numeric(sub(".*,", "", GAR_table$GAR5.AD))

BRUC_table$total_AD <- BRUC_table$corrected_BRUC1.AD + BRUC_table$corrected_BRUC2.AD + BRUC_table$corrected_BRUC3.AD + BRUC_table$corrected_BRUC4.AD + BRUC_table$corrected_BRUC5.AD
BRUC_table$total_DP <- BRUC_table$BRUC1.DP + BRUC_table$BRUC2.DP + BRUC_table$BRUC3.DP + BRUC_table$BRUC4.DP + BRUC_table$BRUC5.DP
BRUC_table$BRUCAF <- BRUC_table$total_AD / BRUC_table$total_DP #creates total BRUC allele frequency values by dividing total AD values from all BRUC populations by the total DP values from all BRUC populations

GAR_table$total_AD <- GAR_table$corrected_GAR1.AD + GAR_table$corrected_GAR2.AD + GAR_table$corrected_GAR3.AD + GAR_table$corrected_GAR4.AD + GAR_table$corrected_GAR5.AD
GAR_table$total_DP <- GAR_table$GAR1.DP + GAR_table$GAR2.DP + GAR_table$GAR3.DP + GAR_table$GAR4.DP + GAR_table$GAR5.DP
GAR_table$GARAF <- GAR_table$total_AD / GAR_table$total_DP #creates total GAR allele frequency values by dividing total AD values from all GAR populations by the total DP values from all GAR populations
GAR_table$AFD <- abs(GAR_table$GARAF - BRUC_table$BRUCAF) #outputs the positive difference between allele frequencies of BRUC and tetraploid populations at each scaffold position

plot(GAR_table$POS, GAR_table$AFD, ylab = "Absolute AFD", xlab = "Scaffold Position", xlim=c(70000,75000), xaxs="i") #plots the absolute allele frequency differences between paired populations
