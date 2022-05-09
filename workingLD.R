setwd("C:/Users/matth/OneDrive/Documents/Bioinformatics")
install.packages("readxl")
library("readxl")
table <- read_excel("long.xlsx") #excel document in three columns: column 1 is contig identity, column 2 is distance between SNPs and column 3 is rsquared
table[is.na(table)] <- 0 #removes NA values
library(ggplot2)
library(dplyr)
install.packages("devtools")
library(devtools)
devtools::install_github("kassambara/easyGgplot2", force = TRUE)
library(easyGgplot2)
table <- as.data.frame(table)
ggplot2.scatterplot(data=table, xName="Distance",yName="rsquared", 
                    size=3, groupName="Contig", groupColors=c('black','pink','orange', 'yellow', 'turquoise', 'red'), addRegLine=TRUE, smoothingMethod = "loess", regLineSize = 1.5) #plots LD for each contig on same graph

