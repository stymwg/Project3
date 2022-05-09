#===============================================================================
#                               Welcome!
# This script basically does a PCA plot for genomic data. It was written for
# plants so it can deal with different ploidies (which is very useful!)
# If a comment doesn't have name, it was made by Levi, otherwise it has the
# name of the person that has done it between hash tags
#
# Script by: Levi Yant
# altered by: Matthew Gaskins and Ana C. da Silva
# Date: 21st April 2022
#===============================================================================

install.packages("adegenet", dep=TRUE)
options(warn=1)

library(adegenet)
library(adegraphics) #not strictly necessary for all of this (homebrew R installs will interfere)
library(vcfR)
library(pegas)
library(StAMPP)
library(ade4)
library(MASS)
library(ggplot2)

######################=========MODIFIED FUNCTIONS=========######################

# a function for conversion from vcfR object to genlight in tetraploids
##Levi##: note not all of this is necessary for LIFE4136 project, but some is helpful
vcfR2genlight.tetra <- function (x, n.cores = 1) 
{
  bi <- is.biallelic(x)
  if (sum(!bi) > 0) {
    msg <- paste("Found", sum(!bi), "loci with more than two alleles.")
    msg <- c(msg, "\n", paste("Objects of class genlight only support loci with two alleles."))
    msg <- c(msg, "\n", paste(sum(!bi), "loci will be omitted from the genlight object."))
    warning(msg)
    x <- x[bi, ]
  }
  x <- addID(x)
  CHROM <- x@fix[, "CHROM"]
  POS <- x@fix[, "POS"]
  ID <- x@fix[, "ID"]
  x <- extract.gt(x)
  x[x == "0|0"] <- 0     #Ana# for diploids it's only lines 43 to 50
  x[x == "0|1"] <- 1
  x[x == "1|0"] <- 1
  x[x == "1|1"] <- 2
  x[x == "0/0"] <- 0
  x[x == "0/1"] <- 1
  x[x == "1/0"] <- 1
  x[x == "1/1"] <- 2
  x[x == "1/1/1/1"] <- 4
  x[x == "0/1/1/1"] <- 3
  x[x == "0/0/1/1"] <- 2
  x[x == "0/0/0/1"] <- 1
  x[x == "0/0/0/0"] <- 0
  x[x == "0/0/0/0/0/0"] <- 0
  x[x == "0/0/0/0/0/1"] <- 1
  x[x == "0/0/0/0/1/1"] <- 2
  x[x == "0/0/0/1/1/1"] <- 3
  x[x == "0/0/1/1/1/1"] <- 4
  x[x == "0/1/1/1/1/1"] <- 5
  x[x == "1/1/1/1/1/1"] <- 6
  if (requireNamespace("adegenet")) {
    x <- new("genlight", t(x), n.cores = n.cores)
  }
  else {
    warning("adegenet not installed")
  }
  adegenet::chromosome(x) <- CHROM
  adegenet::position(x) <- POS
  adegenet::locNames(x) <- ID
  return(x)
}

glPcaFast <- function(x,
                      center=TRUE,
                      scale=FALSE,
                      nf=NULL,
                      loadings=TRUE,
                      alleleAsUnit=FALSE,
                      returnDotProd=FALSE){
  
  if(!inherits(x, "genlight")) stop("x is not a genlight object")
  # keep the original mean / var code, as it's used further down
  # and has some NA checks...
  if(center) {
    vecMeans <- glMean(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecMeans))) stop("NAs detected in the vector of means")
  }
  if(scale){
    vecVar <- glVar(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecVar))) stop("NAs detected in the vector of variances")
  }
  # convert to full data, try to keep the NA handling as similar to the original as possible
  # - dividing by ploidy keeps the NAs
  mx <- t(sapply(x$gen, as.integer)) / ploidy(x)
  # handle NAs
  NAidx <- which(is.na(mx), arr.ind = T)
  if (center) {
    mx[NAidx] <- vecMeans[NAidx[,2]]
  } else {
    mx[NAidx] <- 0
  }
  # center and scale
  mx <- scale(mx,
              center = if (center) vecMeans else F,
              scale = if (scale) vecVar else F)
  # all dot products at once using underlying BLAS to support thousands of samples,
  # this could be replaced by 'Truncated SVD', but it would require more changes
  # in the code around
  allProd <- tcrossprod(mx) / nInd(x) # assume uniform weights
  
  ## PERFORM THE ANALYSIS ## ---------------------------------------------------
  # eigen analysis
  eigRes <- eigen(allProd, symmetric=TRUE, only.values=FALSE)
  rank <- sum(eigRes$values > 1e-12)
  eigRes$values <- eigRes$values[1:rank]
  eigRes$vectors <- eigRes$vectors[, 1:rank, drop=FALSE]
  # scan nb of axes retained
  if(is.null(nf)){
    barplot(eigRes$values, main="Eigenvalues", col=heat.colors(rank))
    cat("Select the number of axes: ")
    nf <- as.integer(readLines(n = 1))
  }
  # rescale PCs
  res <- list()
  res$eig <- eigRes$values
  nf <- min(nf, sum(res$eig>1e-10))
  ##res$matprod <- allProd # for debugging
  ## use: li = XQU = V\Lambda^(1/2)
  eigRes$vectors <- eigRes$vectors * sqrt(nInd(x)) # D-normalize vectors
  res$scores <- sweep(eigRes$vectors[, 1:nf, drop=FALSE],2, sqrt(eigRes$values[1:nf]), FUN="*")
  
  ## GET LOADINGS ## -----------------------------------------------------------
  # need to decompose X^TDV into a sum of n matrices of dim p*r
  # but only two such matrices are represented at a time
  if(loadings){
    if(scale) {
      vecSd <- sqrt(vecVar)
    }
    res$loadings <- matrix(0, nrow=nLoc(x), ncol=nf) # create empty matrix
    ## use: c1 = X^TDV
    ## and X^TV = A_1 + ... + A_n
    ## with A_k = X_[k-]^T v[k-]
    myPloidy <- ploidy(x)
    for(k in 1:nInd(x)){
      temp <- as.integer(x@gen[[k]]) / myPloidy[k]
      if(center) {
        temp[is.na(temp)] <- vecMeans[is.na(temp)]
        temp <- temp - vecMeans
      } else {
        temp[is.na(temp)] <- 0
      }
      if(scale){
        temp <- temp/vecSd
      }
      res$loadings <- res$loadings + matrix(temp) %*% eigRes$vectors[k, 1:nf, drop=FALSE]
    }
    res$loadings <- res$loadings / nInd(x) # don't forget the /n of X_tDV
    res$loadings <- sweep(res$loadings, 2, sqrt(eigRes$values[1:nf]), FUN="/")
  }
  ## FORMAT OUTPUT ## ----------------------------------------------------------
  colnames(res$scores) <- paste("PC", 1:nf, sep="")
  if(!is.null(indNames(x))){
    rownames(res$scores) <- indNames(x)
  } else {
    rownames(res$scores) <- 1:nInd(x)
  }
  if(!is.null(res$loadings)){
    colnames(res$loadings) <- paste("Axis", 1:nf, sep="")
    if(!is.null(locNames(x)) & !is.null(alleles(x))){
      rownames(res$loadings) <- paste(locNames(x),alleles(x), sep=".")
    } else {
      rownames(res$loadings) <- 1:nLoc(x)
    }
  }
  if(returnDotProd){
    res$dotProd <- allProd
    rownames(res$dotProd) <- colnames(res$dotProd) <- indNames(x)
  }
  res$call <- match.call()
  class(res) <- "glPca"
  return(res)
}
######################====================================######################

# IMPORT SNP data from VCF
#Ana# here Levi used the .gz file, but we are going to use the vcf directly
#vcf <- read.vcfR("BF_65mb_F4.4dg_no_new_PA_5.vcf.gz", nrows=10000)      # nrow = number of SNPs to read in
#setwd("C:/Users/acs/project3group2")  #Ana# it didn't let me do this for some reason!
vcf <- read.vcfR("4dg_noPA5.vcf")   #read in all data

# convert to genlight 	
aa.genlight <- vcfR2genlight.tetra(vcf)
#Ana# here the bit "tetra" just means that this could be used with tetraploids

locNames(aa.genlight) <- paste(vcf@fix[,1],vcf@fix[,2],sep="_")  # add real SNP.names
#Ana# ---
# We can use this feature to group the populations by their characteristics!
# However to do that, I had to rename the populations (on the vcf) so that the
# first two characters now mean inland (IN) or coastal (CO), and then there are 
# three characters that identify location, and one character that identifies sample nr
# So if you choose 2, you get Inland vs Coastal, if 5 grouping by location
# -------
pop(aa.genlight)<-substr(indNames(aa.genlight),1,2)  # add pop names: here pop names are first chars of ind name

#check    ====VERY IMPORTANT====
aa.genlight
indNames(aa.genlight)
ploidy(aa.genlight)

toRemove <- is.na(glMean(aa.genlight, alleleAsUnit = FALSE)) # TRUE where NA
which(toRemove) # position of entirely non-typed loci
aa.genlight_correct <- aa.genlight[, !toRemove]

pca.1 <- glPcaFast(aa.genlight_correct, nf=300)

scatter(pca.1, posi="bottomright")  # plot scatter with the individual labels
loadingplot(pca.1)

# proportion of explained variance by each axes
pca.1$eig[1]/sum(pca.1$eig) # proportion of variation explained by 1st axis
pca.1$eig[2]/sum(pca.1$eig) # proportion of variation explained by 2nd axis
pca.1$eig[3]/sum(pca.1$eig) # proportion of variation explained by 3rd axis
pca.1$eig[4]/sum(pca.1$eig) # proportion of variation explained by 4th axis
pca.1$eig[5]/sum(pca.1$eig) # proportion of variation explained by 5th axis
pca.1$eig[6]/sum(pca.1$eig) # proportion of variation explained by 6th axis
pca.1$eig[7]/sum(pca.1$eig) # proportion of variation explained by 7th axis

#Ana# ---
# there is a way to check which axis you need to consider (see below)
# In a random sample with 9 components the proportion of variance will be 1/9 = 0.111
# so every component that has a proportion of variance higher than 0.111 needs to be considered
# but here I'm not really sure if we want 1/29 = 0.034 (from populations) (PCA1 to PCA7)
# or 1/8=0.125 (from genotype - see line 43) - which only gives PCA1 as important
# -------

# just to see pops coloured in a palette
col <- funky(6) #Ana# adjust the number here to the total of locations you have
#Ana# you can adjust "xax" and "yax" for the PCAs you want to compare
s.class(pca.1$scores, pop(aa.genlight),  xax=1, yax=2, col=transp(col,.6))

# save nice figs
pdf ("PCA_all_SNPs.pdf", width=14, height=7)
#Ana# I've been trying to move that data label from the middle of the group for the last hour '(-_-)' if you know how to do it, please tell me!
g1 <- s.class(pca.1$scores, pop(aa.genlight), ppoints.col="black", xax=1, yax=2, leg = T, col=transp(col,.6), plabels = list(box = list(draw = FALSE)))
g1 #Ana# I like to see the graphs asap
g2 <- s.label (pca.1$scores, xax=1, yax=2, ppoints.col = "red", plabels = list(box = list(draw = FALSE), 
                                                                               optim = TRUE), paxes.draw=T, pgrid.draw =F, plabels.cex=1, plot = FALSE)
g2 #Ana# I like to see the graphs asap
ADEgS(c(g1, g2), layout = c(1, 2))
adev.off()

#Ana# this next part doesn't matter as ploidy is fixed here, so I commented everything!

# #ploidy - differentiated plots
# pdf ("PCA_all_ploidycol_SNPs_ax12_1K_less.pdf", width=14, height=7)
# g3 <- s.class(pca.1$scores, as.factor(as.vector(ploidy(aa.genlight))), xax=1, yax=2, col=transp(c("#FF0000", "#0000FF", "#00FF00")), 
#               ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=T, pgrid.draw =F, plab.cex = 0 , plot = FALSE)
# g4 <- s.label (pca.1$scores, xax=1, yax=2, ppoints.col = "red", plabels = list(box = list(draw = FALSE), 
#                                                                                optim = TRUE), paxes.draw=T, pgrid.draw =F, plabels.cex=1, plot = FALSE)
# ADEgS(c(g3, g4), layout = c(1, 2))
# dev.off()


