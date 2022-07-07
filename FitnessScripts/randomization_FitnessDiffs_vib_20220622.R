library(tidyverse)
library(data.table)
library(matrixStats)
library(ggplot2)
library(rstatix)

# Set seed for reproducibility
set.seed(25)

# read in mean fitness data for the log2foldchange, read in the raw fitness values that were input to calculate mean fitness, read in data containing rpom genes and gene names
meanFitness <- fread('Data/meanFitness_20220701.csv')
fitness <- fread('Data/fitness_raw_20220701.csv')
spoFunction <- fread('DSS3_annotations_April2020.csv', header=TRUE)


# Take the mean fitness file and just select the fold changes as this is what we care about testing, join this data with the raw fitness data to make randomization code easy
meanFitness_fold <- meanFitness %>% 
  select(spo, Vib_fold_diff, Mari_fold_diff, CC_fold_diff)
fit_rando <- inner_join(fitness, meanFitness_fold, by='spo')


# These correspond with the data columns in the fit_rando file
rpomAlone_col <- 2:5
vib_col <- 6:9
mari_col <- 10:13
all_col <- 14:17

# Set empty tables to be used for the pvalue calculated and the table of genes tested
spoListLength <- nrow(fit_rando)
pVal <- as.numeric(rep(NA, spoListLength))
spoTable <- as.numeric(rep(NA, spoListLength))


# Function for calculating randomized log2fold change in means
# This code was adapted from http://strata.uga.edu/8370/lecturenotes/resampling.html
randomizedDifferenceInMeans <- function(rpomAlone_W=rpomAlone_W_test, treatment_W=treatment_W_test) {
  n_rpomAlone_W <- length(rpomAlone_W)
  n_treatment_W <- length(treatment_W)
  
  # combine the data
  combined <- c(rpomAlone_W, treatment_W)
  ncom <- length(combined)
  indices <- 1:ncom
  
  # initially assign all observations to y group
  group <- rep('treatment', ncom)
  
  # assign a subsample to to x group
  rpomAlone_sub <- sample(indices, n_rpomAlone_W)
  group[rpomAlone_sub] <- 'rpomAlone'
  
  # calculate the means
  meanW_rpomAlone <- mean(combined[group=='rpomAlone'], na.rm=TRUE)
  meanW_treatment <- mean(combined[group=='treatment'], na.rm=TRUE)
  LFC_diff <- round(log2((meanW_treatment+0.01)/(meanW_rpomAlone+0.01)),3)
  LFC_diff
}


# function to calculate pvalue from null distribution built up by randomization, at the end is multiplied by 2 to represent a 2-tailed test
# This code was adapted from http://strata.uga.edu/8370/lecturenotes/resampling.html
pvalue <- function(x, observed){
  # if the observed L2FC is greater than 0, take the length of all null observations equal to or more extreme (greater than) than the value observed
  # if the observed L2FC is less than 0, take the length of all null observations equal to or more extreme (smaller than) than the value observed
  if(observed > 0){
    nExtreme <- length(x[ x>= observed])
  } else{
    nExtreme <- length(x[ x<= observed])
  }
  # how many observations are there total? should be 10000
  # how many null observations was our observation more extreme than 
  nTotal <- length(x)
  p <- (nExtreme / nTotal)*2
  p
}



# loop to do randomization 10000 times, returns a p value for each gene
for(i in 1:nrow(fit_rando)){
  # assign a gene to our sampling
  spoFit <- fit_rando$spo[i]
  spoSampling <- fit_rando$spo == spoFit
  # if the observed L2FC of the treatment we are testing, Mari in this case, is not NA, then assign the fitness values from rpomAlone and the treatment columns to the test objects, then replicate the randomized difference in means function 10000 times to return the distribution, record the gene tested and record the pvalue by using the distribution tested against the observed L2FC of the treatment for this gene. If L2FC is NA, then assign a p-value of NA to that gene
  if(!is.na(fit_rando$Vib_fold_diff[fit_rando$spo==spoFit])){
    rpomAlone_W_test <- as.numeric(fit_rando[spoSampling, 2:5])
    treatment_W_test <- as.numeric(fit_rando[spoSampling, 6:9])
      LFC_dist <- replicate(10000, randomizedDifferenceInMeans(rpomAlone_W = rpomAlone_W_test, treatment_W=treatment_W_test))
    spoTable[i] <- spoFit
    pVal[i] <- pvalue(LFC_dist, fit_rando$Vib_fold_diff[fit_rando$spo==spoFit])
  } else {
    spoTable[i] <- spoFit
    pVal[i] <- NA
  }
}


# creates an empty dataframe to store the genes and associated p-values
fitnessSig <- data.frame(matrix(ncol=2, nrow=spoListLength))
colnames(fitnessSig) <- c('spo', 'pval')
fitnessSig$spo <- spoTable
fitnessSig$pval <- pVal


# combine the significance with the function of the genes, then do a benjamini-hochberg adjustment to the p-values for multiple corrections
fitnessSig_final <- inner_join(fitnessSig, spoFunction, by='spo')
fitnessSig_final$padj <- p.adjust(fitnessSig_final$pval, method="BH")

# write out the data as a table of the genes with associated significance
write.csv(fitnessSig_final, 'vib_W_rando_20220705.csv', row.names=FALSE)
