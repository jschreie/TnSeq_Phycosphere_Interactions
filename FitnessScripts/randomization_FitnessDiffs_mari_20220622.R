library(tidyverse)
library(data.table)
library(matrixStats)
library(ggplot2)
library(rstatix)

# Puts all the mean fitness data into one file called meanFitness, then adds the SPO functions

set.seed(20)

meanFitness <- fread('Data/meanFitness_20220622.csv')

fitness <- fread('Data/fitness_raw_20220622.csv')
spoFunction <- fread('DSS3_annotations_April2020.csv', header=TRUE)

rpomAlone_col <- 2:5
vib_col <- 6:9
mari_col <- 10:13
all_col <- 14:17



meanFitness_fold <- meanFitness %>% 
  select(spo, Vib_fold_diff, Mari_fold_diff, CC_fold_diff)
fit_rando <- inner_join(fitness, meanFitness_fold, by='spo')

spoListLength <- nrow(fit_rando)
pVal <- as.numeric(rep(NA, spoListLength))
spoTable <- as.numeric(rep(NA, spoListLength))

# Function for calculating one randomized difference in means
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

# function to calculate pvalue from null distribution built up by randomization
pvalue <- function(x, observed){
  if(observed > 0){
    nExtreme <- length(x[ x>= observed])
  } else{
    nExtreme <- length(x[ x<= observed])
  }
  nTotal <- length(x)
  p <- (nExtreme / nTotal)* 2
  p
}


for(i in 1:nrow(fit_rando)){
  spoFit <- fit_rando$spo[i]
  spoSampling <- fit_rando$spo == spoFit
  if(!is.na(fit_rando$Mari_fold_diff[fit_rando$spo==spoFit])){
    rpomAlone_W_test <- as.numeric(fit_rando[spoSampling, 2:5])
    treatment_W_test <- as.numeric(fit_rando[spoSampling, 10:13])
      LFC_dist <- replicate(10000, randomizedDifferenceInMeans(rpomAlone_W = rpomAlone_W_test, treatment_W=treatment_W_test))
    spoTable[i] <- spoFit
    pVal[i] <- pvalue(LFC_dist, fit_rando$Mari_fold_diff[fit_rando$spo==spoFit])
  } else {
    spoTable[i] <- spoFit
    pVal[i] <- NA
  }
}


fitnessSig <- data.frame(matrix(ncol=2, nrow=spoListLength))
colnames(fitnessSig) <- c('spo', 'pval')

fitnessSig$spo <- spoTable
fitnessSig$pval <- pVal

fitnessSig_final <- inner_join(fitnessSig, spoFunction, by='spo')

fitnessSig_final$padj <- p.adjust(fitnessSig_final$pval, method="BH")


write.csv(fitnessSig_final, 'mari_W_rando.csv', row.names=FALSE)
