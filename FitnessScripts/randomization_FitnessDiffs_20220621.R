library(tidyverse)
library(data.table)
library(matrixStats)
library(ggplot2)
library(rstatix)

# Puts all the mean fitness data into one file called meanFitness, then adds the SPO functions

rpomAloneFitness <- fread('Data/rpomAloneFitness_raw_20220621.csv', header=TRUE)
vibFitness <- fread('Data/vibFitness_raw_20220621.csv', header=TRUE)
mariFitness <- fread('Data/mariFitness_raw_20220621.csv', header=TRUE)
compCommFitness <- fread('Data/compCommFitness_raw_20220621.csv', header=TRUE) 
spoFunction <- fread('DSS3_annotations_April2020.csv', header=TRUE)

Fitness <- rpomAloneFitness %>% 
  full_join(vibFitness, by='spo', suffix=c('_rpomAlone', '_vib'))  %>%   full_join(mariFitness, by='spo', suffix=c('_vib','_mari')) %>% 
  full_join(compCommFitness, by='spo', suffix=c('_mari', '_CC')) 


rpomAlone_W_test <- as.numeric(Fitness[Fitness$spo=='SPO0397', 2:5])
treatment_W_test <- as.numeric(Fitness[Fitness$spo=='SPO0397', 14:17])

# Function for calculating one randomized difference in means
randomizedDifferenceInMeans <- function(rpomAlone_W, treatment_W) {
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
  meanW_rpomAlone <- mean(combined[group=='rpomAlone'])
  meanW_treatment <- mean(combined[group=='treatment'])
  LFC_diff <- round(log2((meanW_rpomAlone+0.01)/(meanW_treatment+0.01)),3)
  LFC_diff
}

# Repeat that function many times to build the distribution
LFC_diff <- replicate(10000, randomizedDifferenceInMeans(rpomAlone_W = rpomAlone_W_test, treatment_W=treatment_W_test))

pvalue <- function(x, observed, tails=2){
  if(observed > 0){
    nExtreme <- length(x[ x>= observed])
  } else{
    nExtreme <- length(x[ x<= observed])
  }
  nTotal <- length(x)
  p <- nExtreme / nTotal * tails
  p
}

pvalue(x=LFC_diff, observed=5)
