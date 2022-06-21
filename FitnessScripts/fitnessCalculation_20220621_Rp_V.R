library(tidyverse)
library(data.table)
library(matrixStats)

insertionDensity <- fread('20211004_normalized_terminusRM_spo.csv', header=TRUE) %>% 
  select(-V1)

# Table of spo functions, attach to insertion Density
spoFunction <- fread('DSS3_annotations_April2020.csv', header=TRUE) %>% 
  select(spo,funct) %>% 
  filter(spo %in% insertionDensity$spo)
insertionDensity <- merge(insertionDensity, spoFunction, by='spo')


# Sum of all insertions grouped by spo (i.e. total reads per gene mutant)
geneSum <- insertionDensity %>% 
  select(-pos) %>% 
  pivot_longer(cols= -c(spo, funct), names_to='condition', values_to='counts') %>%
  group_by(condition, spo, funct) %>% 
  summarise(sum = sum(counts)) %>% 
  pivot_wider(names_from=condition, values_from=sum) %>% 
  ungroup()


# set # of rows in spoFunction table as it is # of genes in genome

spoListLength <- nrow(spoFunction)

#initialize columns of table of final fitness values
spoFitnessTable <- as.numeric(rep(NA, spoListLength))
w_1 <- as.numeric(rep(NA, spoListLength))
w_2 <- as.numeric(rep(NA, spoListLength))
w_3 <- as.numeric(rep(NA, spoListLength))
w_4 <- as.numeric(rep(NA, spoListLength))

# treatment's corresponding expansion factor and column totals
rpomAloneExpansion <- 38
vibExpansion <- 60
mariExpansion <- 91
completeCommExpansion <-65

treatmentTotals <- colSums(geneSum[, c(-1,-2)])


#column selection

all <- as.integer(3:6)
mari <- as.integer(7:10)
initial <- as.integer(11:14)
rpomAlone <- as.integer(15:18)
vib <- as.integer(19:22)

# Replace objects here to calculate fitness for different treatments

expansion_fitnessCalc <- vibExpansion

initial_treatmentColumns_meanFitnessCalc <- initial
treatmentColumns_meanFitnessCalc <- vib

outputName_raw <- 'vibFitness_raw_20220621.csv'
outputName <- 'vibAloneFitness_20220621.csv'


# calculating a fitness of a gene from sum of insertion reads

# set sum of insertion reads from T0 and TF 
# If reads1 is zero, it didn't exist, and average reads less than 20 are less reliable
# mt_freq_t1 = mutant frequenceny in the population initially
# mt_freq_t2 = mutant frequency in the population after selection
# pop_freq_t1 = frequency of the rest of the population initially
# pop_freq_t2 = frequency of the rest of the population after selection

# expansion = expansion factor of data at time zero vs the flow cytometry data at time final(corrected for percentage of Rpom vs others, see paper_building document for calculations)
# we do a log(x+1) to top and bottom so that there are no -Inf/Inf values
# w = fitness

# fitnessCalc is looesly following calc_fitness.py, credit to van Opjinen paper
# total1 is the total reads in sample 1
# total2 is the total reads in sample 2

fitnessCalc <- function(geneSumI, geneSumF, T0_total, Tf_total, expansion=expansion_fitnessCalc){
  reads1 <- geneSumI
  reads2 <- geneSumF
  if(reads1+reads2/2 >10) {
    mt_freq_t1 <- reads1/T0_total
    mt_freq_t2 <- reads2/Tf_total
    pop_freq_t1 <- 1 - mt_freq_t1
    pop_freq_t2 <- 1 - mt_freq_t2
    
    top <- log((mt_freq_t2*(expansion/mt_freq_t1))+1)
    bottom <- log((pop_freq_t2*(expansion/pop_freq_t1))+1)
    w <- top/bottom 
  } else {
    NA
  }
}


# Takes an initial and final time point for gene of interest and calculates 4 values of fitness for the gene, stores these values
 # spoSampling is the spo of interest, filled in by the loop below
 # treatmentColumns_meanFitnessCalc is the final 4 reps of interest
 # initial_treatmentColumns_meanFitnessCalc is filled in manually and is often fixed
 # to the 1/2YTSS generated library
 meanFitnessCalc <- function(spoSample = spoSampling, treatmentColumns = treatmentColumns_meanFitnessCalc){
       w <- as.numeric(rep(NA, 4))
       
       initialRep <- initial_treatmentColumns_meanFitnessCalc
       finalRep <- treatmentColumns
       initialRep1 <- initialRep[1]
       initialRep2 <- initialRep[2]
       initialRep3 <- initialRep[3]
       initialRep4 <- initialRep[4]
       finalRep1 <- finalRep[1]
       finalRep2 <- finalRep[2]
       finalRep3 <- finalRep[3]
       finalRep4 <- finalRep[4]
         
           
       geneSumI_1 <- as.numeric(geneSum[spoSample, initialRep1])
       geneSumI_2 <- as.numeric(geneSum[spoSample, initialRep2])
       geneSumI_3 <- as.numeric(geneSum[spoSample, initialRep3])
       geneSumI_4 <- as.numeric(geneSum[spoSample, initialRep4])
           
       geneSumF_1 <- as.numeric(geneSum[spoSample, finalRep1])
       geneSumF_2 <- as.numeric(geneSum[spoSample, finalRep2])
       geneSumF_3 <- as.numeric(geneSum[spoSample, finalRep3])
       geneSumF_4 <- as.numeric(geneSum[spoSample, finalRep4])
             
       w[1] <- fitnessCalc(geneSumI_1, geneSumF_1, treatmentTotals[initialRep1 -2], treatmentTotals[finalRep1 -2])
       w[2] <- fitnessCalc(geneSumI_2, geneSumF_2, treatmentTotals[initialRep2 -2], treatmentTotals[finalRep2 -2])
       w[3] <- fitnessCalc(geneSumI_3, geneSumF_3, treatmentTotals[initialRep3 -2], treatmentTotals[finalRep3 -2])
       w[4] <- fitnessCalc(geneSumI_4, geneSumF_4, treatmentTotals[initialRep4 -2], treatmentTotals[finalRep4 -2])
       return(w)
          
 }
 
 
for(i in 1:nrow(spoFunction)){
   spoFit <-  as.character(spoFunction$spo[i])
   spoSampling <- geneSum$spo == spoFit
   W <- meanFitnessCalc()
     
   spoFitnessTable[i] <- spoFit
   w_1[i] <- W[1]
   w_2[i] <- W[2]
   w_3[i] <- W[3]
   w_4[i] <- W[4]
}
 
 
fitness <- data.frame(matrix(ncol=5, nrow=nrow(spoFunction)))
colnames(fitness) <- c('spo', 'w_1', 'w_2', 'w_3', 'w_4')
 
fitness$spo <- spoFitnessTable
fitness$w_1 <- w_1
fitness$w_2 <- w_2
fitness$w_3 <- w_3
fitness$w_4 <- w_4

fitness$w_1[is.infinite(fitness$w_1)] <- NA 
fitness$w_2[is.infinite(fitness$w_2)] <- NA 
fitness$w_3[is.infinite(fitness$w_3)] <- NA 
fitness$w_4[is.infinite(fitness$w_4)] <- NA 

write.csv(fitness, outputName_raw, row.names=FALSE)

meanFit <- fitness %>% 
  pivot_longer(., cols=2:5, names_to='fit', values_to='W') %>% 
  group_by(spo) %>% 
  na.omit() %>% 
  summarize(meanW=mean(W), semW=sd(W)/sqrt(length(W)))

meanFit$meanW <- round(meanFit$meanW,3)
meanFit$semW <- round(meanFit$semW,3)

write.csv(meanFit, outputName, row.names=FALSE)

