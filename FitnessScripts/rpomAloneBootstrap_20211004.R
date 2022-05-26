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


# set # of trials and # of rows in spoFunction table as it is # of genes in genome
trials <- 10000
spoListLength <- nrow(spoFunction)
set.seed(15)

#initialize columns of table of final fitness values
spoFitnessTable <- as.numeric(rep(NA, spoListLength))
meanFitnessTable <- as.numeric(rep(NA, spoListLength))
lowerCL <- as.numeric(rep(NA, spoListLength))
upperCL <- as.numeric(rep(NA, spoListLength))
naCount <- as.numeric(rep(NA, spoListLength))

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

expansion_fitnessCalc <- rpomAloneExpansion

initial_treatmentColumns_meanFitnessCalc <- initial
treatmentColumns_meanFitnessCalc <- rpomAlone

outputName <- 'rpomAloneFitness_20211004.csv'

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

# fitnessCalc is looesly following calc_fitness.py, McCoy, K. M., Antonio, M. L., & van Opijnen, T. (2017). MAGenTA: a Galaxy implemented tool for complete Tn-Seq analysis and data visualization. Bioinformatics, 33(17), 2781-2783.
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


# Randomly chooses 4 reps with replacement for initial time, and for final time for the current gene of interest, calculates 4 different fitness values for the gene, averages them and stores the value. This function replicated x times for bootstrap
# spoSampling is the spo of interest, filled in by the loop below
# treatmentColumns_meanFitnessCalc is the final 4 reps of interest
# initial_treatmentColumns_meanFitnessCalc is filled in manually and is often fixed
# to the 1/2YTSS generated library
meanFitnessCalc <- function(spoSample = spoSampling, treatmentColumns = treatmentColumns_meanFitnessCalc){
    w <- as.numeric(rep(NA, 4))
    
    initialRep <- sample(initial_treatmentColumns_meanFitnessCalc, 4, replace=TRUE)
    finalRep <- sample(treatmentColumns, 4, replace=TRUE)
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
    
    w <- as.numeric(w[!is.infinite(w) & !is.na(w)])
    meanW <- mean(w, na.rm=TRUE)
}

#meanW <- replicate(10000, meanFitnessCalc(spoSample = geneSum$spo == 'SPO2849'))

# Bootstrap/confidence intervals calculated using method in Crawley, Michael J. The R book. John Wiley & Sons, 2012.
# Bootstrap code is adapted from http://strata.uga.edu/8370/lecturenotes/resampling.html

for(i in 1:nrow(spoFunction)){
    spoFit <-  as.character(spoFunction$spo[i])
    meanW <- as.numeric(rep(NA, trials))
    spoSampling <- geneSum$spo == spoFit
    meanW <- replicate(trials, meanFitnessCalc())
    meanWNA <- count(is.na(meanW))
    meanW <- meanW[!is.infinite(meanW) & !is.na(meanW)]
    meanFit <- mean(meanW, na.rm=TRUE)
        
    alpha <- 0.05
    lwrCL <- quantile(meanW, alpha/2)
    uprCL <- quantile(meanW, 1-alpha/2)
            
    meanFitnessTable[i] <- meanFit
    lowerCL[i] <- lwrCL
    upperCL[i] <- uprCL
    spoFitnessTable[i] <- spoFit
    naCount[i] <- meanWNA
}

fitness <- data.frame(matrix(ncol=5, nrow=nrow(spoFunction)))
colnames(fitness) <- c('mean', 'lwrCL', 'uprCL', 'spo', 'naCount')

fitness$mean <- round(meanFitnessTable, 3)
fitness$lwrCL <- round(lowerCL, 3)
fitness$uprCL <- round(upperCL, 3)
fitness$spo <- spoFitnessTable
fitness$naCount <- naCount

write.table(fitness, file=outputName, sep=',', quote=FALSE, col.names=TRUE, row.names=FALSE)
#
#fitness <- fread(outputName, header=TRUE)
#fitness 

#y <- fitness %>%  filter(., uprCL-lwrCL < 0.3)
#x <- 1:nrow(y)

##dev.new()
#ysort <- y[order(y$mean), ]
#plot(x, ysort$mean, type='n')
#points(x, ysort$mean[x], col='red', main = 'All')
#segments(y0=ysort$lwrCL[x], y1=ysort$uprCL[x], x0=x)
#abline(h=mean(y$mean), col='blue')
#abline(h=mean(y$mean) - sd(y$mean)*2, col='blue', lty=2)
#abline(h=mean(y$mean) + sd(y$mean)*2, col='blue', lty=2)

##dev.new()
#plot(y$mean[x], x, xlim=c(0, 1.5), type='n', main = 'all')

#segments(x0=y$lwrCL[x], x1=y$uprCL[x], y0=x, col='grey')
#points(y$mean[x], x, pch=16, col='red')


