library(tidyverse)
library(data.table)
library(matrixStats)
library(ggplot2)
library(rstatix)

# Puts all the mean fitness data into one file called meanFitness, then adds the SPO functions

rpomAloneFitness <- fread('Data/rpomAloneFitness_20220701.csv', header=TRUE)
vibFitness <- fread('Data/vibFitness_20220701.csv', header=TRUE)
mariFitness <- fread('Data/mariFitness_20220701.csv', header=TRUE)
compCommFitness <- fread('Data/compCommFitness_20220701.csv', header=TRUE) 
spoFunction <- fread('DSS3_annotations_April2020.csv', header=TRUE)

meanFitness <- rpomAloneFitness %>% 
  full_join(vibFitness, by='spo', suffix=c('_rpomAlone', '_vib'))  %>%   full_join(mariFitness, by='spo', suffix=c('_vib','_mari')) %>% 
  full_join(compCommFitness, by='spo', suffix=c('_mari', '_CC')) 
meanFitness <- merge(meanFitness, spoFunction, by='spo')


# puts all raw fitness values into a file called fitness_raw to be used for randomization
rpomAloneFitness <- fread('Data/rpomAloneFitness_raw_20220701.csv', header=TRUE)
vibFitness <- fread('Data/vibFitness_raw_20220701.csv', header=TRUE)
mariFitness <- fread('Data/mariFitness_raw_20220701.csv', header=TRUE)
compCommFitness <- fread('Data/compCommFitness_raw_20220701.csv', header=TRUE) 
spoFunction <- fread('DSS3_annotations_April2020.csv', header=TRUE)


Fitness <- rpomAloneFitness %>% 
  full_join(vibFitness, by='spo', suffix=c('_rpomAlone', '_vib'))  %>%   full_join(mariFitness, by='spo', suffix=c('_vib','_mari')) %>% 
  full_join(compCommFitness, by='spo', suffix=c('_mari', '_CC')) 



# Shep suggested reporting the log2 fold change instead of the difference in means

meanFitness$Vib_fold_diff <- signif(log2((meanFitness$meanW_vib+0.01)/(meanFitness$meanW_rpomAlone+0.01)),2)

meanFitness$Mari_fold_diff <- signif(log2((meanFitness$meanW_mari+0.01)/(meanFitness$meanW_rpomAlone+0.01)),2)

meanFitness$CC_fold_diff <- signif(log2((meanFitness$meanW_CC+0.01)/(meanFitness$meanW_rpomAlone+0.01)), 2)


# 
write.csv(meanFitness, 'meanFitness_20220701.csv', row.names=FALSE)
# 
missingGenes <- filter(spoFunction, !(spoFunction$spo %in% meanFitness$spo))
# 
write.csv(missingGenes, 'missingMutants20220701.csv', row.names=FALSE)


write.csv(Fitness, 'fitness_raw_20220701.csv', row.names=FALSE)
#
