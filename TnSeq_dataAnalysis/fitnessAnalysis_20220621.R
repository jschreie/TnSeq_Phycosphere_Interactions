library(tidyverse)
library(data.table)
library(matrixStats)
library(ggplot2)
library(rstatix)

# Puts all the mean fitness data into one file called meanFitness, then adds the SPO functions

rpomAloneFitness <- fread('Data/rpomAloneFitness_20220621.csv', header=TRUE)
vibFitness <- fread('Data/vibFitness_20220621.csv', header=TRUE)
mariFitness <- fread('Data/mariFitness_20220621.csv', header=TRUE)
compCommFitness <- fread('Data/compCommFitness_20220621.csv', header=TRUE) 
spoFunction <- fread('DSS3_annotations_April2020.csv', header=TRUE)

meanFitness <- rpomAloneFitness %>% 
  full_join(vibFitness, by='spo', suffix=c('_rpomAlone', '_vib'))  %>%   full_join(mariFitness, by='spo', suffix=c('_vib','_mari')) %>% 
  full_join(compCommFitness, by='spo', suffix=c('_mari', '_CC')) 

meanFitness <- merge(meanFitness, spoFunction, by='spo')
meanFitness[, 2:9] <- signif(meanFitness[, 2:9], 2)

# puts all raw fitness values into a file called fitness_raw to be used for randomization

rpomAloneFitness <- fread('Data/rpomAloneFitness_raw_20220621.csv', header=TRUE)
vibFitness <- fread('Data/vibFitness_raw_20220621.csv', header=TRUE)
mariFitness <- fread('Data/mariFitness_raw_20220621.csv', header=TRUE)
compCommFitness <- fread('Data/compCommFitness_raw_20220621.csv', header=TRUE) 
spoFunction <- fread('DSS3_annotations_April2020.csv', header=TRUE)


Fitness <- rpomAloneFitness %>% 
  full_join(vibFitness, by='spo', suffix=c('_rpomAlone', '_vib'))  %>%   full_join(mariFitness, by='spo', suffix=c('_vib','_mari')) %>% 
  full_join(compCommFitness, by='spo', suffix=c('_mari', '_CC')) 

Fitness[, 2:17] <- signif(Fitness[, 2:17], 2)


# Shep suggested reporting the log2 fold change instead of the difference in means

meanFitness$Vib_fold_diff <- signif(log2((meanFitness$meanW_vib+0.01)/(meanFitness$meanW_rpomAlone+0.01)),2)

meanFitness$Mari_fold_diff <- signif(log2((meanFitness$meanW_mari+0.01)/(meanFitness$meanW_rpomAlone+0.01)),2)

meanFitness$CC_fold_diff <- signif(log2((meanFitness$meanW_CC+0.01)/(meanFitness$meanW_rpomAlone+0.01)), 2)


# 
write.csv(meanFitness, 'meanFitness_20220622.csv', row.names=FALSE)
# write.csv(significantFitness, 'significant_genes_20220601.csv', row.names=FALSE)
# 
missingGenes <- filter(spoFunction, !(spoFunction$spo %in% meanFitness$spo))
# 
write.csv(missingGenes, 'missingMutants20220622.csv', row.names=FALSE)


write.csv(Fitness, 'fitness_raw_20220622.csv', row.names=FALSE)
# 


vibSig <- fread('vib_W_rando.csv') %>% 
  subset(., padj <= 0.05) 
vibSig$sig_vib <- 'yes'
colnames(vibSig)[c(2,6)] <- paste(colnames(vibSig)[c(2,6)],"vib",sep="_")


mariSig <- fread('mari_W_rando.csv') %>% 
  subset(., padj <= 0.05)
mariSig$sig_mari <- 'yes'
colnames(mariSig)[c(2,6)] <- paste(colnames(mariSig)[c(2,6)],"mari",sep="_")

ccSig <- fread('comp_W_rando.csv') %>% 
  subset(., padj <= 0.05)
ccSig$sig_cc <- 'yes'
colnames(ccSig)[c(2,6)] <- paste(colnames(ccSig)[c(2,6)],"cc",sep="_")


significantW <- vibSig %>% 
  full_join(mariSig, by=c('spo', 'gene', 'funct', 'NCBI')) %>% 
  full_join(ccSig, by=c('spo', 'gene', 'funct', 'NCBI')) %>% 
  inner_join(., meanFitness) %>% 
  filter(., abs(Vib_fold_diff) >0.1 | abs(Mari_fold_diff) >0.1 | abs(CC_fold_diff) >0.1)

write.csv(significantW, 'significantW_randomization.csv', row.names=FALSE)
# sigFit <- fread('significant_genes_grouped_Metabolism_filtered_20220413.csv') %>%
#     select(spo, gene, Group, Category, Subcategory, Metabolism, COG, COG_Category, COG_Name)
# significant <- inner_join(significantFitness, sigFit, by=c('spo', 'gene'))
# 
# write.csv(significant, 'significant_genes_grouped_Metabolism_filtered_20220601.csv', row.names=FALSE)
