library(tidyverse)
library(data.table)
library(matrixStats)
library(ggplot2)
library(rstatix)

# Puts all the bootstrap estimates into one file called meanFitness, then adds the SPO functions, removes bootstrap estimates that have a mean NA count across treatments >=5000 meaning the bootstrap estimate failed more than 5000 times due to a low number of mutant counts for that gene

rpomAloneFitness <- fread('Data/rpomAloneFitness_20220531.csv', header=TRUE)
rpomAloneFitness$Treatment <- 'rpomAlone'
vibFitness <- fread('Data/vibFitness_20220531.csv', header=TRUE)
vibFitness$Treatment <- 'vib'
mariFitness <- fread('Data/mariFitness_20220531.csv', header=TRUE)
compCommFitness <- fread('Data/compCommFitness_20220531.csv', header=TRUE) 
spoFunction <- fread('DSS3_annotations_April2020.csv', header=TRUE)

meanFitness <- rpomAloneFitness %>% 
  full_join(vibFitness, by='spo', suffix=c('_rpomAlone', '_vib'))  %>%   full_join(mariFitness, by='spo', suffix=c('_vib','_mari')) %>% 
  full_join(compCommFitness, by='spo', suffix=c('_mari', '_CC')) 

meanFitness <- merge(meanFitness, spoFunction, by='spo') %>% 
  group_by(spo) %>% 
  mutate(meanNA=mean(c(naCount_rpomAlone, naCount_vib, naCount_mari, naCount_CC))) %>% 
  subset(meanNA <= 5000) %>% 
  ungroup()

# Shep suggested reporting the log2 fold change instead of the difference in means

meanFitness$Vib_fold_diff <- round(log2((meanFitness$mean_vib+0.01)/(meanFitness$mean_rpomAlone+0.01)),3)

meanFitness$Mari_fold_diff <- round(log2((meanFitness$mean_mari+0.01)/(meanFitness$mean_rpomAlone+0.01)),3)

meanFitness$CC_fold_diff <- round(log2((meanFitness$mean_CC+0.01)/(meanFitness$mean_rpomAlone+0.01)), 3)


# Reports which data are statistically significant thru non-overlapping Bonferroni adjusted confidence intervals, puts them into a table together of significant genes

rpomAlone_Vib <- meanFitness$spo[meanFitness$spo %in% meanFitness$spo[meanFitness$uprCL_adj_rpomAlone < meanFitness$lwrCL_adj_vib | meanFitness$lwrCL_adj_rpomAlone > meanFitness$uprCL_adj_vib]]
rpomAlone_Mari <-  meanFitness$spo[meanFitness$spo %in% meanFitness$spo[meanFitness$uprCL_adj_rpomAlone < meanFitness$lwrCL_adj_mari | meanFitness$lwrCL_adj_rpomAlone > meanFitness$uprCL_adj_mari]]
rpomAlone_CC <-  meanFitness$spo[meanFitness$spo %in% meanFitness$spo[meanFitness$uprCL_adj_rpomAlone < meanFitness$lwrCL_adj_CC | meanFitness$lwrCL_adj_rpomAlone > meanFitness$uprCL_adj_CC]]

Vib <- meanFitness %>% filter(., meanFitness$spo %in% rpomAlone_Vib) %>% 
    subset(., Vib_fold_diff>0.1 | Vib_fold_diff < -0.1)
Vib$Vib_sig <- 'yes'
Mari <- meanFitness %>% filter(., meanFitness$spo %in% rpomAlone_Mari) %>% 
    subset(., Mari_fold_diff>0.1 | Mari_fold_diff < -0.1)
Mari$Mari_sig <- 'yes'

compComm <- meanFitness %>% filter(., meanFitness$spo %in% rpomAlone_CC) %>% 
    subset(., CC_fold_diff>0.1 | CC_fold_diff < -0.1)
compComm$CC_sig <- 'yes'

significantFitness <- full_join(Vib, Mari) %>% 
  full_join(., compComm) 
# 
# write.csv(meanFitness, 'meanFitness_20220601.csv', row.names=FALSE)
# write.csv(significantFitness, 'significant_genes_20220601.csv', row.names=FALSE)
# 
# missingGenes <- filter(spoFunction, !(spoFunction$spo %in% meanFitness$spo))
# 
# write.csv(missingGenes, 'missingMutants20220601.csv')



# 
# sigFit <- fread('significant_genes_grouped_Metabolism_filtered_20220413.csv') %>%
#     select(spo, gene, Group, Category, Subcategory, Metabolism, COG, COG_Category, COG_Name)
# significant <- inner_join(significantFitness, sigFit, by=c('spo', 'gene'))
# 
# write.csv(significant, 'significant_genes_grouped_Metabolism_filtered_20220601.csv', row.names=FALSE)
