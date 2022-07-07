library(tidyverse)
library(data.table)

meanFitness <- fread('Data/meanFitness_20220701.csv')

fold <- meanFitness %>% 
  select(spo, Vib_fold_diff, Mari_fold_diff, CC_fold_diff)

vibSig <- fread('Data/vib_W_rando_20220705.csv') %>% 
  subset(., padj <= 0.05) 
vibSig$sig_vib <- 'yes'
colnames(vibSig)[c(2,6)] <- paste(colnames(vibSig)[c(2,6)],"vib",sep="_")
vibSig <- inner_join(vibSig, fold) %>%
  subset(., abs(Vib_fold_diff) >0.100)


mariSig <- fread('Data/mari_W_rando_20220705.csv') %>% 
  subset(., padj <= 0.05)
mariSig$sig_mari <- 'yes'
colnames(mariSig)[c(2,6)] <- paste(colnames(mariSig)[c(2,6)],"mari",sep="_")
mariSig <- inner_join(mariSig, fold) %>%
  subset(., abs(Mari_fold_diff) >0.100)

ccSig <- fread('Data/comp_W_rando_20220705.csv') %>% 
  subset(., padj <= 0.05)
ccSig$sig_cc <- 'yes'
colnames(ccSig)[c(2,6)] <- paste(colnames(ccSig)[c(2,6)],"cc",sep="_")
ccSig <- inner_join(ccSig, fold) %>%
  subset(., abs(CC_fold_diff) >0.100)


significantW <- vibSig %>% 
  full_join(.,mariSig, by=c('spo', 'gene', 'funct', 'NCBI', 'Vib_fold_diff', 'Mari_fold_diff', 'CC_fold_diff' )) %>% 
  full_join(.,ccSig, by=c('spo', 'gene', 'funct', 'NCBI', 'Vib_fold_diff', 'Mari_fold_diff', 'CC_fold_diff' ))

sig <- inner_join(significantW, meanFitness) 


write.csv(sig, 'Data/significantW_randomization_20220705.csv', row.names=FALSE)

