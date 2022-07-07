library(tidyverse)
library(data.table)


 # significantFitness <- fread('Data/significantW_randomization_20220705.csv')
 # 
 # AA_COGs <- fread('COGs/AAtransportMetabolismCOG.txt') %>%
 #     inner_join(., significantFitness, by='spo') %>%
 #     mutate(., COG='AA')
 # 
 # Carbon_COGs <- fread('COGs/CarbontransportMetabolismCOG.txt') %>%
 #     inner_join(., significantFitness, by='spo') %>%
 #     mutate(., COG='Carbon')
 # 
 # Coenz_COGs <- fread('COGs/CoenzymeTransportMetabolismCOG.txt') %>%
 #     inner_join(., significantFitness, by='spo') %>%
 #     mutate(., COG='Coenzyme')
 # 
 # Energy_COGs <- fread('COGs/EnergyproductionConversion.txt') %>%
 #     inner_join(., significantFitness, by='spo') %>%
 #     mutate(., COG='Energy')
 # 
 # Inorg_COGs <- fread('COGs/InorganicTransportMetabolism.txt') %>%
 #     inner_join(., significantFitness, by='spo') %>%
 #     mutate(., COG='Inorganics')
 # 
 # Lipid_COGs <- fread('COGs/lipidTransportMetabolism.txt')  %>%
 #     inner_join(., significantFitness, by='spo') %>%
 #     mutate(., COG='lipids')
 # 
 # nucleotide_COGs <- fread('COGs/nucleotideTransportMetabolismCOG.txt') %>%
 #     inner_join(., significantFitness, by='spo') %>%
 #     mutate(., COG='nucleotide')
 # 
 # translation_COGs <- fread('COGs/translationRibosomalstructureBiogenesisCOG.txt') %>%
 #     inner_join(., significantFitness, by='spo') %>%
 #     mutate(., COG='translation')
 # 
 # 
 # COGtable <- AA_COGs %>%
 #     rbind(Carbon_COGs) %>%
 #     rbind(Coenz_COGs) %>%
 #     rbind(Energy_COGs) %>%
 #     rbind(Inorg_COGs) %>%
 #     rbind(Lipid_COGs) %>%
 #     rbind(translation_COGs) %>%
 #     rbind(nucleotide_COGs)
 # 
 # write.csv(COGtable, 'significantCOGfitness_20220705.csv')

# x <- merge(significantFitness, metabolism, all=TRUE)

metabolism <- fread('significantFitness_Metabolism_20220706.csv')


anabolic <- subset(metabolism, Metabolism == 'Anabolic') 
anabolicfold <- c(anabolic$Mari_fold_diff[anabolic$sig_mari=='yes'], anabolic$Vib_fold_diff[anabolic$sig_vib=='yes'], anabolic$CC_fold_diff[anabolic$sig_cc=='yes']) %>% 
    na.omit()

anabolicW <- c(anabolic$meanW_mari[anabolic$sig_mari=='yes'], anabolic$meanW_vib[anabolic$sig_vib=='yes'], anabolic$meanW_CC[anabolic$sig_cc=='yes']) %>% 
    na.omit()
anabolicWAlone <- anabolic$meanW_rpomAlone


catabolic <- subset(metabolism, Metabolism=='Catabolic')

catabolicfold <- c(catabolic$Mari_fold_diff[catabolic$sig_mari=='yes'], catabolic$Vib_fold_diff[catabolic$sig_vib=='yes'], catabolic$CC_fold_diff[catabolic$sig_cc=='yes']) %>% 
    na.omit()

catabolicW <- c(catabolic$meanW_mari[catabolic$sig_mari=='yes'], catabolic$meanW_vib[catabolic$sig_vib=='yes'], catabolic$meanW_CC[catabolic$sig_cc=='yes']) %>% 
    na.omit()
catabolicWAlone <- catabolic$meanW_rpomAlone

both <- subset(metabolism, Metabolism=='Amphibolic')
bothfold <- c(both$Mari_fold_diff[both$sig_mari=='yes'], both$Vib_fold_diff[both$sig_vib=='yes'], both$CC_fold_diff[both$sig_cc=='yes']) %>% 
    na.omit()
bothW <- c(both$meanW_mari[both$sig_mari=='yes'], both$meanW_vib[both$sig_vib=='yes'], both$meanW_CC[both$sig_cc=='yes']) %>% 
    na.omit()
bothWAlone <- both$meanW_rpomAlone


transport <- subset(metabolism, Metabolism=='Transport')
transportfold <- c(transport$Mari_fold_diff[transport$sig_mari=='yes'], transport$Vib_fold_diff[transport$sig_vib=='yes'], transport$CC_fold_diff[transport$sig_cc=='yes']) %>% 
    na.omit()
transportW <- c(transport$meanW_mari[transport$sig_mari=='yes'], transport$meanW_vib[transport$sig_vib=='yes'], transport$meanW_CC[transport$sig_cc=='yes']) %>% 
    na.omit()
transportWAlone <- transport$meanW_rpomAlone

anabolicTest <- t.test(anabolicfold)
anabolicWTest <- t.test(anabolicW)
anabolicWAloneTest <- t.test(anabolicWAlone)
catabolicTest <- t.test(catabolicfold)
catabolicWTest <- t.test(catabolicW)
catabolicWAloneTest <- t.test(catabolicWAlone)
bothTest <- t.test(bothfold)
bothWTest <- t.test(bothW)
bothWAloneTest <- t.test(bothWAlone)
transportTest <- t.test(transportfold)
transportWTest <- t.test(transportW)
transportWAloneTest <- t.test(transportWAlone)

combodW <- c(anabolicfold, catabolicfold, transportfold, bothfold)
comboW <- c(anabolicW, anabolicWAlone, catabolicW, catabolicWAlone)


dotPlotter <- function(metabolism, metabolismTest, X, meanColor='black'){
    points(runif(length(metabolism), X-0.01, X+0.01), metabolism, col='grey85')
    points(X, metabolismTest$estimate, pch=16, cex=1.2, col=meanColor)
    segments(x0=X, y0=metabolismTest$conf.int[1], y1=metabolismTest$conf.int[2], lwd=1.2, col=meanColor)
    
}


#pdf('metabolic_FoldDifferences_20220706.pdf', useDingbats=FALSE)


plot(c(0.8, 1.1), c(-8,8), type='n', xlim=c(0,0.8), xlab='', ylab='Log2_FC', xaxt='n', bty='n', las=2)
axis(1, at=c(0.1, 0.3, 0.5, 0.7), labels=c('Tport', 'Cata', 'Amp', 'Ana'))
dotPlotter(transportfold, transportTest, 0.1)
dotPlotter(catabolicfold, catabolicTest, 0.3)
dotPlotter(bothfold, bothTest, 0.5)
dotPlotter(anabolicfold, anabolicTest, 0.7)



#dev.off()
