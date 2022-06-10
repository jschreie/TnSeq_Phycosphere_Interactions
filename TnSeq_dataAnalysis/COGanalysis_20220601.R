library(tidyverse)
library(data.table)


metabolism <- fread('significant_genes_grouped_Metabolism_filtered_20220601.csv')


anabolic <- subset(metabolism, Metabolism == 'Anabolic') 
anabolicfold <- c(anabolic$Mari_fold_diff[anabolic$Mari_sig=='yes'], anabolic$Vib_fold_diff[anabolic$Vib_sig=='yes'], anabolic$CC_fold_diff[anabolic$CC_sig=='yes']) %>% 
    na.omit()

anabolicW <- c(anabolic$mean_mari[anabolic$Mari_sig=='yes'], anabolic$mean_vib[anabolic$Vib_sig=='yes'], anabolic$mean_CC[anabolic$CC_sig=='yes']) %>% 
    na.omit()
anabolicWAlone <- anabolic$mean_rpomAlone


catabolic <- subset(metabolism, Metabolism=='Catabolic')

catabolicfold <- c(catabolic$Mari_fold_diff[catabolic$Mari_sig=='yes'], catabolic$Vib_fold_diff[catabolic$Vib_sig=='yes'], catabolic$CC_fold_diff[catabolic$CC_sig=='yes']) %>% 
    na.omit()

catabolicW <- c(catabolic$mean_mari[catabolic$Mari_sig=='yes'], catabolic$mean_vib[catabolic$Vib_sig=='yes'], catabolic$mean_CC[catabolic$CC_sig=='yes']) %>% 
    na.omit()
catabolicWAlone <- catabolic$mean_rpomAlone

both <- subset(metabolism, Metabolism=='Amphibolic')
bothfold <- c(both$Mari_fold_diff[both$Mari_sig=='yes'], both$Vib_fold_diff[both$Vib_sig=='yes'], both$CC_fold_diff[both$CC_sig=='yes']) %>% 
    na.omit()
bothW <- c(both$mean_mari[both$Mari_sig=='yes'], both$mean_vib[both$Vib_sig=='yes'], both$mean_CC[both$CC_sig=='yes']) %>% 
    na.omit()
bothWAlone <- both$mean_rpomAlone


transport <- subset(metabolism, Metabolism=='Transport')
transportfold <- c(transport$Mari_fold_diff[transport$Mari_sig=='yes'], transport$Vib_fold_diff[transport$Vib_sig=='yes'], transport$CC_fold_diff[transport$CC_sig=='yes']) %>% 
    na.omit()
transportW <- c(transport$mean_mari[transport$Mari_sig=='yes'], transport$mean_vib[transport$Vib_sig=='yes'], transport$mean_CC[transport$CC_sig=='yes']) %>% 
    na.omit()
transportWAlone <- transport$mean_rpomAlone

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


#pdf('metabolic_FoldDifferences_20220603.pdf', useDingbats=FALSE)


plot(c(0.8, 1.1), c(-8,8), type='n', xlim=c(0,0.8), xlab='', ylab='Log2_FC', xaxt='n', bty='n', las=2)
axis(1, at=c(0.1, 0.3, 0.5, 0.7), labels=c('Tport', 'Cata', 'Amp', 'Ana'))
dotPlotter(transportfold, transportTest, 0.1)
dotPlotter(catabolicfold, catabolicTest, 0.3)
dotPlotter(bothfold, bothTest, 0.5)
dotPlotter(anabolicfold, anabolicTest, 0.7)



#dev.off()
