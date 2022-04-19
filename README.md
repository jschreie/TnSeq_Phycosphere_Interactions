# TnSeq_Phycosphere_Interactions

Sequence files (.fastq) are downloadable through the release titled 'Transposon insertion fastq files'. These were processed using TPP as part of TRANSIT. See the file 'all_samples_tn_stats.txt' in the 'TnSeq_dataAnalysis/Data' folder of this repository to see a summary of commands used and output of TPP.

Transposon insertion counts (as a 20200218_combined_all.wig), normalized counts (20211004_normalized_terminusRM_spo.wig), as well as the scripts used to normalize the transposon insertion counts and focus on the central 90% of the gene are downloadable through the release titled 'Transposon insertion counts - raw and normalized'


FitnessAnalysis contains scripts and files for  calculating bootstrap estimates of fitness from a normalized .wig file. 

TnSeq_dataAnalysis contains scripts and files for analyzing and plotting fitness data that has been output from the FitnessAnalysis (these data are combined and are present in this folder as meanFitness_20220316.csv, which is then subset into statistically significant mutants in the file significant_genes_grouped_Metabolism_filtered_20220413.csv)
