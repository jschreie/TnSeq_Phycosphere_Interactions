# TnSeq_Phycosphere_Interactions

Sequence files (.fastq) are available from NCBI SRA PRJNA910220. These read files were processed using TPP as part of TRANSIT using the script tpp_loop.sh, which uses the 'jeremy_sample_names.txt' file and outputs statistics in the 'all_samples_tn_stats.txt'. DSS3_samples_metadata.xlsx ties the name of the fastq files to the treatment., DSS3mega.fna is a nucleotide file containing _R.pomeroyi_ chromosome and megaplasmid, used for aligning insertion reads to the genome. The output of TPP was merged using the combine_wig function in TRANSIT.

Transposon insertion counts (as a 20200218_combined_all.wig), normalized counts (20211004_normalized_terminusRM_spo.wig), as well as the scripts used to normalize the transposon insertion counts and focus on the central 90% of the gene are downloadable through the release titled 'Transposon insertion counts - raw and normalized'

FitnessAnalysis contains scripts and files for calculating fitness from a normalized transposon insertion count data. 

TnSeq_dataAnalysis contains scripts and files for analyzing and plotting fitness data that has been output from the FitnessAnalysis (these data are combined and are present in this folder as meanFitness_20220701.csv, which is then subset into statistically significant mutants in the file SuppTable5_significantFitness_Metabolism_Grouped_Organized_20220815.xlsx)
