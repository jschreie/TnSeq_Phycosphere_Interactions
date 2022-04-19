# TnSeq_Phycosphere_Interactions

Sequence files (.fastq) are in the zip file titled 'MergedFASTQ' in the RawData folder. These were processed using TPP as part of TRANSIT using the script. This fodler also contains the 'tpp_loop.sh' script which uses the 'jeremy_sample_names.txt'. Output/files processed using TPP are in the folder fastq_tpp_processed.zip. See the file 'all_samples_tn_stats.txt' to see a summary of commands used and output of TPP. DSS3_samples_metadata.xlsx ties the name of the fastq files to the treatment., DSS3mega.fna is a nucleotide file containing _R.pomeroyi_ chromosome and megaplasmid, used for aligning insertion reads to the genome.

Transposon insertion counts (as a 20200218_combined_all.wig), normalized counts (20211004_normalized_terminusRM_spo.wig), as well as the scripts used to normalize the transposon insertion counts and focus on the central 90% of the gene are downloadable through the release titled 'Transposon insertion counts - raw and normalized'


FitnessAnalysis contains scripts and files for  calculating bootstrap estimates of fitness from a normalized .wig file. 

TnSeq_dataAnalysis contains scripts and files for analyzing and plotting fitness data that has been output from the FitnessAnalysis (these data are combined and are present in this folder as meanFitness_20220316.csv, which is then subset into statistically significant mutants in the file significant_genes_grouped_Metabolism_filtered_20220413.csv)
