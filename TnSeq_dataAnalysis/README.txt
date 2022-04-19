This folder contains all scripts used for data analysis/for producing the figures seen in the paper.

fitnessAnalysis_20220413.R takes the bootstrap estimates of fitness in the 'Data' folder and combines them to create the file 'meanFitness_20220316.csv' and a file of significant fitness. This significant fitness file was then combined with COG categories in the script COGanalysis_20220413.R by accessing the 'COGs' folder to create the file 'significant_genes_grouped_Metabolism_filtered_20220413.csv', which has been added to this file for ease of code usage. The COGanalysis_20220413.R script was also used to generate supplemental figure 2. 


20220419_project_Summary.Rmd is an rmarkdown file that was used for creating the main figures seen in the text. This script accesses files in 'Data' folder, as well as meanFitness_20220316.csv, significant_genes_grouped_Metabolism_filtered_20220413.csv, and Table1_20220328.csv. A version of this markdown file has been 'knit' to 20220419_project_Summary.html to show the output of this code.


The 'Data' folder contains the raw values from growth curves (flow cytometry, CFUs, OD measurements), and bootstrap estimates of fitness. This folder also contains 'readStats_v2.csv' which details the number of reads in each treatment, as well as all_samples_tn_stats.txt, which contains the commands and outputs from TPP and TRANSIT. See the 'sequence_reads' release to download the raw fastq files used.
