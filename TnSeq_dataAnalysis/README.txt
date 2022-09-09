This folder contains all scripts used for data analysis/for producing the figures seen in the paper.

20220818_project_Summary.Rmd is an rmarkdown file that was used for creating the main figures seen in the text. This script accesses files in 'Data' folder. A version of this markdown file has been 'knit' to 20220818_project_Summary.html to show the output of this code, and should be reproducible.




Significant fitness phenotypes were combined with COG categories in the script COGanalysis_20220718.R by accessing the 'COGs' folder to create the file significantCOGfitness, which was combined with all significant data to create significantFitness_Metabolism_Grouped. The COGanalysis_20220718.R script was also used to generate supplemental figure 3 using the significantFitness_Metabolism_Grouped file. 


The 'Data' folder contains the raw values from growth curves (flow cytometry, CFUs, OD measurements), and bootstrap estimates of mean relative fitness. This folder also contains 'readStats_v2.csv' which details the number of reads in each treatment, as well as all_samples_tn_stats.txt, which contains the commands and outputs from TPP and TRANSIT. See the 'sequence_reads' release to download the raw fastq files used.
