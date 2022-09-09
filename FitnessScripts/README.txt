

These are scripts used to calculate the mean fitness for each experimental treatment. Inupt for these scripts are the 20211004_normalized_terminusRM_spo.csv and the DSS3_annotations_April2020.csv

fitnessCalculation_20220621_Rp_V.R= used to calculate W for treatment Rp+M Rp+V
fitnessCalculation_20220621_Rp_M.R= used to calculate W for treatment Rp+M
fitnessCalculation_20220621_Rp_V_M.R= used to calculate W for treatment Rp+V+M

These scripts read in insertion count data, sum them per gene, then calculates mean relative fitness (W). Each script is written to be used to calculate mean relative fitness for the respective treatment. Returns both the mean relative fitness as well as the relative fitness values for each replicate to be used for the randomization resampling

fitnessAnalysis script combines the output from fitness calculations into files that can be used for input to the randomization procedure
randomization_FitnessDiffs scripts are used to perform 10000 permutations of a randomization test for statistical significance.



Data folder contains the output of the fitness calcuation (_raw), mean fitness values, results from randomization (_rando) tests per treatment, finally combining all data into meanFitness_pval and significantW_randomization
