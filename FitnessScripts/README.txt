

These are scripts used to calculate the mean fitness for each experimental treatment. Inupt for these scripts are the normalized_terminusRM_spo.csv and the DSS3_annotations_April2020.csv

fitnessCalculation_20220621_Rp_V.R= used to calculate W for treatment Rp+M Rp+V
fitnessCalculation_20220621_Rp_M.R= used to calculate W for treatment Rp+M
fitnessCalculation_20220621_Rp_V_M.R= used to calculate W for treatment Rp+V+M

These scripts read in insertion count data, sum them per gene, then calculates mean relative fitness (W). Each script is written to be used to calculate mean relative fitness for the respective treatment. Returns both the mean relative fitness as well as the relative fitness values for each replicate to be used for the randomization resampling