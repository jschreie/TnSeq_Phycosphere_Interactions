

These are scripts used to calculate the bootstrap distribution of fitness for each experimental treatment. Inupt for these scripts are the normalized_terminusRM_spo.csv and the DSS3_annotations_April2020.csv


rpomAloneBootstrap = Rp
vibrioBootstrap = Rp+V
marivivensBootstrap = Rp+M
completeCommBootstrap = Rp+V+M

These scripts read in insertion count data, sum them per gene, then calculates mean relative fitness (W) by resampling 10000 means of 4 replicates. Each script is written to be used to calculate bootstrap estimates of mean relative fitness for the respective treatment.