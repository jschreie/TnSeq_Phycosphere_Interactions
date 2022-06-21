

These are scripts used to calculate the mean fitness for each experimental treatment. Inupt for these scripts are the normalized_terminusRM_spo.csv and the DSS3_annotations_April2020.csv


rpomAloneBootstrap = Rp
vibrioBootstrap = Rp+V
marivivensBootstrap = Rp+M
completeCommBootstrap = Rp+V+M

These scripts read in insertion count data, sum them per gene, then calculates mean relative fitness (W). Each script is written to be used to calculate mean relative fitness for the respective treatment.