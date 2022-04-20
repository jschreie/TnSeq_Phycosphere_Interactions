#!/bin/bash

for i in $(cat jeremy_sample_names.txt);
python3.6 /home/jeremy/Desktop/transit-master/src/tpp.py -bwa /usr/local/bin/bwa -bwa-alg 'aln' -ref /home/jeremy/Desktop/Sequencing_Info/DSS3mega.fna -reads1 /home/jeremy/Desktop/TnSeqTransit/${i}_merged.fastq -output /home/jeremy/Desktop/TnSeqTransit/$i -protocol Tn5 -primer AGATGTGTATAAGAGACAG -primer-start-window 0,60 -mismatches 3
done


cat */*.tn_stats > all_samples_tn_stats.txt
