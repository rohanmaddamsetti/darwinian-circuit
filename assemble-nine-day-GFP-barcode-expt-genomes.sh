#!/usr/bin/env bash

## assemble-nine-day-GFP-barcode-expt-genomes.sh by Rohan Maddamsetti.
## This shell script does a bunch of genome assembly tasks using breseq and gdtools.
## COMMENT OUT all tasks that don't need to be done.

#########################################################################
## Assemble the ancestral clones using the K12-MG1655-NC_000913.gb reference genome.
## sample names are not contiguous because I threw away clones with bad barcodes.

## The no plasmid ancestor clones are RM7.107.3,4,5,6,7,8,9,10,11,13.
#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_3_S550_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_3_S550_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-4 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_4_S551_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_4_S551_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-5 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_5_S552_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_5_S552_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-6 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_6_S553_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_6_S553_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-7 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_7_S554_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_7_S554_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-8 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_8_S555_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_8_S555_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-9 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_9_S556_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_9_S556_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-10 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_10_S557_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_10_S557_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-11 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_11_S558_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_11_S558_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-13 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_13_S559_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_13_S559_R2_001.fastq.gz"


## The p15A ancestor clones are RM7.106.3,5,6,7,8,10,11,12,13,15.
#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A31-p15A.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_3_S560_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_3_S560_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-5 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A31-p15A.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_5_S561_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_5_S561_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-6 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A31-p15A.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_6_S562_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_6_S562_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-7 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A31-p15A.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_7_S563_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_7_S563_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-8 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A31-p15A.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_8_S564_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_8_S564_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-10 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A31-p15A.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_10_S565_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_10_S565_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-11 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A31-p15A.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_11_S566_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_11_S566_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-12 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A31-p15A.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_12_S567_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_12_S567_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-13 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A31-p15A.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_13_S568_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_13_S568_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-15 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A31-p15A.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_15_S569_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_15_S569_R2_001.fastq.gz"


## The pUC ancestor clones are RM7.107.38,39,41,42,44,45,46,48,49,56.
#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-38 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A18-pUC.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_38_S570_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_38_S570_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-39 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A18-pUC.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_39_S571_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_39_S571_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-41 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A18-pUC.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_41_S572_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_41_S572_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-42 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A18-pUC.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_42_S573_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_42_S573_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-44 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A18-pUC.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_44_S574_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_44_S574_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-45 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A18-pUC.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_45_S575_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_45_S575_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-46 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A18-pUC.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_46_S576_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_46_S576_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-48 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A18-pUC.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_48_S577_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_48_S577_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-49 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A18-pUC.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_49_S578_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_49_S578_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-56 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A18-pUC.gb ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_56_S579_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_56_S579_R2_001.fastq.gz"

#########################################################################
## Generate experiment-specific reference sequences by applying mutations to the K-12 reference sequence.

## make references for the ancestral no plasmid clones.
#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-3.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-3/output/output.gd

#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-4.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-4/output/output.gd

#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-5.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-5/output/output.gd

#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-6.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-6/output/output.gd

#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-7.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-7/output/output.gd

#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-8.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-8/output/output.gd

#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-9.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-9/output/output.gd

#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-10.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-10/output/output.gd

#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-11.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-11/output/output.gd

#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-13.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-13/output/output.gd

## make references for the ancestral p15A clones.
#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-3.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A31-p15A.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-3/output/output.gd

#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-5.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A31-p15A.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-5/output/output.gd

#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-6.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A31-p15A.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-6/output/output.gd

#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-7.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A31-p15A.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-7/output/output.gd

#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-8.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A31-p15A.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-8/output/output.gd

#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-10.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A31-p15A.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-10/output/output.gd

#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-11.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A31-p15A.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-11/output/output.gd

#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-12.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A31-p15A.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-12/output/output.gd

#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-13.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A31-p15A.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-13/output/output.gd

#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-15.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A31-p15A.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-15/output/output.gd

## make references for the ancestral pUC clones.
#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-38.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A18-pUC.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-38/output/output.gd

#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-39.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A18-pUC.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-39/output/output.gd

#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-41.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A18-pUC.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-41/output/output.gd

#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-42.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A18-pUC.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-42/output/output.gd

#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-44.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A18-pUC.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-44/output/output.gd

#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-45.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A18-pUC.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-45/output/output.gd

#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-46.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A18-pUC.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-46/output/output.gd

#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-48.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A18-pUC.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-48/output/output.gd

#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-49.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A18-pUC.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-49/output/output.gd

#gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-56.gff3 -f GFF3 -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb -r ../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb -r ../data/genome-sequencing/reference-genome/A18-pUC.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-56/output/output.gd

#########################################################################
## Now remap reads to the ancestral genomes for quality control.

## ancestral no plasmid clones.
#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-107-3 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-3.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_3_S550_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_3_S550_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-107-4 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-4.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_4_S551_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_4_S551_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-107-5 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-5.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_5_S552_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_5_S552_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-107-6 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-6.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_6_S553_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_6_S553_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-107-7 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-7.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_7_S554_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_7_S554_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-107-8 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-8.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_8_S555_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_8_S555_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-107-9 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-9.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_9_S556_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_9_S556_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-107-10 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-10.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_10_S557_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_10_S557_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-107-11 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-11.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_11_S558_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_11_S558_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-107-13 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-13.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_13_S559_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_13_S559_R2_001.fastq.gz"


## ancestral p15A clones.
#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-106-3 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-3.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_3_S560_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_3_S560_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-106-5 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-5.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_5_S561_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_5_S561_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-106-6 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-6.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_6_S562_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_6_S562_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-106-7 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-7.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_7_S563_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_7_S563_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-106-8 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-8.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_8_S564_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_8_S564_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-106-10 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-10.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_10_S565_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_10_S565_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-106-11 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-11.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_11_S566_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_11_S566_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-106-12 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-12.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_12_S567_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_12_S567_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-106-13 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-13.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_13_S568_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_13_S568_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-106-15 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-15.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_15_S569_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_106_15_S569_R2_001.fastq.gz"


## ancestral pUC clones.
#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-107-38 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-38.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_38_S570_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_38_S570_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-107-39 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-39.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_39_S571_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_39_S571_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-107-41 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-41.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_41_S572_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_41_S572_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-107-42 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-42.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_42_S573_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_42_S573_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-107-44 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-44.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_44_S574_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_44_S574_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-107-45 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-45.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_45_S575_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_45_S575_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-107-46 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-46.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_46_S576_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_46_S576_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-107-48 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-48.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_48_S577_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_48_S577_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-107-49 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-49.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_49_S578_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_49_S578_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/remapped-RM7-107-56 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-56.gff3  ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_56_S579_R1_001.fastq.gz ../data/genome-sequencing/SeqCenter18042023_ancestral-tetA-GFP-barcode-clones/RM7_107_56_S579_R2_001.fastq.gz"

#########################################################################
## Assemble evolved genomes using barcode-specific reference genomes.

#####
## Assemble K12 + B31 no plasmid clones 1-10, using barcode-specific references.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-31 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-3.gff3 ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S17_S17_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S17_S17_L001_R2_001.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-32 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-4.gff3 ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S18_S18_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S18_S18_L001_R2_001.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-33 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-5.gff3 ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S19_S19_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S19_S19_L001_R2_001.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-34 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-6.gff3 ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S20_S20_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S20_S20_L001_R2_001.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-35 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-7.gff3 ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S21_S21_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S21_S21_L001_R2_001.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-36 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-8.gff3 ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S22_S22_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S22_S22_L001_R2_001.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-37 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-9.gff3 ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S23_S23_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S23_S23_L001_R2_001.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-38 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-10.gff3 ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S24_S24_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S24_S24_L001_R2_001.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-39 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-11.gff3 ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S25_S25_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S25_S25_L001_R2_001.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-40 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-13.gff3 ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S26_S26_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S26_S26_L001_R2_001.fastq.gz"

## Assemble K12 + B31 + p15A clones 1-6, using barcode-specific references.
##These were included in run 8336 for the Duke Sequencing Voucher.

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-41 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-3.gff3 ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S27_S27_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S27_S27_L001_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-42 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-5.gff3 ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S28_S28_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S28_S28_L001_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-43 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-6.gff3 ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S29_S29_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S29_S29_L001_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-44 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-7.gff3 ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S30_S30_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S30_S30_L001_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-45 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-8.gff3 ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S31_S31_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S31_S31_L001_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-46 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-10.gff3 ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S32_S32_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8336_tetA-GFP-barcode-clones/8336-S32_S32_L001_R2_001.fastq.gz"

## Assemble K12 + B31 + p15A clones 7-10, using K-12 as reference.
## These were included in run 8337 for the Duke Sequencing Voucher.

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-47 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-11.gff3 ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S115_S115_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S115_S115_L001_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-48 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-12.gff3 ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S116_S116_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S116_S116_L001_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-49 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-13.gff3 ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S117_S117_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S117_S117_L001_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-50 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-106-15.gff3 ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S118_S118_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S118_S118_L001_R2_001.fastq.gz"

## Assemble K12 + B31 + pUC clones using K-12 as reference.
## These were included in run 8337 for the Duke Sequencing Voucher.

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-51 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-38.gff3 ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S119_S119_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S119_S119_L001_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-52 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-39.gff3 ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S120_S120_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S120_S120_L001_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-53 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-41.gff3 ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S121_S121_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S121_S121_L001_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-54 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-42.gff3 ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S122_S122_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S122_S122_L001_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-55 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-44.gff3 ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S123_S123_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S123_S123_L001_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-56 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-45.gff3 ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S124_S124_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S124_S124_L001_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-57 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-46.gff3 ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S125_S125_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S125_S125_L001_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-58 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-48.gff3 ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S126_S126_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S126_S126_L001_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-59 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-49.gff3 ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S127_S127_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S127_S127_L001_R2_001.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -o ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-140-60 -r ../results/nine-day-GFP-barcode-expt-genome-analysis/RM7-107-56.gff3 ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S128_S128_L001_R1_001.fastq.gz ../data/genome-sequencing/Maddamsetti_8337_tetA-GFP-barcode-clones/8337-S128_S128_L001_R2_001.fastq.gz"
