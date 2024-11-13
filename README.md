# darwinian-circuit

## Software requirements
## Python 3.10+

## COMPUTATIONAL PROTOCOL FOR TRANSPOSON-PLASMID ANALYSIS

Make a top-level directory with three directories inside, named "data", "results", and "src".  
Now copy all source code files in this repository into "src".

Now, download prokaryotes.txt into ../data/:  

wget https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt  

Then, filter the prokaryotes.txt genome data for those that have complete genomes,
and replace "GCA" with "GCF" throughout this file, so that RefSeq data and not Genbank data
is accessed in all downstream steps:  

python filter-genome-reports.py > ../results/complete-prokaryotes-with-plasmids.txt  

Then, fetch genome annotation for each row in complete-prokaryotes-with-plasmids.txt,
fetch the protein-coding genes for all chromosomes and plasmids for
each row in best-prokaryotes.txt,

These steps can be done at the same time on the Duke Compute Cluster (DCC).
And make sure these scripts are called from the src directory.
fetch-gbk-annotation runs for several hours.  

First, make some output directories.  
mkdir ../results/gbk-annotation/

Then:  
sbatch --mem=16G -t 36:00:00 --wrap="python fetch-gbk-annotation.py"  