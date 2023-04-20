#!/bin/bash

## make a GFF3 reference file for the K12-MG1655-NC_000913.gb reference genome,
## for downstream analysis with the copy-number-analysis.R script.
gdtools APPLY -o ../results/nine-day-GFP-barcode-expt-genome-analysis/LCA.gff3 -f GFF3  -r ../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb ../results/nine-day-GFP-barcode-expt-genome-analysis/LCA.gd
