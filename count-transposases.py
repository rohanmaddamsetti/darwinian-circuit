#!/usr/bin/env python

"""
 count-transposases.py by Rohan Maddamsetti.

This script counts the number of transposases on plasmids versus chromosomes
in each genome in ../results/gbk-annotation.

IMPORTANT: filter out genomes with plasmids larger than chromosome.

proposed solution: for each genome, write relevant data per replicon into a list (buffer).
then, check that the buffer doesn't contain a plasmid larger than chromosome.
if the assertion passes, then write the buffer to file.

Usage: python count-transposases.py > ../results/transposase-counts.csv
"""

import os
import re
from Bio import SeqIO


def is_transposase(product_annotation, debugging=False):
    transposon_keywords = "^IS|transpos\S*|insertion|conjugate transposon|Transpos\S*|Tn[0-9]|tranposase|Tnp|Ins|ins"
     ## check if the gene is a transposase.
    is_transposase = True if re.search(transposon_keywords, product_annotation) else False
    if debugging and is_transposase:
        print(product_annotation)
    return is_transposase


def longest_replicon_is_chromosome(genome_buffer):
    is_longest_replicon_a_chromosome = False
    longest_replicon_type = "NA"
    longest_replicon_length = 0
    for line in genome_buffer:
        AnnotationAccession, SeqID, SeqType, SeqLength, CDS_count, CDS_length, transposase_count, transposase_length = line.split(",")
        if int(SeqLength) > longest_replicon_length:
            longest_replicon_length = int(SeqLength)
            longest_replicon_type = SeqType
    if longest_replicon_type == "chromosome":
        is_longest_replicon_a_chromosome = True
    return is_longest_replicon_a_chromosome


def main():
    gbk_annotation_path = "../results/gbk-annotation/"
    gbk_files = [x for x in os.listdir(gbk_annotation_path) if x.endswith(".gbff")]

    ## print the header of the file
    print("AnnotationAccession,SeqID,SeqType,SeqLength,CDS_count,CDS_length,transposase_count, transposase_length")
    for gbk in gbk_files:
        AnnotationAccession = gbk.split("_genomic.gbff")[0]
        gbk_path = os.path.join(gbk_annotation_path, gbk)

        genome_buffer = list() ## buffer output for this genome-- only print if the chromosome is larger than all the plasmids.
        
        replicon_records = SeqIO.parse(gbk_path, 'genbank')
        ## Iterate through each SeqRecord in the GenBank file
        for replicon in replicon_records:
            replicon_length = str(len(replicon.seq))
            cds_count = 0
            cds_length = 0
            transposase_count = 0
            transposase_length = 0

            SeqID = replicon.id
            if "plasmid" in replicon.description:
                SeqType = "plasmid"
            elif "complete" in replicon.description or "chromosome" in replicon.description:
                SeqType = "chromosome"
            else:
                SeqType = "unknown"
            
            ## Iterate through features and calculate CDS length
            for feature in replicon.features:
                my_seq_length = len(feature.location)
                if feature.type == 'CDS' and "translation" in feature.qualifiers:
                    cds_count += 1
                    cds_length += my_seq_length
                    ## now count transposases among CDS features.
                    if "product" in feature.qualifiers and is_transposase(feature.qualifiers["product"][0]):
                        transposase_count += 1
                        transposase_length += my_seq_length

            ## now add a line for this replicon to the genome buffer.
            genome_buffer.append(",".join([AnnotationAccession, replicon.id, SeqType, replicon_length,
                                           str(cds_count), str(cds_length),str(transposase_count),str(transposase_length)]))
            ## check this genome to make sure that there are no weird plasmid/chromosome annotation errors.
        if longest_replicon_is_chromosome(genome_buffer):
            for line in genome_buffer:
                print(line)
                
    return


if __name__ == "__main__":
    main()
