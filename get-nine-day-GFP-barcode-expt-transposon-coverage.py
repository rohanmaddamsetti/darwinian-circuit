#!/usr/bin/env python

""" get-nine-day-GFP-barcode-expt-transposon-coverage.py by Rohan Maddamsetti.

This script calculates the coverage for the actual B31-N20 transposons
in the evolved clones from my nine-day tetA-GFP-barcode evolution experiment

Usage: python get-nine-day-GFP-barcode-expt-transposon-coverage.py > ../results/nine-day-GFP-barcode-expt-genome-analysis/transposon-coverage.csv

"""

import os

## Global constants from B31-N20-miniTn5-TetA-GFP-barcode.gb.
B31_N20_TRANSPOSON_COORDS = (465, 2567)


def get_B31_N20_transposon_coverage(transposon_coverage_f):
    """ This function can only handle B31_N20 coverage files. """
    if "B31-N20" in transposon_coverage_f:
        Tn_start, Tn_end = B31_N20_TRANSPOSON_COORDS
    else:
        raise AssertionError("ERROR: Unknown transposon in coverage file name!")

    total_coverage = 0.0
    positions_examined = 0
    with open(transposon_coverage_f, "r") as transposon_coverage_fh:
        ## the header is line number 0! so line 1 is the first line with data.
        for i, line in enumerate(transposon_coverage_fh):
            if (i < Tn_start): continue
            if (i >= Tn_end): break ## don't count the last position of the transposon. some coverage oddity?
            positions_examined += 1
            fields = line.split("\t") ## tab-delimited data.
            top_coverage = float(fields[0])
            bottom_coverage = float(fields[1])
            total_position_coverage = top_coverage + bottom_coverage
            total_coverage += total_position_coverage

        ## round to 2 decimal places.
        my_transposon_coverage = str(round(float(total_coverage)/float(positions_examined),2))
    return my_transposon_coverage


def get_transposon_coverage_for_sample(breseq_outpath, transposon_coverage_prefix="B31-N20"):
    my_sample = os.path.basename(breseq_outpath)
    coverage_dir = os.path.join(breseq_outpath, "08_mutation_identification")
    transposon_coverage_f = [os.path.join(coverage_dir, x) for x in os.listdir(coverage_dir) if (x.startswith(transposon_coverage_prefix))][0]
    my_transposon_coverage = get_B31_N20_transposon_coverage(transposon_coverage_f)
    return (my_sample, my_transposon_coverage)


def main():

    breseq_dir = "../results/nine-day-GFP-barcode-expt-genome-analysis"

    ## note: this includes both evolved clones and ancestral clones.
    clone_paths = [os.path.join(breseq_dir, x) for x in os.listdir(breseq_dir) if x.startswith("RM7") and not x.endswith("gff3")]

    sample_vec = []
    transposon_vec = []
    transposon_coverage_vec = []
            
    ## walk through the file structure for the evolved samples.
    for p in sorted(clone_paths):
        my_sample, my_transposon_coverage = get_transposon_coverage_for_sample(p)
        sample_vec.append(my_sample)
        transposon_coverage_vec.append(my_transposon_coverage)

    ## now print the data to file.
    rownum = len(sample_vec)
    assert len(transposon_coverage_vec) == rownum

    print("Sample,TransposonCoverage")
    for i in range(rownum):
        myrow = ",".join([sample_vec[i], transposon_coverage_vec[i]])
        print(myrow)

        
main()
