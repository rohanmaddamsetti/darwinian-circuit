#!/usr/bin/env python

""" get-DH5a-expt-transposon-coverage.py by Rohan Maddamsetti.

This script calculates the coverage for the actual transposons
in the evolved clones from my nine-day DH5a evolution experiment.

Usage: python get-DH5a-expt-transposon-coverage.py > ../results/DH5a-expt-genome-analysis/transposon-coverage.csv

"""

import os

## Global constants from B30-miniTn5-TetA.gb. (strong tetA promoter)
B30_TRANSPOSON_COORDS = (2990, 4367)

## Global constants from B59-miniTn5-TetA.gb.
## (strong tetA promoter, no Tn5 in these strains).
B59_TRANSPOSON_COORDS = (1459, 2836)

## Global constants from B20-miniTn5-TetA.gb.
## (medium tetA promoter).
B20_TRANSPOSON_COORDS = (2990, 4367)


def get_transposon_coverage(transposon_coverage_f):
    """ This function can only handle B30, B20, and B59 coverage files. """
    if "B30" in transposon_coverage_f:
        Tn_start, Tn_end = B30_TRANSPOSON_COORDS
    elif "B59" in transposon_coverage_f:
        Tn_start, Tn_end = B59_TRANSPOSON_COORDS
    elif "B20" in transposon_coverage_f:
        Tn_start, Tn_end = B20_TRANSPOSON_COORDS
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


def get_transposon_coverage_file(coverage_dir, transposon_coverage_prefix_list = ["B30","B59","B20"]):
    for prefix in transposon_coverage_prefix_list:
        matched_coverage_file_list = [x for x in os.listdir(coverage_dir) if x.startswith(prefix)]
        if len(matched_coverage_file_list):
            return os.path.join(coverage_dir, matched_coverage_file_list[0])
    return None


def get_transposon_coverage_for_sample(breseq_outpath, transposon_coverage_prefix="B31"):
    my_sample = os.path.basename(breseq_outpath)
    coverage_dir = os.path.join(breseq_outpath, "08_mutation_identification")
    transposon_coverage_f = get_transposon_coverage_file(coverage_dir)
    my_transposon_coverage = get_transposon_coverage(transposon_coverage_f)
    return (my_sample, my_transposon_coverage)


def main():

    breseq_mixedpop_dir = "../results/DH5a-expt-genome-analysis/mixed-pops"
    mixedpop_paths = [os.path.join(breseq_mixedpop_dir, x) for x in os.listdir(breseq_mixedpop_dir) if x.startswith("RM")]

    sample_vec = []
    transposon_vec = []
    transposon_coverage_vec = []
            
    ## walk through the file structure for the evolved samples.
    for p in sorted(mixedpop_paths):
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
