#!/usr/bin/env python

"""
subtract_LCA_from_gdiffs.py by Rohan Maddamsetti.

1) get the intersection of all annotated gdiffs to infer the LCA.
2) subtract the LCA_intersection.gd file from all evolved *.gd files
to make *_subtracted_LCA.gd files.

## the output files are then turned into a CSV file by the script process_subtracted_gdiffs.py.

"""

import os.path
import subprocess
from itertools import chain
from os import listdir

RESULTSDIR = "../results/nine-day-GFP-barcode-expt-genome-analysis"

# full method with doctests
def interleave(*iters):
    """
    From https://gist.github.com/smdabdoub/5213405, but izip is now zip in Python 3.
    Given two or more iterables, return a list containing 
    the elements of the input list interleaved.
    
    >>> x = [1, 2, 3, 4]
    >>> y = ('a', 'b', 'c', 'd')
    >>> interleave(x, x)
    [1, 1, 2, 2, 3, 3, 4, 4]
    >>> interleave(x, y, x)
    [1, 'a', 1, 2, 'b', 2, 3, 'c', 3, 4, 'd', 4]
    
    On a list of lists:
    >>> interleave(*[x, x])
    [1, 1, 2, 2, 3, 3, 4, 4]
    
    Note that inputs of different lengths will cause the 
    result to be truncated at the length of the shortest iterable.
    >>> z = [9, 8, 7]
    >>> interleave(x, z)
    [1, 9, 2, 8, 3, 7]
    
    On single iterable, or nothing:
    >>> interleave(x)
    [1, 2, 3, 4]
    >>> interleave()
    []
    """
    return list(chain(*zip(*iters)))


def write_LCA_file(genome_samples, LCA_outputpath, resultsdir=RESULTSDIR):
    gd_paths = [os.path.join(resultsdir, x, "output/evidence/annotated.gd") for x in genome_samples]
    ## now let's run gdtools INTERSECT on these gd files.
    gdtools_args = ["gdtools", "INTERSECT", "-o", LCA_outputpath] + gd_paths
    if not os.path.exists(LCA_outputpath):
        subprocess.run(gdtools_args)
    return None


def annotate_GDfile(gd_path, refgenomepath_list, gd_format="HTML"):
    my_dirname, my_basename = os.path.split(gd_path)
    just_basename, my_extension = my_basename.split(".")

    if gd_format == "HTML":
        my_new_basename = just_basename + ".html"
    else:
        my_new_basename = just_basename + ".gd"
    my_annotated_basename = "annotated_" + my_new_basename
    annotated_outputfile = os.path.join(my_dirname, my_annotated_basename)

    ref_flag_list = ["-r" for x in refgenomepath_list]
    refgenome_args = interleave(ref_flag_list, refgenomepath_list)

    gdtools_args = ["gdtools", "ANNOTATE"] + refgenome_args + ["-f", gd_format, "-o", annotated_outputfile, gd_path]
    if not os.path.exists(annotated_outputfile):
        subprocess.run(gdtools_args)
    return None


def subtract_LCA(genome_sample_id, LCA_path, resultsdir=RESULTSDIR):
    gd_path = os.path.join(resultsdir, genome_sample_id, "output/evidence/annotated.gd")
    output_filename = "subtracted_"  + genome_sample_id + ".gd"
    subtracted_output_file = os.path.join(resultsdir, output_filename)
    ## gdtools SUBTRACT [-o output.gd] input.gd subtract1.gd [subtract2.gd ...]
    gdtools_args = ["gdtools", "SUBTRACT", "-o", subtracted_output_file, gd_path, LCA_path]
    if not os.path.exists(subtracted_output_file):
        subprocess.run(gdtools_args)
    return None


def main():

    ## 1) get the intersection of all gdiffs to infer the LCA.
    ## let's first make a list of all the gdiffs.
    ## The samples range from RM7-140-31 to RM7-140-60, inclusive.
    all_samples = ["RM7-140-" + str(x) for x in range(31,60+1)]
    LCA_outputpath = os.path.join(RESULTSDIR, "LCA.gd")
    write_LCA_file(all_samples, LCA_outputpath)
    
    ## let's check the LCA for each treatment as a control.
    no_plasmid_samples = ["RM7-140-" + str(x) for x in range(31,40+1)]
    no_plasmid_LCA_outputpath = os.path.join(RESULTSDIR, "no_plasmid_LCA.gd")
    write_LCA_file(no_plasmid_samples, no_plasmid_LCA_outputpath)
    
    p15A_samples = ["RM7-140-" + str(x) for x in range(41,50+1)]
    p15A_LCA_outputpath = os.path.join(RESULTSDIR, "p15A_LCA.gd")
    write_LCA_file(p15A_samples, p15A_LCA_outputpath)
    
    pUC_samples = ["RM7-140-" + str(x) for x in range(51,60+1)]
    pUC_LCA_outputpath = os.path.join(RESULTSDIR, "pUC_LCA.gd")
    write_LCA_file(pUC_samples, pUC_LCA_outputpath)

    ## let's annotate these gd files, and take a look.
    no_plasmid_refgenomes = ["../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb", "../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb"]
    p15A_refgenomes = ["../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb", "../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb", "../data/genome-sequencing/reference-genome/A31-p15A.gb"]
    pUC_refgenomes = ["../data/genome-sequencing/reference-genome/K12-MG1655-NC_000913.gb", "../data/genome-sequencing/reference-genome/B31-N20-miniTn5-TetA-GFP-barcode.gb", "../data/genome-sequencing/reference-genome/A18-pUC.gb"]

    ## make HTML output.
    annotate_GDfile(LCA_outputpath, no_plasmid_refgenomes)
    annotate_GDfile(no_plasmid_LCA_outputpath, no_plasmid_refgenomes)
    annotate_GDfile(p15A_LCA_outputpath, p15A_refgenomes)
    annotate_GDfile(pUC_LCA_outputpath, pUC_refgenomes)

    ## now, subtract the treatment-specific LCA_intersection.gd file from the 
    ## the set of evolved annotated.gd files in that treatment, and write to RESULTSDIR.
    for my_sample in no_plasmid_samples:
        subtract_LCA(my_sample, no_plasmid_LCA_outputpath)

    for my_sample in p15A_samples:
        subtract_LCA(my_sample, p15A_LCA_outputpath)

    for my_sample in pUC_samples:
        subtract_LCA(my_sample, pUC_LCA_outputpath)

            
main()
