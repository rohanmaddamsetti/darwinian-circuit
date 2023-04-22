#!/usr/bin/env python

"""
process_gdiffs.py by Rohan Maddamsetti.

I take the evolved genomes, and write out evolved-mutations.csv for
downstream analysis in R.

Usage: python process_gdiffs.py

"""

import os.path
## add genomediff parser to path.
import sys
sys.path.append("genomediff-python")
import genomediff


def get_plasmid(sample):
    ## get the sample ID, and subtract 30 to get a sample between 1 and 30.
    pop_num = int(sample.split("RM7-140-")[-1]) - 30
    if pop_num <= 10:
        plasmid = "None"
    elif pop_num <= 20:
        plasmid = "p15A"
    elif pop_num <= 30:
        plasmid = "pUC"
    else:
        raise AssertionError("ERROR:  out-of-bounds sample ID!")
    return plasmid


def get_pop(sample):
    ## get the sample ID, and subtract 30 to get a sample between 1 and 30.
    pop_num = int(sample.split("RM7-140-")[-1]) - 30
    pop = (pop_num % 10) + 1
    ## then convert to a text string.
    return str(pop)
    

def write_evolved_mutations(gdiff_paths, outf):
    outfh = open(outf, "w")
    outfh.write("Sample,Plasmid,Population,Mutation,Mutation_Category,Gene,Position,Allele\n")
    
    for gdiff in gdiff_paths:
        infh = open(gdiff, 'r', encoding='utf-8')
        ## hacky solution to get the sample name from the full path.
        sample = os.path.dirname(gdiff).split("/output/evidence")[0].split('/')[-1]
        plasmid = get_plasmid(sample)
        population = get_pop(sample)
        gd = genomediff.GenomeDiff.read(infh)
        muts = gd.mutations
        muttype = ""
        allele = ""
        for rec in muts:
            pos = str(rec.attributes['position'])
            mut_category = rec.attributes['mutation_category']
            gene = rec.attributes['gene_name']
            ## handle transpositions of the miniTn5 transposon.
            if "repeat_name" in rec.attributes and rec.attributes["repeat_name"] == "miniTn5":
                gene = "tetA-Tn5-" + gene

            if 'new_seq' in rec.attributes:
                allele = rec.attributes['ref_seq'] + "->" + rec.attributes['new_seq']
            else:
                allele = rec.type
            if rec.type == 'SNP':
                muttype = rec.attributes['snp_type']
            else:
                muttype = rec.type
            ## skip barcode "mutations"
            if allele.startswith("NNNNNNNNNN"):
                continue
            mutation_row_data = ','.join([sample, plasmid, population, muttype, mut_category, gene, pos, allele])
            outfh.write(mutation_row_data + "\n")
    outfh.close()
    return None


def main():
    srcdir = os.path.dirname(os.path.realpath(__file__))
    assert srcdir.endswith("src")
    projdir = os.path.dirname(srcdir)
    assert projdir.endswith("darwinian-circuit")
    genome_results_dir = os.path.join(projdir,"results","nine-day-GFP-barcode-expt-genome-analysis")

    evolved_clone_ids = ["RM7-140-" + str(x) for x in range(31,61)]
    evolved_genomediffs = [x + "/output/evidence/annotated.gd" for x in evolved_clone_ids]
    gdiff_paths = [os.path.join(genome_results_dir, x) for x in evolved_genomediffs]
    
    ''' tabulate the mutations in the evolved populations.'''
    mut_table_outf = os.path.join(genome_results_dir, "evolved_mutations.csv")
    write_evolved_mutations(gdiff_paths, mut_table_outf)

    
main()
    
