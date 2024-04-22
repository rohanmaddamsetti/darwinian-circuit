#!/usr/bin/env python

"""
high-PCN-plasmid-SV-pipeline.py by Rohan Maddamsetti.

For this pipeline to work, ncbi datasets and pysradb must be in the $PATH.

In a separate project, I calculated a large dataset of plasmid copy numbers.
I filtered those data for high confidence plasmids with copy numbers > 100.
These data are in ../data/high-PCN-plasmids.csv

This pipeline runs the following steps:
1) download long-read data for these genomes.
2) use minimap2 -L mode, with either ONT or PacBio parameters for minimap2
to make BAM alignments for Sniffles2.
3) Use Sniffles2 to call structural variants on these HCN plasmids.

Then, I will report how common structural variation is on these plasmids,
and report details of those structural variants.

"""

import subprocess
import json
import os
import sys
import gzip
import re
import csv
import math
from Bio import SeqIO
from os.path import basename, exists
import urllib.request
import time
import logging
import glob
import pprint


################################################################################
## Functions.

def get_SRA_ID_from_RefSeqID(refseq_id):
    ## datasets must be in $PATH.
    bash_command = f'datasets summary genome accession {refseq_id}'
    cmd_output = subprocess.check_output(bash_command, shell=True)
    json_output = cmd_output.decode('utf-8')
    json_data = json.loads(json_output)
    sra_id = "NA"
    reports = json_data.get("reports")
    if reports:
        sample_ids = reports[0].get("assembly_info", {}).get("biosample", {}).get("sample_ids")
        if sample_ids:
            for sample_id in sample_ids:
                if sample_id.get("db") == "SRA":
                    sra_id = sample_id.get("value")
                    break
    return(sra_id)


def get_RunID_and_longread_datatype_tuples(sra_id):
    ## pysradb must be in $PATH.
    pysradb_command = f'pysradb metadata {sra_id}'
    pysradb_attempts = 5
    pysra_command_worked = False
    while pysradb_attempts:
        try:
            pysradb_output = subprocess.check_output(pysradb_command, shell=True)
            ## assuming the previous line worked--
            pysradb_attempts = 0
            pysra_command_worked = True
        except subprocess.CalledProcessError:
            pysradb_attempts -= 1
    ## check to see if pysradb_output is meaningful.
    if pysra_command_worked:
        pysradb_output_str = pysradb_output.decode('utf-8')
        # Splits the metadata of the SRA ID into respective rows. 
        # And isolates the rows that use the Illumina instrument.
        rows = pysradb_output_str.strip().split('\n')
        run_accession_longreadtype_tuples = list()
        for row in rows:
            if ("OXFORD_NANOPORE" in row or "PACBIO_SMRT" in row) and ("WGS" in row):
                ## the Run_ID is the 3rd field from the end.
                my_run_accession = row.split("\t")[-3]
                if "OXFORD_NANOPORE" in row:
                    my_longreadtype = "OXFORD_NANOPORE"
                elif "PACBIO_SMRT" in row:
                    my_longreadtype = "PACBIO_SMRT"
                else:
                    raise AssertionError("ERROR: Unrecognized datatype")
                my_tuple = (my_run_accession, my_longreadtype)
                run_accession_longreadtype_tuples.append(my_tuple)
        if len(run_accession_longreadtype_tuples) == 0:
            print(f"No ONT or PacBio data found for accession {sra_id}")
            print("Non-Illumina PYSRADB DATA DUMP:")
            for row in rows:
                if "ILLUMINA" in row:
                    continue
                else:
                    print(row)
    return(run_accession_longreadtype_tuples)


def create_HighPCN_plasmid_SRA_RunID_table(high_PCN_plasmid_csv, RunID_table_outfile):
    ## first, get all RefSeq IDs in the ../data/high-PCN-plasmids.csv file.
    with open(high_PCN_plasmid_csv, "r") as plasmids_fh:
        plasmids_lines = plasmids_fh.read().splitlines()
    ## skip the header.
    plasmids_data = plasmids_lines[1:]
    ## get the first column containing RefSeqIDs_AssemblyAccession strings.
    AnnotationAccession_column = [line.split(",")[0] for line in plasmids_data]
    ## split to get the RefSeq IDs.
    refseq_ids = ["_".join(x.split("_")[:2]) for x in AnnotationAccession_column]
  
    ## now make the RunID csv file.
    with open(RunID_table_outfile, "w") as RunID_table_fh:
        header = "RefSeq_ID,SRA_ID,Run_ID,LongReadDataType\n"
        RunID_table_fh.write(header) 
        for RefSeq_accession in refseq_ids:
            my_SRA_ID = get_SRA_ID_from_RefSeqID(RefSeq_accession)
            if my_SRA_ID == "NA": continue ## skip genomes without SRA data.
            RunID_longread_datatype_tuples = get_RunID_and_longread_datatype_tuples(my_SRA_ID)
            for my_Run_ID, my_longread_datatype in RunID_longread_datatype_tuples:
                row = f"{RefSeq_accession},{my_SRA_ID},{my_Run_ID},{my_longread_datatype}\n"
                print(row) ## just to show that the program is running properly.
                RunID_table_fh.write(row)
    return


def create_refseq_accession_to_ftp_path_dict(prokaryotes_with_plasmids_file):
    refseq_accession_to_ftp_path_dict = dict()
    ## first, get all RefSeq IDs in the prokaryotes-with-plasmids.txt file.
    with open(prokaryotes_with_plasmids_file, "r") as prok_with_plasmids_file_obj:
        for i, line in enumerate(prok_with_plasmids_file_obj):
            if i == 0: continue ## skip the header.
            ## get the accession field (5th from end) and turn GCA Genbank IDs into GCF RefSeq IDs.
            refseq_id = line.split("\t")[-5].replace("GCA", "GCF")        
            ## get the ftp_url field (3rd from end) and make sure that we turn the GCA Genbank URL
            ## into the GCF RefSeq FTP URL.
            ftp_url = line.split("\t")[-3].replace("GCA", "GCF")
            ## check for for valid IDs and URLs (some rows have a '-' as a blank placeholder).
            if refseq_id.startswith("GCF") and refseq_id in ftp_url:
                refseq_accession_to_ftp_path_dict[refseq_id] = ftp_url
    return refseq_accession_to_ftp_path_dict


def reference_genome_passes_md5_checksum(gbff_gz_file, md5_file):
    with open(md5_file, "r") as checksum_fh:
        target_string = "_genomic.gbff.gz"
        for line in checksum_fh:
            if target_string in line:          
                my_target_checksum, my_target_filename = line.strip().split()
                break
    if sys.platform == "darwin":
        my_md5_cmd = "md5" ## run md5 on my mac,
    elif sys.platform == "linux":
        my_md5_cmd = "md5sum" ## but run md5sum on DCC (linux)
    else:
        raise AssertionError("UNKNOWN PLATFORM")
    ## run md5 on the local file and get the output.
    md5_call = subprocess.run([my_md5_cmd, gbff_gz_file], capture_output=True, text=True)
    my_md5_checksum = md5_call.stdout.split("=")[-1].strip()
    ## verify that the checksums match.
    return my_md5_checksum == my_target_checksum


def fetch_reference_genomes(RunID_table_file, refseq_accession_to_ftp_path_dict, reference_genome_dir):
    ## we get RefSeq IDs from the RunID table because this file *only* contains those RefSeq IDs 
    ## for which we could download raw Illumina short reads from the NCBI Short Read Archive.

    with open(RunID_table_file, "r") as RunID_file_obj:
        RunID_table_lines = RunID_file_obj.read().splitlines()

    ## remove the header from the imported data.
    RunID_table_data = RunID_table_lines[1:]
    ## get the first column to get all refseq_ids of interest.
    refseq_ids = [line.split(",")[0] for line in RunID_table_data]
    ## now look up the FTP URLs for each refseq id.
    ftp_paths = [refseq_accession_to_ftp_path_dict[x] for x in refseq_ids]

    for ftp_path in ftp_paths:
        ## note that the format of this accession is {refseqid}_{assemblyid}.
        my_full_accession = basename(ftp_path)
        my_base_filename = my_full_accession + "_genomic.gbff.gz"
        ## files on the NCBI FTP site to download
        gbff_ftp_path = os.path.join(ftp_path, my_base_filename)
        md5_ftp_path = os.path.join(ftp_path, "md5checksums.txt")
        ## local paths to download these files
        gbff_gz_file = os.path.join(reference_genome_dir, my_base_filename)
        md5_file = os.path.join(reference_genome_dir, my_full_accession + "_md5checksums.txt")
        
        if exists(gbff_gz_file) and exists(md5_file): ## then check whether the reference genome is OK.
            if reference_genome_passes_md5_checksum(gbff_gz_file, md5_file):
                continue
            else:
                os.remove(gbff_gz_file)
                os.remove(md5_file)
        
        gbff_fetch_attempts = 5
        gbff_fetched = False
        
        while not gbff_fetched and gbff_fetch_attempts:
            try:
                urllib.request.urlretrieve(gbff_ftp_path, filename=gbff_gz_file)
                urllib.request.urlretrieve(md5_ftp_path, filename=md5_file)                
            except urllib.error.URLError:
                ## if some problem happens, try again.
                gbff_fetch_attempts -= 1
                ## delete the corrupted files if they exist.
                if exists(gbff_gz_file):
                    os.remove(gbff_gz_file)
                if exists(md5_file):
                    os.remove(md5_file)
            ## if we are here, then assume the try block worked.
            if exists(gbff_gz_file) and exists(md5_file): ## then check whether the reference genome is OK.
                if reference_genome_passes_md5_checksum(gbff_gz_file, md5_file):
                    gbff_fetched = True  ## assume success if the checksum matches,
                    gbff_fetch_attempts = 0  ## and don't try again.
                else:
                    os.remove(gbff_gz_file)
                    os.remove(md5_file)
    return
 

def download_fastq_long_reads(SRA_data_dir, RunID_table_file):
        """
        the Run_ID has to be the last part of the directory.
        see documentation here:
        https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump
        """
        Run_IDs = []
        with open(RunID_table_file, "r") as RunID_table_file_obj:
            table_csv = csv.DictReader(RunID_table_file_obj)
            Run_IDs = [row["Run_ID"] for row in table_csv]
        for Run_ID in Run_IDs:
            prefetch_dir_path = os.path.join(SRA_data_dir, Run_ID)
            if os.path.exists(prefetch_dir_path): ## skip if we have already prefetched the read data.
                continue
            ## prefetch will create the prefetch_dir_path automatically-- give it the SRA_data_dir.
            prefetch_args = ["prefetch", Run_ID, "-O", SRA_data_dir]
            print (" ".join(prefetch_args))
            subprocess.run(prefetch_args)
        print("prefetch step completed.")
        my_cwd = os.getcwd()
        os.chdir(SRA_data_dir)
        for Run_ID in Run_IDs:
            sra_fastq_file = Run_ID + ".fastq"
            if os.path.exists(sra_fastq_file):
                continue
            else:
                print ("Generating fastq for: " + Run_ID)
                fasterq_dump_args = ["fasterq-dump", "--threads", "10", Run_ID]
                print(" ".join(fasterq_dump_args))
                subprocess.run(fasterq_dump_args)
        ## now change back to original working directory.
        os.chdir(my_cwd)
        return


def make_RefSeq_to_SRA_RunList_dict(RunID_table_csv):
    RefSeq_to_SRA_RunList_dict = dict()
    with open(RunID_table_csv, "r") as csv_fh:
        for i, line in enumerate(csv_fh):
            if i == 0: continue ## skip the header.
            line = line.strip() 
            RefSeqID, SRA_ID, RunID = line.split(',')
            if RefSeqID in RefSeq_to_SRA_RunList_dict:
                RefSeq_to_SRA_RunList_dict[RefSeqID].append(RunID)
            else:
                RefSeq_to_SRA_RunList_dict[RefSeqID] = [RunID]
    return RefSeq_to_SRA_RunList_dict


def generate_replicon_fasta_references_for_themisto(gbk_gz_path, fasta_outdir):
    print("reading in as input:", gbk_gz_path)
    ## open the input reference genome file.
    with gzip.open(gbk_gz_path, 'rt') as gbk_gz_fh:
        SeqID = None
        SeqType = None
        for i, record in enumerate(SeqIO.parse(gbk_gz_fh, "genbank")):
            SeqID = record.id
            if "chromosome" in record.description or i == 0:
                ## IMPORTANT: we assume here that the first record is a chromosome.
                SeqType = "chromosome"
            elif "plasmid" in record.description:
                SeqType = "plasmid"
            else:
                continue
            ## replace spaces with underscores in the replicon annotation field.
            replicon_description = record.description.replace(" ","_")
            header = ">" + "|".join(["SeqID="+SeqID,"SeqType="+SeqType,"replicon="+replicon_description])
            my_replicon_fastafile = SeqID + ".fna"
            my_replicon_outfilepath = os.path.join(fasta_outdir, my_replicon_fastafile)
            with open(my_replicon_outfilepath, "w") as outfh:
                outfh.write(header + "\n")
                outfh.write(str(record.seq) + "\n")


def generate_replicon_fasta_reference_list_file_for_themisto(fasta_outdir):
    genome_id = os.path.basename(fasta_outdir)
    replicon_fasta_filelist = [x for x in os.listdir(fasta_outdir) if x.endswith(".fna")]
    replicon_listfile = os.path.join(fasta_outdir, genome_id + ".txt")
    with open(replicon_listfile, "w") as fastatxtfile_fh:
        for fastafile in replicon_fasta_filelist:
            my_replicon_fasta_path = os.path.join(fasta_outdir, fastafile)
            fastatxtfile_fh.write(my_replicon_fasta_path + "\n")
    return


def make_NCBI_replicon_fasta_refs_for_themisto(refgenomes_dir, themisto_fasta_ref_outdir):
    ## this function makes a genome directory for each genome.
    ## each directory contains separate fasta files for each replicon.

    ## make the output directory if it does not exist.
    if not exists(themisto_fasta_ref_outdir):
        os.mkdir(themisto_fasta_ref_outdir)

    gzfilelist = [x for x in os.listdir(refgenomes_dir) if x.endswith("gbff.gz")]
    for gzfile in gzfilelist:
        gzpath = os.path.join(refgenomes_dir, gzfile)
        genome_id = gzfile.split(".gbff.gz")[0]
        fasta_outdir = os.path.join(themisto_fasta_ref_outdir, genome_id)
        ## make the fasta output directory if it does not exist.
        if not exists(fasta_outdir):
            os.mkdir(fasta_outdir)
        generate_replicon_fasta_references_for_themisto(gzpath, fasta_outdir)
        generate_replicon_fasta_reference_list_file_for_themisto(fasta_outdir)
    return


################################################################################
## Run the pipeline.
def main():

    run_log_file = "../results/plasmid-SV-pipeline-log.txt"
    ## Configure logging
    logging.basicConfig(filename=run_log_file, level=logging.INFO)

    high_PCN_plasmid_csv = "../data/high-PCN-plasmids.csv"
    RunID_table_csv = "../results/high-PCN-plasmid-RunID_table.csv"
    prokaryotes_with_plasmids_file = "../data/prokaryotes-with-chromosomes-and-plasmids.txt"
    reference_genome_dir = "../data/NCBI-reference-genomes/"
    SRA_data_dir = "../data/SRA/"


    #####################################################################################
    ## Stage 1: get SRA IDs and Run IDs for the RefSeq bacterial genomes containing high copy number plasmids.
    if exists(RunID_table_csv):
        Stage1DoneMessage = f"{RunID_table_csv} exists on disk-- skipping stage 1."
        print(Stage1DoneMessage)
        logging.info(Stage1DoneMessage)
    else:
        RunID_table_start_time = time.time()  # Record the start time
        create_HighPCN_plasmid_SRA_RunID_table(high_PCN_plasmid_csv, RunID_table_csv)        
        RunID_table_end_time = time.time()  # Record the end time
        RunID_table_execution_time = RunID_table_end_time - RunID_table_start_time
        Stage1TimeMessage = f"Stage 1 execution time: {RunID_table_execution_time} seconds"
        print(Stage1TimeMessage)
        logging.info(Stage1TimeMessage)
    
    #####################################################################################
    ## Stage 2: download reference genomes for each of the bacterial genomes containing plasmids,
    ## for which we can download Oxford Nanopore long-read data from the NCBI Sequencing Read Archive.
    ## first, make a dictionary from RefSeq accessions to ftp paths using the
    ## ../data/prokaryotes-with-chromosomes-and-plasmids.txt file.
    stage_2_complete_file = "../results/stage2.done"
    if exists(stage_2_complete_file):
        print(f"{stage_2_complete_file} exists on disk-- skipping stage 2.")
    else:
        refseq_accession_to_ftp_path_dict = create_refseq_accession_to_ftp_path_dict(prokaryotes_with_plasmids_file)
        ## now download the reference genomes for the high PCN plasmid genomes in the RunID_table_csv.
        fetch_reference_genomes(RunID_table_csv, refseq_accession_to_ftp_path_dict, reference_genome_dir)        
        with open(stage_2_complete_file, "w") as stage_2_complete_log:
            stage_2_complete_log.write("reference genomes downloaded successfully.\n")
            
    #####################################################################################
    ## Stage 3: download long reads for the genomes from the NCBI Short Read Archive (SRA).
    stage_3_complete_file = "../results/stage3.done"
    if exists(stage_3_complete_file):
        print(f"{stage_3_complete_file} exists on disk-- skipping stage 3.")
    else:
        SRA_download_start_time = time.time()  # Record the start time
        download_fastq_long_reads(SRA_data_dir, RunID_table_csv)
        SRA_download_end_time = time.time()  # Record the end time
        SRA_download_execution_time = SRA_download_end_time - SRA_download_start_time
        Stage3TimeMessage = f"Stage 3 (SRA download) execution time: {SRA_download_execution_time} seconds"
        print(Stage3TimeMessage)
        logging.info(Stage3TimeMessage)
        with open(stage_3_complete_file, "w") as stage_3_complete_log:
            stage_3_complete_log.write("SRA read data downloaded successfully.\n")

    #####################################################################################
    ## Stage 4: use minimap2 to align long reads to the high copy number plasmids. 
    stage_4_complete_file = "../results/stage4.done"
    if exists(stage_4_complete_file):
        print(f"{stage_4_complete_file} exists on disk-- skipping stage 4.")
    else:
        stage4_start_time = time.time()  # Record the start time

        
        ## CODE GOES HERE
        quit()
        
        stage4_end_time = time.time()  # Record the end time
        stage4_execution_time = stage4_end_time - stage4_start_time
        Stage4TimeMessage = f"Stage 4 execution time: {stage4_execution_time} seconds"
        print(Stage4TimeMessage)
        logging.info(Stage4TimeMessage)
        with open(stage_4_complete_file, "w") as stage_4_complete_log:
            stage_4_complete_log.write("Stage 4 complete.\n")

            
    return


main()


