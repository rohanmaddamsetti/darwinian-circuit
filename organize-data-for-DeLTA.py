#!/usr/bin/env python

"""
organize-data-for-DeLTA.py by Rohan Maddamsetti.

This scripts takes a directory with the following structure:
--- KeyenceDataDir
    --- XY01
         --- T0001
              --- Image_T0001_CH4.tif
         --- T0002
              --- Image_T0002_CH4.tif
    --- XY02
         --- T0001
              --- Image_T0001_CH4.tif
         --- T0002
              --- Image_T0002_CH4.tif

And copies the images into a new directory with the following names:
--- DeLTADataDir
     --- XY01T0001.tif
     --- XY01T0002.tif
     --- XY02T0001.tif
     --- XY02T0002.tif

This can be parsed by DeLTA2 (https://gitlab.com/dunloplab/delta/) as follows:
# Init xpreader:
# (make sure you update prototype parameters if you use your own movie)
reader = delta.utilities.xpreader(
    extract_folder,
    prototype='DeLTADataDir/XY%02dT%04d.tif'
    fileorder='pt',
    filenamesindexing=1


Usage: python organize-data-for-DeLTA.py -o DeLTADataDir -i KeyenceDataDir

python organize-data-for-DeLTA.py -o ../results/line-hotel-results/20230410_RM7-140-1_RM7-140-21_LB_37C -i ../data/line-hotel-data/20230410_RM7-140-1_RM7-140-21_LB_37C/bestfocus
python organize-data-for-DeLTA.py -o ../results/line-hotel-results/20230419_RM7-140-1_RM7-140-21_LB_37C -i ../data/line-hotel-data/20230419_RM7-140-1_RM7-140-21_LB_37C/bestfocus


"""

import argparse
import os
import shutil

def main():
    parser = argparse.ArgumentParser(description="Copies a directory of images from Keyence Microscope into a directory formatted for DeLTA\nUsage: python organize-data-for-DeLTA.py -o DeLTADataDir -i KeyenceDataDir")
    parser.add_argument('-i', dest='inputdir', action='store', required=True, help="input directory.")
    parser.add_argument('-o', dest='outputdir', action='store', required=True, help="output directory.")
    args = parser.parse_args()
    
    ## print arguments to show user we got the right things.
    print("input directory =>", args.inputdir)
    print("output directory =>", args.outputdir)

    assert os.path.isdir(args.inputdir), "Input directory not found: " + args.inputdir
    
    if not os.path.exists(args.outputdir):
        ## then make the output directory.
        os.mkdir(args.outputdir)

    ## walk through the input directory.
    for my_path, my_directories, my_files in os.walk(args.inputdir):
        my_tif_files = [x for x in my_files if x.endswith(".tif") or x.endswith(".tiff")]
        if len(my_tif_files) == 0: continue ## skip if no tif files.
        pos_time_path = my_path.split(args.inputdir)[-1]
        empty_str, pos_str, time_str = pos_time_path.split('/')

        new_tif_filename = pos_str + time_str + ".tif"
        new_tif_path  = os.path.join(args.outputdir, new_tif_filename)
        ## important assumption: one tif file per position and time in this subdirectory.
        assert len(my_tif_files) == 1
        input_tif_file = my_tif_files[0]
        input_tif_path = os.path.join(my_path, input_tif_file)
        if not os.path.exists(new_tif_path):
            print("Copying: ",input_tif_path)
            print("to Destination: ", new_tif_path)
            shutil.copyfile(input_tif_path, new_tif_path)
        else:
            print(new_tif_path, " already exists-- no action taken.")

    return None

main()
