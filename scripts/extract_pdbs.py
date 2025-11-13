#!/usr/bin/env python

from pyrosetta import init
import pyrosetta.rosetta.core.io.silent as silent
from pyrosetta.rosetta.core.pose import Pose


import argparse
import pandas as pd
import os

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--input-file", type=str, help="A csv file that has file name and description tag")
    p.add_argument("--silentfile-dir", type=str, help="Path to where silentfiles are housed")
    p.add_argument("--outfile-dir", type=str, help="Path to dump pdbs")
    p.add_argument("--single-input", action="store_true", help="Specify if only passing a csv with one silentfile")
    args = p.parse_args()

    init(extra_options="-in:file:fullatom true")

    # Make the outdir if it doesnt exist yet.
    if not os.path.exists(args.outfile_dir):
        os.makedirs(args.outfile_dir)

    # store our output
    df = pd.read_csv(args.input_file)

    # get list of unique silentfiles
    silentfiles = df.tag.unique()

    # iter through groups and make subdfs
    for sf in silentfiles:
        sub_df = df[df.tag == sf].reset_index()

        for row in range(sub_df.shape[0]):
            print("Output:", row)
            # extract row
            out = sub_df.iloc[row,:]
            # Extract data
            file = out.file
            tag = out.tag

    #         # Setup the silentFile output
            if row == 0 and args.single_input:
                opts = silent.SilentFileOptions()
                opts.in_fullatom(True)
                opts.set_binary_output(True)
                silentfile = silent.SilentFileData(opts)
                silentfile._read_file(os.path.join(args.silentfile_dir, file+".silent"))
#             elif not args.single_input:
#                 opts = silent.SilentFileOptions()
#                 opts.in_fullatom(True)
#                 opts.set_binary_output(True)
#                 silentfile = silent.SilentFileData(opts)
#                 silentfile._read_file(os.path.join(args.silentfile_dir, file+".silent"))

            # grab the structure
            struct = silentfile.get_structure(tag)
            # set and fill pose
            pose = Pose()
            struct.fill_pose(pose)

            filename = tag + ".pdb"

            pose.dump_pdb(os.path.join(args.outfile_dir,filename))
            pose.clear()
