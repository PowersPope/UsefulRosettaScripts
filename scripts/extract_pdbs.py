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
    args = p.parse_args()

    init(extra_options="-in:file:fullatom")

    # store our output
    df = pd.read_csv(args.input_file)

    for row in range(df.shape[0]):
        print("Output:", row)
        # extract row
        out = df.iloc[row,:]
        # Extract data
        file = out.file
        tag = out.tag

        # Setup the silentFile output
        opts = silent.SilentFileOptions()
        opts.in_fullatom(True)
        opts.set_binary_output(True)
        silentfile = silent.SilentFileData(opts)
        silentfile._read_file(os.path.join(args.silentfile_dir, file+".silent"))

        # grab the structure
        struct = silentfile.get_structure(tag)
        # set and fill pose
        pose = Pose()
        struct.fill_pose(pose)

        filename = tag + ".pdb"

        pose.dump_pdb(os.path.join(args.outfile_dir,filename))
        pose.clear()
