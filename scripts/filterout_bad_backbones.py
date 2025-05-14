#!/usr/bin/env python
#
# @author: Andrew Powers (apowers4@uoregon.edu)
# @brief: Filter out bad backbone outputs from whatever your generation process is (cyclicCAE or GenKIC)

# import Packages
import os
import glob
import utils.helper as helpfunc
import utils.filters as filterfuncs
import argparse

from pyrosetta import init
import pyrosetta.rosetta.core as core
import pyrosetta.rosetta.core.io.silent as silent
import pyrosetta.rosetta.core.select.residue_selector as rs
import pyrosetta.io as io

def silentfile_process(sf: silent.SilentFileData, args) -> int:
    """Process the outputs for a silent file
    
    PARAMS
    ------
    :sf: Our input silentfile specified and loaded before this
    :args: The args from our argparse

    RETURNS
    -------
    iters and processes the outputs based on our filter workflow
    """
    return 0

def pdb_process(pdb_list: list[str], args) -> int:
    """Process the outputs for a pdb list
    
    PARAMS
    ------
    :pdb_list: Our input pdb list
    :args: The args from our argparse

    RETURNS
    -------
    iters and processes the outputs based on our filter workflow
    """
    n = 1
    # Loop through files and run our script through
    for f in pdb_list:
        pose = io.pose_from_pdb(f)

        # Clone our pose
        init_pose = pose.clone()

        if args.relax:
            # get the number of bb-bb hbonds (before relaxing if relax is set)
            bb_hbonds = filterfuncs.determine_internal_bb_hbonds(
                    pose,
                    rs.ChainSelector(1),
                    scorefxn,
                    "prerelax_bb_hbonds",
                    )
    return 0


if __name__ == "__main__":
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--inpath", type=str, default=None, help="path to PDB files")
    p.add_argument("--insilent", type=str, default=None, help="path to input silentfile")
    p.add_argument("--silentoutname", type=str, help="Output silent file name")
    p.add_argument("--relax", action="store_true", help="Perform relax on the full structure (chi & bb)")
    args = p.parse_args()

    # Setup our initial Rosetta instance with our presets
    init(extra_options="-mute all -out:file:silent_struct_type binary -ex1 -ex2aro -score:weights ref2015_cart")

    # Setup our silentfile
    outSilentFile = helpfunc.SilentFileWrite(outname=args.silentoutname)

    # setup our scorefxn
    scorefxn = helpfunc.init_scorefunction()

    if args.insilent == None:
        # Grab our pdb files
        files = glob.glob(os.path.join(args.inpath, "*.pdb"))

        # filter the outputs in our pdb list
        pdb_process(files, args)
    else:
        # A silent file is passed instead
        opts = silent.SilentFileOptions()
        opts.in_fullatom(True)
        opts.set_binary_output(True)
        files = silent.SilentFileData(opts)
        files._read_file(args.insilent)

        # filter the outputs in our silentfile
        silentfile_process(files, args)

