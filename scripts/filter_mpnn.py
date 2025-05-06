#!/usr/bin/env python
#
# @author: Andrew Powers (apowers4@uoregon.edu)
# @brief: Filter the MPNN outputs based on similar metrics I use for my own Rosetta based designs
#

# import Packages
import os, sys
import glob
import utils.helper as helpfunc
import utils.filters as filterfuncs
import argparse

from pyrosetta import init
import pyrosetta.rosetta.core as core
import pyrosetta.rosetta.core.select.residue_selector as rs
import pyrosetta.io as io


if __name__ == "__main__":
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--inpath", type=str, help="path to PDB files")
    p.add_argument("--silentoutname", type=str, help="Output silent file name")
    args = p.parse_args()

    init(extra_options="-mute all -out:file:silent_struct_type binary -ex1 -ex2aro -score:weights ref2015_cart")

    # Setup our silentfile
    outSilentFile = helpfunc.SilentFileWrite(outname=args.silentoutname)

    # setup our scorefxn
    scorefxn = helpfunc.init_scorefunction()

    # Grab our files
    files = glob.glob(os.path.join(args.inpath, "*.pdb"))

    n = 1
    # Loop through files and run our script through
    for f in files:
        pose = io.pose_from_pdb(f)

        # Clone our pose
        init_pose = pose.clone()

        # get the number of bb-bb hbonds
        bb_hbonds = filterfuncs.determine_internal_bb_hbonds(
                pose,
                rs.ChainSelector(1),
                scorefxn,
                "prerelax_bb_hbonds",
                )

        # Setup Movemap
        mmf = helpfunc.setup_mmf(bb=True, chi=True, bondangles=True, bondlengths=True,
                                 cartesian=True,)
        # Relax our output
        helpfunc.relax_selection(
                pose,
                rs.ChainSelector(1),
                scorefxn,
                "monomer",
                True,
                )

        rmsd = filterfuncs.compare_rmsd_pose(
                init_pose,
                pose,
                rs.ChainSelector(1),
                rs.ChainSelector(1),
                "after_relax_rmsd",
                )

        outSilentFile.generate_plus_add_structure(pose, f"struct_relax_{str(n).zfill(6)}")
        n+=1

    outSilentFile.write_all()
        








