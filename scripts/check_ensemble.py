#!/usr/bin/env python
#
# @author: Andrew Powers (apowers4@uoregon.edu)
# @brief: Check to see the diversity of a generated macrocycle ensemble

# import Packages
import os
import glob
import utils.helper as helpfunc
import utils.filters as filterfuncs
import argparse
import numpy as np

from pyrosetta import init
import pyrosetta.rosetta.core as core
import pyrosetta.rosetta.core.io.silent as silent
import pyrosetta.rosetta.core.select.residue_selector as rs
import pyrosetta.io as io
from pyrosetta.rosetta.core.scoring import ScoreFunction


def silentfile_process(args, 
                       insf: silent.SilentFileData,
#                        outsf: silent.SilentFileData,
                       decoy_names: list[str],
                       ) -> list[int]:
    """Process the outputs for a silent file
    
    PARAMS
    ------
    :args: The args from our argparse
    :scorefxn: Our init'd scorefunction for relax/scoring/oversat
    :outsf: The silentfile of filtered outputs
    :start_idx: The index to represent the scaffold we are using to compare against
    :rmsd_mat: RMSD full matrix. This will be returned at the end
    :decoy_names: rmsd list of decoy names in our SilentDataFile

    RETURNS
    -------
    iters and processes the outputs based on our filter workflow
    """
    nstructs = len(decoy_names)

    # list RMSD
    rmsd_matrix = np.zeros((nstructs, nstructs), dtype=torch.float32)

    # take a pose (i) and (j) 
    pose_i = core.pose.Pose()
    pose_j = core.pose.Pose()
    for i in range(nstructs):
        pose_i.clear()
        silentstruct_i = insf.get_structure(decoy_names[i])
        silentstruct_i.fill_pose(pose_i)

        for j in range(1, nstructs):
            # init empty poses
            silentstruct_j = insf.get_structure(decoy_names[j])
            silentstruct_j.fill_pose(pose_j)

            # Grab the Atom Map to compute RMSD
            bb_heavy_atommap = helpfuncs.grab_atomid_map(
                pose_reference = pose_i,
                pose_target = pose_j,
                residue_anchors = list(range(1, pose_j.size()+1)),
                target_start_resi = 1,
                ref_chain = 1,
                bb_heavy_only = True,
            )
            # Superimpose
            rmsd_value = core.scoring.superimpose_pose(
                    pose_j,
                    pose_i,
                    bb_heavy_atommap,
                    )
            rmsd_matrix[i, j] = rmsd_value
            rmsd_matrix[j, i] = rmsd_value
            pose_j.clear()
        rmsd_check = rmsd_matrix[i, :] <= 0.3
        print("RMSD CHECK <= 0.3:", rmsd_check)
        print(rmsd_check.sum())

    return 0


if __name__ == "__main__":
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--insilent", type=str, default=None, help="path to input silentfile")
    p.add_argument("--outpath", type=str, default="./", help="path to output our filtered silentfile")
    p.add_argument("--silentoutname", type=str, help="Output silent file name")
    p.add_argument("--peptide-chain", type=int, default=1, help="Which chain is our peptide? Needed for Oversat Filter")
    args = p.parse_args()

    # Setup our initial Rosetta instance with our presets
    init(extra_options="-mute all -in:file:fullatom true -out:file:silent_struct_type binary -ex1 -ex2aro -score:weights ref2015_cart")

    # Setup our silentfile
    outSilentFile = helpfunc.SilentFileWrite(outname=os.path.join(args.outpath,args.silentoutname))

    # setup our scorefxn
    scorefxn = helpfunc.init_scorefunction()

    # A silent file is passed instead
    sf_data = silent.SilentFileData(silent.SilentFileOptions())
    sf_data.read_file(args.insilent)
    decoy_names = sf_data.structure_list()

    # filter the outputs in our silentfile
    silentfile_process(args, sf_data, decoy_names)

