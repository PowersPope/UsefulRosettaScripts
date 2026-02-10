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
                       ) -> list[int]:
    """Cluster the structures in our SilentFile, to see how many unique clusters
    are present in our ensemble.
    
    PARAMS
    ------
    :args: The args from our argparse
    :insf: Our passed in SilentDataFile with our structures
    """
    decoy_names = insf.tags()
    nstructs = len(decoy_names)
    rmsd_check_list = np.zeros(nstructs, dtype=np.int32)
    cluster_map = dict()

    # list RMSD
#     rmsd_matrix = np.zeros((nstructs, nstructs), dtype=np.float32)
    skip_scaffolds = list()

    # take a pose (i) and (j) 
    pose_i = core.pose.Pose()
    pose_j = core.pose.Pose()
    for i in range(1, nstructs):
        # Init our vec and list for tracking this round
        if i in skip_scaffolds: continue
        rmsd_vec = np.zeros(nstructs - len(skip_scaffolds) - 2, dtype=np.float32)
        scaffolds_accumulated = list()  
        skip_scaffolds.append(i)
        scaffold_count = -1
        pose_i.clear()

        # Fill our pose
        silentstruct_i = insf.get_structure(decoy_names[i])
        silentstruct_i.fill_pose(pose_i)
        score = silentstruct_i.get_energy("score")
        print("Decoy:", decoy_names[i], "Score:", score)

        # Iter through our other availabel structures
        for j in range(2, nstructs):
            if j in skip_scaffolds: continue
            scaffold_count += 1

            # init empty poses
            silentstruct_j = insf.get_structure(decoy_names[j])
            silentstruct_j.fill_pose(pose_j)

            # Grab the Atom Map to compute RMSD
            bb_heavy_atommap = helpfunc.grab_atomid_map(
                pose_reference = pose_i,
                pose_target = pose_j,
                residue_anchors = list(range(1, pose_j.size()+1)),
                target_start_resi = 1,
                ref_chain = args.peptide_chain,
                bb_heavy_only = True,
            )
            # Superimpose
            rmsd_value = core.scoring.superimpose_pose(
                    pose_j,
                    pose_i,
                    bb_heavy_atommap,
                    )
            rmsd_vec[scaffold_count] = rmsd_value
#             rmsd_matrix[i-1, j-1] = rmsd_value
#             rmsd_matrix[j-1, i-1] = rmsd_value
            scaffolds_accumulated.append(j)
            pose_j.clear()
        rmsd_check = rmsd_vec <= args.rmsd_cutoff
        cluster_map[decoy_names[i]] = rmsd_check.sum()
        skip_scaffolds.extend(np.array(scaffolds_accumulated)[rmsd_check].tolist())

    print(f"RMSD Dict of Matches <= {args.rmsd_cutoff}:")
    print(cluster_map)
    print("Number of Unique Clusters:", len(cluster_map.keys()))

    return 0


if __name__ == "__main__":
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--insilent", type=str, default=None, help="path to input silentfile")
    p.add_argument("--rmsd-cutoff", type=float, default=0.3, help="How close in RMSD (A), do we want to count as something being in the same cluster.")
    p.add_argument("--peptide-chain", type=int, default=1, help="Which chain is our peptide? Needed for Oversat Filter")
    p.add_argument("--order-by-energy", "-obe", action="store_true", help="If you scored your structural ensemble, then use this flag to do energy based clustering. \
            This will segfault if you did not score the outputs!")
    args = p.parse_args()

    # Setup our initial Rosetta instance with our presets
    init(extra_options="-mute all -in:file:fullatom true")

    # A silent file is passed instead
    sfo = silent.SilentFileOptions()
    sfo.in_fullatom(True)
    sf_data = silent.SilentFileData(sfo)
    sf_data.read_file(args.insilent)
    if args.order_by_energy: 
        sf_data.order_by_energy()
        sf_data.renumber_all_decoys()

    # filter the outputs in our silentfile
    silentfile_process(args, sf_data)

