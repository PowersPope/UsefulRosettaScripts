#!/usr/bin/env python
#
# @author: Andrew Powers (apowers4@uoregon.edu)
# @brief: Check to see the diversity of a generated macrocycle ensemble

# import Packages
import os
from typing import TextIO
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
                       outfile: TextIO,
                       ) -> int:
    """Cluster the structures in our SilentFile, to see how many unique clusters
    are present in our ensemble.

    Logic:
        1. Grab i, j poses
        2. i is a cluster center
        3. now we loop through our list of structs (j)
            a. Check RMSD of i and j
                if <= rmsd_cutoff:
                    1. remove from list
                    2. Add to cluster center list of values
                else:
                    1. continue
    
    PARAMS
    ------
    :args: The args from our argparse
    :insf: Our passed in SilentDataFile with our structures
    """
    structure_list = list(insf.structure_list())
    decoy_names = list(insf.tags())

    # take a pose (i) and (j) 
    pose_i = core.pose.Pose()
    pose_j = core.pose.Pose()

    # Track variables
    cluster_count = 0

    # Loop to cluster everything
    while len(structure_list) != 0:
        pose_i.clear()

        # Fill our pose
        silentstruct_i = structure_list.pop(0)
        name = decoy_names.pop(0)
        silentstruct_i.fill_pose(pose_i)
        score = silentstruct_i.get_energy("score")
        print("Decoy:", name, "Score:", score)

        cluster_count += 1
        if args.output_cluster_centers:
            pose_i.dump_pdb(os.path.join(args.output_path, f"{name}-clustercenter-{cluster_count}.pdb"))

        rmsd_count = np.zeros(1, dtype=np.int32)
        for j in range(len(structure_list) - 1, -1, -1):
            silentstruct_j = structure_list[j]
            j_name = decoy_names[j]
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
            if rmsd_value <= args.rmsd_cutoff:
                rmsd_count += 1
                structure_list.pop(j)
                decoy_names.pop(j)
            pose_j.clear()

        print(f"{cluster_count},{rmsd_count},{score},{len(structure_list)}", file=outfile)

    outfile.close()
    return 0


if __name__ == "__main__":
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--insilent", type=str, default=None, help="path to input silentfile")
    p.add_argument("--rmsd-cutoff", type=float, default=0.3, help="How close in RMSD (A), do we want to count as something being in the same cluster.")
    p.add_argument("--peptide-chain", type=int, default=1, help="Which chain is our peptide? Needed for Oversat Filter")
    p.add_argument("--order-by-energy", "-obe", action="store_true", help="If you scored your structural ensemble, then use this flag to do energy based clustering. \
            This will segfault if you did not score the outputs!")
    p.add_argument("--output-cluster-centers", action="store_true", help="Output the cluster centers in pdbs, so that you can visualize the peptides.")
    p.add_argument("--output-path", default="./", type=str, help="Write .pdb and .txt files to this location. Will create the folders if it doesnt exist.")
    args = p.parse_args()

    # Create our out path if it doesnt exist
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)

    # Setup our initial Rosetta instance with our presets
    init(extra_options="-mute all -in:file:fullatom true")

    # A silent file is passed instead
    sfo = silent.SilentFileOptions()
    sfo.in_fullatom(True)
    sf_data = silent.SilentFileData(sfo)
    sf_data.read_file(args.insilent)
    if args.order_by_energy: 
        sf_data.order_by_energy()

    # filter the outputs in our silentfile
    cluster_out = open(os.path.join(args.output_path, "cluster_count.txt"), "w")
    print(f"Clustered with order by energy: {args.order_by_energy}, and a RMSD cutoff of: {args.rmsd_cutoff}", file=cluster_out)
    print(f"Cluster,Count,Score,NumOfStructuresLeft", file=cluster_out)
    silentfile_process(args, sf_data, cluster_out)

