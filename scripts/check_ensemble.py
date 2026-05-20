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
import pyrosetta.rosetta.protocols as protocols
import pyrosetta.rosetta.core as core
import pyrosetta.rosetta.core.io.silent as silent
import pyrosetta.rosetta.core.select.residue_selector as rs
import pyrosetta.io as io
from pyrosetta.rosetta.core.scoring import ScoreFunction

def read_in_and_score_data(
        args, 
        insf: silent.SilentFileData,
        ) -> tuple[list, list, int]:
    """Read in our data, apply thioether constraints, and score.
    Store the scores and structures separately for processing.

    PARAMS
    ------
    :args: The passed in arguments
    :insf: The input silentfile dataset

    RETURNS
    -------
    :out_scores: Scores of each pose in out_structures with their index mapped to the
        same struct in out_structures
    :out_structures: All poses with proper constraints applied
    :lowest_Eindex: The lowest energy structure index we can use to start our search
    """
    # Setup our return values
    out_scores, out_structures = list(), list()
    count = 0

    # Setup our scorefunction
    sfxn = core.scoring.get_score_function(is_fullatom=True)
    sfxn.set_weight(core.scoring.dihedral_constraint, 1.0)
    sfxn.set_weight(core.scoring.angle_constraint, 1.0)
    sfxn.set_weight(core.scoring.atom_pair_constraint, 1.0)

    print("Extracting score and structures from our silentfile, storing lowest energy index as our starting point...")
    for tag in insf.tags():
        count+= 1
        tmp_pose = core.pose.Pose()
        tag_silentstruct = insf.get_structure(tag)
        tag_silentstruct.fill_pose(tmp_pose)

        # Setup thioether constraints
        core.pose.remove_variant_type_from_pose_residue(
                tmp_pose,
                core.chemical.VariantType.UPPER_TERMINUS_VARIANT,
                tmp_pose.total_residue(),
                )

        core.pose.remove_lower_terminus_type_from_pose_residue(
                tmp_pose,
                tmp_pose.total_residue(),
                )
        thioether_sidechain_index = helpfunc.set_up_terminal_thioether_lariat_variants(tmp_pose)
        termini = protocols.simple_moves.DeclareBond()
        termini_updated = helpfunc.declare_thioether_bond_mover(tmp_pose, termini)
        termini_updated.apply(tmp_pose)
        helpfunc.declare_thioether_constraints(tmp_pose)
        tmp_pose.update_residue_neighbors()

        # Now score our pose
        sfxn(tmp_pose)
        tmp_score = tmp_pose.energies().total_energy()
        if (count==1) or (tmp_score < lowest_energy):
            lowest_energy = tmp_score
            lowest_Eindex = count
        out_scores.append( tmp_score )
        out_structures.append( tmp_pose.clone() )

    return out_scores, out_structures, lowest_Eindex

def cluster_structures(args,
                       score_list: list,
                       struct_list: list,
                       lowest_Eindex: int,
                       outfile: TextIO,
                       ) -> int:
    """Take in our extracted scores and poses, starting from the lowest energy
    struct build clusters

    PARAMS
    ------
    :score_list: Scores of each pose in struct_list with their index mapped to the
        same struct 
    :struct_list: All poses with proper constraints applied
    :lowest_Eindex: The lowest energy structure index we can use to start our search
    :outfile: THe file we will be writing our output to a TextIO object
    """
    # The number of structures we have left to cluster
    unclustered_structs = len(struct_list)
    clusters_dict = dict()
    # Amount of clusters we will use
    cluster_count = 0
    cluster_assignment = np.zeros(unclustered_structs)

    while unclustered_structs > 0:
        cluster_count += 1
        print(f"Starting on cluster: {cluster_count}...")

        # setup our cluster first
        cluster_assignment[lowest_Eindex] = cluster_count
        cluster_sortedbyenergy = [lowest_Eindex]
        unclustered_structs -= 1
        cluster_center = struct_list[lowest_Eindex]
        if args.output_cluster_centers:
            cluster_center.dump_pdb(os.path.join(args.output_path, f"clustercenter-{cluster_count}.pdb"))

        # generate a list of unassigned candidates
        unassigned_candidates = np.argwhere(cluster_assignment == 0).reshape(-1)

        # Determine distance from cluster center
        for s in unassigned_candidates:
            unassign_pose = struct_list[s]

            # Grab the Atom Map to compute RMSD
            bb_heavy_atommap = helpfunc.grab_atomid_map(
                pose_reference = cluster_center,
                pose_target = unassign_pose,
                residue_anchors = list(range(1, unassign_pose.size()+1)),
                target_start_resi = 1,
                ref_chain = args.peptide_chain,
                bb_heavy_only = True,
            )
            # Superimpose
            rmsd_value = core.scoring.superimpose_pose(
                    cluster_center,
                    unassign_pose,
                    bb_heavy_atommap,
                    )
            if rmsd_value <= args.rmsd_cutoff:
                print(f"Adding struct {s} to cluster center: {cluster_count}, RMSD distance of {rmsd_value:.2f}")
                cluster_assignment[s] = cluster_count
                cluster_sortedbyenergy.append(s)
                unclustered_structs -= 1
        # Store our information
        clusters_dict[cluster_count] = len(cluster_sortedbyenergy)
        print(f"Finished cluster: {cluster_count} with {len(cluster_sortedbyenergy)} matching structs...")

        # Search for next lowest score
        started_search = True
        lowest_energy = np.inf
        for j, score in enumerate(score_list):
            if cluster_assignment[j] != 0: continue
            elif started_search or score < lowest_energy:
                started_search = False
                lowest_energy = score
                lowest_Eindex = j
        print(f"{cluster_count},{args.rmsd_cutoff},{score_list[lowest_Eindex]},{len(cluster_sortedbyenergy)}", file=outfile)

    # Exit While loop
    print("--- Cluster Search Completed ---")
    outfile.close()
    return 0 


if __name__ == "__main__":
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--insilent", type=str, default=None, help="path to input silentfile")
    p.add_argument("--rmsd-cutoff", type=float, default=0.3, help="How close in RMSD (A), do we want to count as something being in the same cluster.")
    p.add_argument("--peptide-chain", type=int, default=1, help="Which chain is our peptide? Needed for Oversat Filter")
    p.add_argument("--output-cluster-centers", action="store_true", help="Output the cluster centers in pdbs, so that you can visualize the peptides.")
    p.add_argument("--output-path", default="./", type=str, help="Write .pdb and .txt files to this location. Will create the folders if it doesnt exist.")
    args = p.parse_args()

    # Create our out path if it doesnt exist
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)

    # Setup our initial Rosetta instance with our presets
    init(extra_options="-mute all -in:file:fullatom true -score:weights ref2015")

    # A silent file is passed instead
    sfo = silent.SilentFileOptions()
    sfo.in_fullatom(True)
    sf_data = silent.SilentFileData(sfo)
    sf_data.read_file(args.insilent)

    # Grab out and score everything
    out_scores, out_structs, lowest_Eindex = read_in_and_score_data(args, sf_data)

    # Setup the outputfile and the headers
    cluster_out = open(os.path.join(args.output_path, "cluster_count-new.txt"), "w")
    print(f"Clustered with order by energy: {args.order_by_energy}, and a RMSD cutoff of: {args.rmsd_cutoff}", file=cluster_out)
    print("Cluster,RMSDCutoff,Score,ClusterSize", file=cluster_out)

    # Cluster our data
    cluster_structures(args, out_scores, out_structs, lowest_Eindex, cluster_out)

