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
from pyrosetta.rosetta.core.scoring import ScoreFunction
from pyrosetta.rosetta.core.import_pose.pose_stream import SilentFilePoseInputStream

def silentfile_process(sf: SilentFilePoseInputStream, args, 
                       scorefxn: ScoreFunction, outsf: silent.SilentFileData,
                       ) -> int:
    """Process the outputs for a silent file
    
    PARAMS
    ------
    :sf: Our input silentfile specified and loaded before this
    :args: The args from our argparse
    :scorefxn: Our init'd scorefunction for relax/scoring/oversat
    :outsf: The silentfile of filtered outputs

    RETURNS
    -------
    iters and processes the outputs based on our filter workflow
    """
    n=0
    # init an empty pose
    pose = core.pose.Pose()

    while sf.has_another_pose():
        # Fill our pose 
        sf.fill_pose(pose)

        # grab out metrics
        # Number of bb-bb hbonds
        try:
            bb_hbonds = pose.scores["bb_hbonds"]
        except:
            bb_hbonds = filterfuncs.determine_internal_bb_hbonds(
                    pose,
                    rs.ChainSelector(args.peptide_chain),
                    scorefxn,
                    "bb_hbonds",
                    )
        # Skip this input if it doesnt meet our criteria
        if bb_hbonds < args.bb_hbond_cutoff:
            continue

        # Check the score of our pose
        try:
            score = pose.scores["score"]
        except:
            score = scorefxn.score(pose)

        # Skip this input if it doesnt meet our score
        if score > args.score_cutoff:
            continue

        # Check if there is oversat hbond acceptors
        oversat = filterfuncs.oversat_filter(
                pose,
                scorefxn,
                rs.ChainSelector(args.peptide_chain),
                rs.ChainSelector(args.peptide_chain),
                True,
                )
        
        # Final if there is oversat write to our filtered silentfile
        if oversat:
            if args.pdb_output:
                pose.dump_pdb(
                        os.path.join(args.outpath, "filtered_"+core.pose.extract_tag_from_pose(pose)),
                        )
            else:
                outsf.generate_plus_add_structure(pose, core.pose.extract_tag_from_pose(pose))
            n+=1

    # Write all of the structs to a file
    print(n, "Structs passed the filtering process")
    outsf.write_all()

    return 0

def pdb_process(pdb_list: list[str], args,
                scorefxn: ScoreFunction, outsf: silent.SilentFileData,
                ) -> int:
    """Process the outputs for a pdb list
    
    PARAMS
    ------
    :pdb_list: Our input pdb list
    :args: The args from our argparse
    :scorefxn: Our init'd scorefunction for relax/scoring/oversat
    :outsf: The silentfile of filtered outputs

    RETURNS
    -------
    iters and processes the outputs based on our filter workflow
    """
    # num track
    n = 0

    # Loop through files and run our script through
    for f in pdb_list:
        pose = io.pose_from_pdb(f)

        # grab out metrics
        # Number of bb-bb hbonds
        try:
            bb_hbonds = pose.scores["bb_hbonds"]
        except:
            bb_hbonds = filterfuncs.determine_internal_bb_hbonds(
                    pose,
                    rs.ChainSelector(args.peptide_chain),
                    scorefxn,
                    "bb_hbonds",
                    )
        # Skip this input if it doesnt meet our criteria
        if bb_hbonds < args.bb_hbond_cutoff:
            continue

        # Check the score of our pose
        try:
            score = pose.scores["total_score"]
        except:
            score = scorefxn(pose)

        # Skip this input if it doesnt meet our score
        if score < args.score_cutoff:
            continue

        # Check if there is oversat hbond acceptors
        oversat = filterfuncs.oversat_filter(
                pose,
                scorefxn,
                rs.ChainSelector(args.peptide_chain),
                rs.ChainSelector(args.peptide_chain),
                True,
                )
        
        # Final if there is oversat write to our filtered silentfile
        if oversat:
            if args.pdb_output:
                pose.dump_pdb(
                        os.path.join(args.outpath, "filtered_"+pose.pdb_info().name().split('/')[-1]),
                        )
            else:
                outsf.generate_plus_add_structure(pose, pose.pdb_info().name())


            n+=1

    # Write all of the structs to a file
    print(n, "Structs passed the filtering process")

    if not args.pdb_output:
        outsf.write_all()

    return 0


if __name__ == "__main__":
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--inpath", type=str, default=None, help="path to PDB files")
    p.add_argument("--insilent", type=str, default=None, help="path to input silentfile")
    p.add_argument("--outpath", type=str, default="./", help="path to output our filtered silentfile")
    p.add_argument("--silentoutname", type=str, help="Output silent file name")
    p.add_argument("--relax", action="store_true", help="Perform relax on the full structure (chi & bb)")
    p.add_argument("--peptide-chain", type=int, default=1, help="Which chain is our peptide? Needed for Oversat Filter")
    p.add_argument("--bb-hbond-cutoff", type=int, default=2, help="Desired minimum bb-bb hbonbds")
    p.add_argument("--score-cutoff", type=float, default=0.0, help="Desired maximum Rosetta score")
    p.add_argument("--pdb-output", action="store_true", help="Write out PDBs instead of a silentfile")
    args = p.parse_args()

    # Setup our initial Rosetta instance with our presets
    init(extra_options="-mute all -in:file:fullatom true -out:file:silent_struct_type binary -ex1 -ex2aro -score:weights ref2015_cart")

    # Setup our silentfile
    if not args.pdb_output:
        outSilentFile = helpfunc.SilentFileWrite(outname=os.path.join(args.outpath,args.silentoutname))
    else:
        outSilentFile = None

    # setup our scorefxn
    scorefxn = helpfunc.init_scorefunction()

    if args.insilent == None:
        # Grab our pdb files
        files = glob.glob(os.path.join(args.inpath, "*.pdb"))

        # filter the outputs in our pdb list
        pdb_process(files, args, scorefxn, outSilentFile)
    else:
        # A silent file is passed instead
        files = SilentFilePoseInputStream(args.insilent)

        # filter the outputs in our silentfile
        silentfile_process(files, args, scorefxn, outSilentFile)

