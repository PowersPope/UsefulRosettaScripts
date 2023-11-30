# Import necessary packages from pyrosetta
from pyrosetta import *

from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.core.pack.task import *
from pyrosetta.rosetta.protocols import *
from pyrosetta.rosetta.protocols.geometry import *
from pyrosetta.rosetta.protocols.relax import FastRelax

# Selector
from pyrosetta.rosetta.core.select import residue_selector as selections
from pyrosetta.rosetta.protocols.residue_selectors import StoredResidueSubsetSelector, StoreResidueSubsetMover

# Core
from pyrosetta.rosetta.core.pack.task import operation, TaskFactory

# Protocols
from pyrosetta.rosetta.protocols import minimization_packing as pack_min

# Import base python packages
import os, sys
import glob
import multiprocessing as mp
import argparse
from functools import partial

# Import data management
import pandas as pd
import numpy as np


# Define function to unbind the peptide and target
def unbind(pose, partners):
    """Perform translation movement to form unbound system

    PARAMS
    ------
        pose: pyrosetta.object
            bound complex pose
        partners: str
            Specifying the jump between A_B
    """
    # Set up the foldtree to allow unbinding
    docking.setup_foldtree(pose, partners, Vector1([-1, -1, -1]))
    # Set up the translation mover and select what to move
    trans_mover = rigid.RigidBodyTransMover(pose, 1)
    # set number of steps to move
    trans_mover.step_size(100)
    trans_mover.apply(pose)

def select_repack(pose):
    """Specify a packer that takes into account pre-translated interactions

    PARAMS
    ------
        pose: pyrosetta.object
            bound complex pose

    RETURNS
    -------
        packer: pyrosetta.minimization_mover
            Built packer for use
    """
    # Setup task factory
    tf = TaskFactory()
    tf.push_back(operation.InitializeFromCommandline())

    # Select chainB
    chainB = selections.ChainSelector("B")
    chainA = selections.ChainSelector("A")


    # generat inter selector and groups
    inter_selector = selections.InterGroupInterfaceByVectorSelector()
    inter_selector.group1_selector(chainB)
    inter_selector.group2_selector(chainA)

    # Set default params
    inter_selector.nearby_atom_cut(5.5)
    inter_selector.cb_dist_cut(11.0)
    inter_selector.vector_angle_cut(75.0)
    inter_selector.vector_dist_cut(9.0)

    # store selections
    store_selector = StoredResidueSubsetSelector("neighbor_pre_translation")

    # Generate selection storage & apply
    store_mover = StoreResidueSubsetMover(inter_selector,
                                          "neighbor_pre_translation",
                                          True
                                          )
    store_mover.apply(pose)

    # restrict particular residues to repack only
    tf.push_back(
        operation.OperateOnResidueSubset(
            operation.RestrictToRepackingRLT(),
            store_selector,
            False
        )
    )
    # flip selector to prevent repacking of everything unselected
    tf.push_back(
        operation.OperateOnResidueSubset(
            operation.PreventRepackingRLT(),
            store_selector,
            True
        )
    )

    packer = pack_min.PackRotamersMover()
    packer.task_factory(tf)

    return packer

def score_residues(pose, scorefxn):
    """Score each amino acid in the peptide for a general overview

    PARAMS
    ------
        pose: pyrosetta.object
            bound complex pose
        scorefxn: pyrosetta.score_function
            Set up scorefunction (needs to be changed in the code if you want to)

    RETURNS
    -------
        res_scores: dict[int: float]
            key - residue number of chain B only
            value - total reu score for that residue
    """
    # clone pose
    res_pose = pose.clone()

    # Grab start and end residues for loop
    start = pose.chain_begin(2)
    end = pose.chain_end(2)+1

    # create hash map to hold information
    res_scores = dict(zip(
        list(range(start,end)),
        [0]*(end-start)
    ))

    # score the pose
    scorefxn.score(res_pose)

    # grab the score for each residue
    for resi in range(start, end):
        # total score extraction for resi
        res_scores[resi] = res_pose.energies().residue_total_energy(resi)

    return res_scores

def setup_and_apply(file, path, relax_req, max_cpus, per_residue, cycpep, nstruct):
    """Perform ddg unbinding. Score complex, pull apart, repack, rescore

    PARAMS
    ------
        file: str
            Name of the pdb file
        path: str
            Path to dir that houses pdbs. (Passed in Argparse)
        relax_req: bool
            Determines if the whole complex is relaxed beforehand or not.
        max_cpus: int
            Maximum number of cpus allowed for multiprocessing
        per_residue: bool
            Determines if score_residues() is run or not
        cycpep: bool
            Determines if chains are specified, since bug is present for pyrosetta for cycpeps
        nstruct: int
            Number of repack and move rescore attempts

    RETURNS
    -------
        file: str
            Name of the pdb file
        bound: float
            REU score of initial bound complex
        unbound: float
            REU score of unbound system
        binding_ddg: float
            ddg between bound and unbound (bounnd - unbound)
    """
    # setup the scorefxn
    scorefxn = create_score_function("ref2015_cart")

    # setup partners
    partners = "A_B"

    # Make file pose and fill
    pose = Pose()
    pose_from_file(pose, file)

    # Apply fast relax
    fr = FastRelax()
    # set scorefxn with cart
    fr.set_scorefxn(scorefxn)
    # number of iterations
    fr.max_iter(20)

    # cycpeps need to be fixed
    if cycpep:
        # grab the foldtree as chain id is bugged, i can hack add it
        fold_tree = str(pose.fold_tree())
        pose.conformation().insert_chain_ending(
            int(fold_tree.split("-1")[0].split(' ')[-2])
        )

        # Make sure that the peptide has a bond
        cyclic_mover = pyrosetta.rosetta.protocols.cyclic_peptide.PeptideCyclizeMover()
        cyclic_mover.set_selector(selections.ChainSelector("B"))
        cyclic_mover.apply(pose)

    # Relax the full complex if it hasn't been before
    if relax_req:
        # Apply
        fr.apply(pose)
        # now dump the pose
        pose.dump_pdb(f'{file[:-4]}_bound.pdb')

    # extract per residue scores if requested
    if per_residue:
        # extract information
        res_scores = score_residues(pose, scorefxn)
    else:
        res_score = None

    # define repacking selectors and task object
    packer = select_repack(pose)
    # Set the score function
    packer.score_function(scorefxn)

    # track unbound and binding_ddg
    unbound = np.zeros(nstruct)
    binding_ddg = np.zeros(nstruct)

    for i in range(nstruct):
        # copy the pose
        test_pose = pose.clone()
        # Score the bound complex
        bound = scorefxn(test_pose)
        #unbind the complex
        unbind(test_pose, partners)

        # repack pose
        packer.apply(test_pose)

        # score the unbound
        unbound[i] = scorefxn(test_pose)
        # get ddG difference
        binding_ddg[i] = bound - unbound[i]


    # make unbound path if it doesn't exist
    unbound_path = os.path.join(path,"unbound_pdb")
    if not os.path.exists(unbound_path):
        os.mkdir(unbound_path)

    # If the unbdoung_pdb dir exists then check to see if the file already exists
    # if not then dump the pdb
    else:
        unbound_file_path = os.path.join(unbound_path, f"{file[:-4].split('/')[-1]}_unbound.pdb")
        if not os.path.exists(unbound_file_path):

            # Change directories for this
            os.chdir(unbound_path)

            # dump the unbound
            test_pose.dump_pdb(f"{file[:-4].split('/')[-1]}_unbound.pdb")

    return file, bound, unbound.mean(), unbound.var(), binding_ddg.mean(), binding_ddg.var(), res_scores


# Define main portion
def main():
    # Input args
    p = argparse.ArgumentParser()
    p.add_argument("--path_to_pdb_dir", help="Path to directory that houses, bound target + peptide complexes",
                   type=str, required=True)
    p.add_argument("--relax_pdb", action="store_true", help="Pass flag if you want to relax your pdb \
                   before you compute the ddg of binding")
    p.add_argument("--max_cpus", type=int, help="Set number of CPUs for Multiprocessing (default: os.cpu_count())",
                   default=os.cpu_count())
    p.add_argument("--score_per_residue", action="store_true", help="Determines if you get per residue \
                    scores or just binding information (default: False)")
    p.add_argument("--cycpep", action="store_true", help="Specify if working with cyclic peptides")
    p.add_argument("--nstruct", type=int, help="Number of attempts to compute the binding ddG (default: 10)",
                   default=10)
    args = p.parse_args()

    # init pyrosetta
    pyrosetta.init("-mute all")

    # Get GLOBAL variables from argparse
    PDB_DIR = args.path_to_pdb_dir
    RELAX = args.relax_pdb
    CPU_COUNT = args.max_cpus
    PER_RESIDUE = args.score_per_residue
    CYCPEP = args.cycpep
    NSTRUCT = args.nstruct

    # check path and list dir
    if os.path.exists(PDB_DIR):
        pdbs = glob.glob(os.path.join(PDB_DIR,"*.pdb"))
    else:
        raise Exception("----------------------- \
                        File Path Does Not EXIST \
                        ------------------------")


    # Count the number of available cpus to spin up instances on
    cpus = os.cpu_count()

    # Generate a pool and run the operation
    pool = mp.Pool(cpus)

    # Perform analysis
    file_out, bound_score, unbound_score, unbound_var, ddg_binding, ddg_binding_var, res_scores = zip(*pool.map(
        partial(
            setup_and_apply,
            path=PDB_DIR,
            relax_req=RELAX,
            max_cpus=CPU_COUNT,
            per_residue=PER_RESIDUE,
            cycpep=CYCPEP,
            nstruct=NSTRUCT,
        ),
        pdbs
    ))

    # Combine into dataframe
    ddg_bind_df = pd.DataFrame(
        {
            "file": file_out,
            "bound_score": bound_score,
            "unbound_score": unbound_score,
            "unbound_var": unbound_var,
            "ddg_binding_score": ddg_binding,
            "ddg_binding_var": ddg_binding_var,
        }
    )

    # If res_scores is set and filled then write out info to separate df and concat them
    if res_scores != None:
        # init new dict
        total_dict = dict()

        # loop through the keys
        for k in res_scores[0].keys():
            total_dict[k] = [temp_dict[k] for temp_dict in res_scores]

            #res scores
            res_df = pd.DataFrame(total_dict)

        # concat data
        ddg_bind_df = pd.concat([ddg_bind_df,res_df],axis=1)


    # write dataframe to csv
    ddg_bind_df.to_csv("ddg_top20_peptides.csv", index=False)

# Run function
if __name__ == "__main__":
    main()
