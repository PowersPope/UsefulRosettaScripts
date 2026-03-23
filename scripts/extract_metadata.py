#!/usr/bin/env python
#
# @author: Andrew Powers (apowers4@uoregon.edu)
# @brief: A quick PyRosetta script for grabbing some metadata that is important for plotting and visualizing differences between structures within
#       a given silentfile.


# import
import argparse
import numpy as np
import os

# helper scripts
import utils.helper as helpfuncs

# Rosetta imports
from pyrosetta import init
import pyrosetta.io as io
import pyrosetta.rosetta as rosetta
import pyrosetta.rosetta.core as core
import pyrosetta.rosetta.protocols as protocols
import pyrosetta.rosetta.protocols.generalized_kinematic_closure as genkic
import pyrosetta.rosetta.protocols.simple_moves as sm
import pyrosetta.rosetta.core.select.residue_selector as select
from pyrosetta.rosetta.protocols.cyclic_peptide import PeptideInternalHbondsFilter
import pyrosetta.rosetta.core.io.silent as silent


def main() -> int:
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--insilent", type=str, help="Silentfile with library of cyclic peptide scaffolds.")
    p.add_argument("--output-path", type=str, help="Location to output files. Makes dirs if not already present.")
    p.add_argument("--outfile-prefix", type=str, default="out", help="The outfile name prefix, so that they are consistent across outputs.")
    p.add_argument("--pep-size", type=int, default=8, help="The size of the peptides within our --insilent file. Currently they should all be the same size, so that we can premake our numpy array to store our data in.")
    args = p.parse_args()

    # init rosetta
    init("-mute all -score:weights ref2015 -overwite -in:file:fullatom true")

    # A silent file is passed
    sfo = silent.SilentFileOptions()
    sfo.in_fullatom(True)
    sf_data = silent.SilentFileData(sfo)
    sf_data.read_file(args.insilent)
    tags = sf_data.tags()

    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)

    # init the files we will be writing to
    phi_psi_file = f"{args.outfile_prefix}-phi-psi.npy"

    # Array
    phi_psi = np.zeros( (len(tags), args.pep_size, 2), dtype=np.float)

    # Iter through our tags
    for i, tag in enumerate(tags):
        # load our pose in
        pose = core.pose.Pose()
        silentstruct = sf_data.get_structure(tag)
        silentstruct.fill_pose(pose)

        # Declare our peptide bond
        db = sm.DeclareBond()
        db.set(
            res1 = pose.total_residue(),
            atom1 = 'C',
            res2 = 1,
            atom2 = 'N',
            add_termini = False,
        )
        db.apply(pose)
        pose.update_residue_neighbors()

        # Extract out phi/psi information
        for ir in range(1, pose.total_residue()+1):
            phi_psi[i, ir-1, 0] = pose.phi(ir)
            phi_psi[i, ir-1, 1] = pose.psi(ir)

    np.save(os.path.join(args.output_path,phi_psi_file), phi_psi)

    return 0


if __name__ == "__main__":
    main()
