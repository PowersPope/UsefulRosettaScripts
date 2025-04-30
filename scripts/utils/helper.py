#!/usr/bin/env python
#
# @author: Andrew Powers (apowers4@uoregon.edu) or (apowers@flatironinstitute.org)
# @description: Helper functions for use in PyRosetta scripts. Funcitons should be
#               modular.

# PyRosetta imports
import pyrosetta.io as io
import pyrosetta.rosetta as rosetta
import pyrosetta.rosetta.core as core
import pyrosetta.rosetta.core.select.residue_selector as rs
import pyrosetta.rosetta.protocols.simple_moves as sm
import pyrosetta.rosetta.protocols as protocols

# python imports
import collections
import os


### Classes


### Functions
def peptide_stub_mover(
        pose: core.pose.Pose, 
        residue_type: str,
        sequence_len: int,
        ) -> core.pose.Pose:
    """Add residues to a pose, if you don't use pose_from_sequence (You should use pose_from_sequence)

    PARAMS
    ------
    :pose: Our empty pose
    :residue_type: The residue type yo would like to add
    :sequence_len: The number of repeats that residue_type should be added

    RETURNS
    -------
    An inplace pose with residue_type added sequence_len times
    """
    # init our Mover using XML, since there are extra options that arent needed in XML
    # whereas PyRosetta wants them specified and they dont have the defaults.
    xml = f"""
    <ROSETTASCRIPTS>
        <MOVERS>
            j
            <PeptideStubMover name="extension" reset="false">
                <Append resname="{residue_type}" repeat="{sequence_len}"/>
            </PeptideStubMover>
        </MOVERS>
        <PROTOCOLS>
            <Add mover="extension"/>
        </PROTOCOLS>
    </ROSETTASCRIPTS>
    """
    pepStubMover = protocols.rosetta_scripts.XmlObjects.create_from_string(xml).get_mover("extension")
    pepStubMover.apply(pose)
    return pose.clone()

def set_macrocycle_root(sequence_len: int) -> int:
    """Determine the root residue of a macrocycle, based on AA length

    PARAMS
    ------
    :sequence_len: Number of Amino Acids

    RETURNS
    -------
    :root: The root residue idx (1-index)
    """
    if sequence_len % 2 == 0:
        return int(sequence_len/2)
    else:
        return int((sequence_len-1)/2)

def flip_trans_omegabonds(
        pose: core.pose.Pose,
        sequence_len: int,
        ) -> int:
    """Flip the omega angles in the backbone to trans orientation
    this needs to be done if PeptideStubMover was used first

    PARAMS
    ------
    :pose: Our filled in pose
    :sequence_len: The size of our peptide

    RETURNS
    -------
    A pose with trans omega bonds, as opposed to the default cis.
    """
    # First get our root
    root = set_macrocycle_root(sequence_len)
    # Unfortunately SetTorsion mover (which is used in RosettaScripts XML)
    # Is not fully converted well to PyRosetta, so we will specify it in XML
    xmlSetTorsion = f"""
    <ROSETTASCRIPTS>
        <MOVERS>
            <SetTorsion name="flip_trans" foldtree_root="{root}">
                <Torsion residue="{','.join([str(i) for i in range(1, sequence_len+1)])}" torsion_name="omega" angle="180"/>
            </SetTorsion>
        </MOVERS>
        <PROTOCOLS>
            <Add mover="flip_trans"/>
        </PROTOCOLS>
    </ROSETTASCRIPTS>
    """
    setTorsion = protocols.rosetta_scripts.XmlObjects.create_from_string(xmlSetTorsion).get_mover("flip_trans")
    setTorsion.apply(pose)

    return 0


def modify_termini_to_cutpoints(pose: core.pose.Pose) -> int:
    """
    Apply inplace the correct variant types to the N- & C-termini

    PARAMS
    ------
    :pose: Un cyclized pose
    """
    # Modify the variant types
    modifyvariant_nterm = sm.ModifyVariantTypeMover()
    modifyvariant_nterm.set_additional_type_to_add("CUTPOINT_UPPER")
    modifyvariant_nterm.set_residue_selector(rs.ResidueIndexSelector(1))
    modifyvariant_nterm.apply(pose)

    modifyvariant_cterm = sm.ModifyVariantTypeMover()
    modifyvariant_cterm.set_additional_type_to_add("CUTPOINT_LOWER")
    modifyvariant_cterm.set_residue_selector(rs.ResidueIndexSelector(pose.total_residue()))
    modifyvariant_cterm.apply(pose)

    pose.update_residue_neighbors()
    return 0

def declare_terminal_bond(pose: core.pose.Pose) -> int:
    """
    Fix terminal bond for macrocycles inplace

    PARAMS
    ------
    :pose: Pose object
    """
    # Fix termini bonds, so they are set correctly
    declarebond = sm.DeclareBond()
    declarebond.set(
        res1 = pose.total_residue(),
        atom1 = 'C',
        res2 = 1,
        atom2 ='N',
        add_termini = True,
    )
    declarebond.apply(pose)
    pose.update_residue_neighbors()
    return 0


def load_macrocycle_pdb(
        filename: str,
        chain_idx: int = 1,
        ) -> core.pose.Pose:
    """Load a macrocycle pdb in

    PARAMS
    ------
    :filename: path to file
    :chain_idx: The chain number the peptide is located in (1 is default)

    RETURNS
    -------
    :pose: Rosetta pose object
    """
    # load in our pdb
    pose = io.pose_from_pdb(filename)

    # Adjust termini cutpoints
    modify_termini_to_cutpoints(pose)
    
    # declare terminal bond
    declare_terminal_bond(pose)
    return pose
