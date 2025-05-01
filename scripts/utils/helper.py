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
from pyrosetta.rosetta.core.scoring import ScoreFunction
from pyrosetta.rosetta.std import map_core_id_AtomID_core_id_AtomID
from pyrosetta.rosetta.core.id import AtomID

# python imports
import collections
import os
from typing import Union, List


### Classes (Empty for now, as the funcitons should be modular enough to
### use within a class)

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

def remove_term_variants(pose: core.pose.Pose,
                        upper_idx: int,
                        lower_idx: int,
                        ) -> int:
    """
    Remove terminal variants

    PARAMS
    ------
    :pose: Un cyclized pose
    :upper_idx: The (1-index) residue index for the N-termini
    :lower_idx: The (1-index) residue index for the C-termini
    """
    # Modify the variant types
    modifyvariant_nterm = sm.ModifyVariantTypeMover()
    modifyvariant_nterm.set_additional_type_to_remove("CUTPOINT_UPPER")
    modifyvariant_nterm.set_additional_type_to_remove("LOWER_TERMINUS_VARIANT")
    modifyvariant_nterm.set_additional_type_to_remove("ACETYLATED_NTERMINUS_VARIANT")
    modifyvariant_nterm.set_additional_type_to_remove("ACETYLATED_NTERMINUS_CONNECTION_VARIANT")
    modifyvariant_nterm.set_residue_selector(rs.ResidueIndexSelector(upper_idx))
    modifyvariant_nterm.apply(pose)

    modifyvariant_cterm = sm.ModifyVariantTypeMover()
    modifyvariant_cterm.set_additional_type_to_remove("CUTPOINT_LOWER")
    modifyvariant_cterm.set_additional_type_to_remove("UPPER_TERMINUS_VARIANT")
    modifyvariant_cterm.set_additional_type_to_remove("DIMETHYLATED_CTERMINUS_VARIANT")
    modifyvariant_cterm.set_additional_type_to_remove("CTERM_AMIDATION")
    modifyvariant_cterm.set_residue_selector(rs.ResidueIndexSelector(lower_idx))
    modifyvariant_cterm.apply(pose)
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

def declare_terminal_bond(
        pose: core.pose.Pose,
        n_term: int,
        c_term: int,
        ) -> int:
    """
    Fix terminal bond for macrocycles inplace

    PARAMS
    ------
    :pose: Pose object
    :n_term: The n_termini idx
    :c_term: The c_termini idx
    """
    # Fix termini bonds, so they are set correctly
    declarebond = sm.DeclareBond()
    declarebond.set(
        res1 = c_term,
        atom1 = 'C',
        res2 = n_term,
        atom2 ='N',
        add_termini = True,
    )
    declarebond.apply(pose)
    pose.update_residue_neighbors()
    return 0

def define_root(size: int) -> int:
    """Define the root residue for a macorcycle

    PARAMS
    ------
    :size: Amount of residues in the pose

    RETURNS
    -------
    :root: The root anchor residue
    """
    if size % 2 == 0:
        return int(size / 2)
    else:
        return int((size-1)/2)

def set_foldtree(pose: core.pose.Pose) -> int:
    """Set the fold tree for our macrocycles"""
    # grab our root
    root = define_root(pose.total_residue())
    # grab our fold tree
    ft = core.kinematics.FoldTree()
    ft.clear()
    ft.add_edge(1, root, -1)
    ft.add_edge(root, pose.total_residue(), -1)
    ft.reorder(root)
    pose.fold_tree(ft)
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

def grab_neighbors(
        home_chain: int,
        distance: float = 4.0,
        ) -> rs.NeighborhoodResidueSelector:
    """Grab the residue neighbors that are within some distance cutoff
    of the home_chain

    PARAMS
    ------
    :home_chain: The chain we want to use as our origin
    :distance: How far away from the origin, do we check to see

    RETURNS
    -------
    :neighborhoodSel: The correct setup you would want to use for 
        downstream searches
    """
    # setup our neighborhood selector
    neighborSel = rs.NeighborhoodResidueSelector()
    neighborSel.set_focus_selector(
            rs.ChainSelector(home_chain),
            )
    neighborSel.set_distance(distance)
    return neighborSel

def setup_mmf(
        bb: bool = False,
        chi: bool = False,
        bondangles: bool = False,
        bondlengths: bool = False,
        nu: bool = False,
        jumps: bool = False,
        specific_bb_selector: Union[List, rs] = [],
        specific_chi_selector: Union[List, rs] = [],
        specific_bondangles_selector: Union[List, rs] = [],
        specific_bondlengths_selector: Union[List, rs] = [],
        specific_nu_selector: Union[List, rs] = [],
        ) -> core.select.movemap.MoveMapFactory:
    mmf = core.select.movemap.MoveMapFactory()
    # Set all movements specified
    mmf.all_bb(bb)
    mmf.all_chi(chi)
    mmf.all_bondangles(bondangles)
    mmf.all_bondlengths(bondlengths)
    mmf.all_branches(False)
    mmf.all_jumps(jumps)
    mmf.all_nu(nu)

    # Logic for specific movements
    if specific_bb_selector:
        mmf.add_bb_action(
                core.select.movemap.move_map_action.mm_enable,
                specific_bb_selector,
                )
    if specific_chi_selector:
        mmf.add_chi_action(
                core.select.movemap.move_map_action.mm_enable,
                specific_chi_selector,
                )
    if specific_bondangles_selector:
        mmf.add_bondangles_action(
                core.select.movemap.move_map_action.mm_enable,
                specific_bondangles_selector,
                )
    if specific_bondlengths_selector:
        mmf.add_bondlengths_action(
                core.select.movemap.move_map_action.mm_enable,
                specific_bondlengths_selector,
                )
    if specific_nu_selector:
        mmf.add_nu_action(
                core.select.movemap.move_map_action.mm_enable,
                specific_nu_selector,
                )
    return mmf

def grab_atomid_map(
        pose_reference: core.pose.Pose,
        pose_target: core.pose.Pose,
        residue_anchors: List[int],
        target_start_resi: int = 2,
        ) -> map_core_id_AtomID_core_id_AtomID:
    """Build an AtomID map for superimpose_pose call to function

    PARAMS
    ------
    :pose_reference: The receptor or target we designed macrocycles for. This should have
        the peptide or binder we used as a base/anchor for our peptide generation.
    :pose_target: Our peptide/macrocycle byitself not in the context of the target
    :residue_anchors: A list of residues that we want to align against (1-index)
    :target_start_resi: The starting resi of the targets motif (1-index)

    RETURNS
    -------
    :core_atomid_map: Our map of AtomIDs for superimpose_pose to function
    """
    # init our map id
    atid_map = map_core_id_AtomID_core_id_AtomID()

    # loop through our residues
    for ir in range(0, len(residue_anchors)):
        # Grab the residues
        fzn_posit = ir+target_start_resi
        nat_fzn_posit = ir+residue_anchors[0]
        curres = pose_target.residue(fzn_posit)
        refres = pose_reference.residue(nat_fzn_posit)

        # Assert
        assert curres.natoms() == refres.natoms(), "The number of atoms do not match"

        # Loop through atoms
        for ia in range(1, curres.natoms()+1):
            if curres.atom_is_hydrogen( ia ): continue
            atid_map[AtomID(ia, fzn_posit)] = AtomID(ia, nat_fzn_posit)
    return atid_map

def place_peptide_incontext(
        pose_reference: core.pose.Pose,
        pose_target: core.pose.Pose,
        ref_chain: int,
        residue_anchors: List[int],
        target_start_resi: int = 2,
        ) -> core.pose.Pose:
    """Place a peptide that is by itself in context with the target receptor.
        This will be done along the motif (residue_anchor)

    PARAMS
    ------
    :pose_reference: The receptor or target we designed macrocycles for. This should have
        the peptide or binder we used as a base/anchor for our peptide generation.
    :pose_target: Our peptide/macrocycle byitself not in the context of the target
    :ref_chain: I am going to write this as if there are only 2 chains present
    :residue_anchors: A list of residues that we want to align against (1-index)
    :target_start_resi: The starting resi of the targets motif (1-index)

    RETURNS
    -------
    :target_peptide: A new pose where our receptor and macrocycle/peptide are now
        placed together.
    """
    # generate an atom map 
    atom_map = grab_atomid_map(
            pose_reference,
            pose_target,
            residue_anchors,
            target_start_resi,
            )

    # Superimpose
    core.scoring.superimpose_pose(
            pose_target,
            pose_reference,
            atom_map,
            )

    # Now delete our ref chain
    ref_pose = pose_reference.clone()
    ref_pose.delete_residue_range_slow(
            ref_pose.chain_begin(ref_chain),
            ref_pose.chain_end(ref_chain),
            )

    # Now combine the two poses
    core.pose.append_pose_to_pose(
            pose_target,
            ref_pose,
            new_chain = True,
            )

    pose_target.update_residue_neighbors()
    # This added because disulfides werent detected and found here:
    # https://forum.rosettacommons.org/node/3932
    pose_target.conformation().detect_disulfides()
    return pose_target.clone()


def load_silentfile(filepath: str) -> core.import_pose.pose_stream.PoseInputStream:
    assert os.path.exists(filepath), "The passed silentfile does not exist make sure \
            it is spelled correctly."
    silentfile = core.import_pose.pose_stream.SilentFilePoseInputStream(filepath)
    return silentfile

def relax_sidechains(
        pose: core.pose.Pose,
        res_selection: rs,
        scorefxn: ScoreFunction,
        ) -> int:
    """Take a pose and relax the residues that are selected within
    the passed selector. Backbone will be fixed

    PARAMS
    ------
    :pose: Our non-empty pose
    :res_selection: The residue selector that has already been specified
    :scorefxn: A scorefuntion that is used for scoring the relax movements

    RETURNS
    -------
    performs in place relaxation of the specified residues in a pose
    """
    # First setup a movemapfactory
    mmf = setup_mmf(specific_chi_selector=res_selection)

    # Now set up our FastRelax Mover
    fr = protocols.relax.FastRelax()
    fr.set_movemap_factory(mmf)
    fr.set_enable_design(False)
    fr.set_scorefxn(scorefxn)

    # Now apply our relax to the pose
    fr.apply(pose)
    return 0

def init_scorefunction(default: bool=False,) -> ScoreFunction:
    """Define and init our scorefunction;
    This will be based on a macrocycle, and making sure the right terms
    are defined to keep geometry correct. This should be tailored to your
    use case.

    PARAMS
    ------
    :default: This makes it so that the default scorefunction is made,
        instead of the macrocycle specific one

    RETURNS
    -------
    :scorefxn: Our defined scorefunction object
    """
    # init class global variables
    scorefxn = core.scoring.get_score_function(is_fullatom = True)
    # Stop here if we want default
    if default:
        return scorefxn
    scorefxn.set_weight(core.scoring.hbond_lr_bb, 1.0)
    scorefxn.set_weight(core.scoring.hbond_sr_bb, 1.0)
    scorefxn.set_weight(core.scoring.coordinate_constraint, 1.0)
    scorefxn.set_weight(core.scoring.atom_pair_constraint, 1.0)
    scorefxn.set_weight(core.scoring.dihedral_constraint, 1.0)
    scorefxn.set_weight(core.scoring.angle_constraint, 1.0)
    scorefxn.set_weight(core.scoring.chainbreak, 1.0)
    emopts = core.scoring.methods.EnergyMethodOptions( scorefxn.energy_method_options() ) # make a copy
    emopts.hbond_options().decompose_bb_hb_into_pair_energies( True )
    scorefxn.set_energy_method_options( emopts )
    return scorefxn

