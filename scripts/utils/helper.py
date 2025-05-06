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
import pyrosetta.rosetta.core.pack.task.operation as opt
import pyrosetta.rosetta.core.io.silent as silent
from pyrosetta.rosetta.core.scoring import ScoreFunction
from pyrosetta.rosetta.std import map_core_id_AtomID_core_id_AtomID
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.protocols.denovo_design.movers import FastDesign
from pyrosetta.rosetta.core.select.movemap import MoveMapFactory
from pyrosetta.rosetta.core.pack.palette import CustomBaseTypePackerPalette

# python imports
import os
from typing import Union, List, Tuple


### Classes (Empty for now, as the funcitons should be modular enough to
### use within a class)
class SilentFileWrite:
    def __init__(
            self,
            outname: str = "test",
            ):
        # Setup the silentFile output
        self.opts = silent.SilentFileOptions()
        self.opts.in_fullatom(True)
        self.opts.set_binary_output(True)
        self.silentFile = silent.SilentFileData(self.opts)
        self.outname = outname

    def generate_plus_add_structure(self,
                      pose: core.pose.Pose,
                      structName: str,
                      ) -> int:
        """Take a structure and add it to our silentfiledata

        PARAMS
        ------
        :pose: The filled pose
        """
        struct = silent.BinarySilentStruct(
                self.opts, pose, structName,
                )
        self.silentFile.add_structure(struct)
        return 0

    def generate_plus_write_to_silentfile(
            self,
            pose: core.pose.Pose,
            structName: str,
            ) -> int:
        struct = silent.BinarySilentStruct(
                self.opts, pose, structName,
                )
        self.silentFile.write_silent_struct(struct, self.outname)
        return 0

    def write_all(self) -> int:
        """Write all of the added structues (use in conjunction with generate_plus_add_structure)
        """
        self.silentFile.write_all(self.outname)
        return 0

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

def modify_termini_to_cutpoints(pose: core.pose.Pose,
                                start: int = 1,
                                end: int = 8,
                                ) -> int:
    """
    Apply inplace the correct variant types to the N- & C-termini

    PARAMS
    ------
    :pose: Un cyclized pose
    :start: resi of the N-termini
    :end: resi of the C-termini
    """
    # Modify the variant types
    modifyvariant_nterm = sm.ModifyVariantTypeMover()
    modifyvariant_nterm.set_additional_type_to_add("CUTPOINT_UPPER")
    modifyvariant_nterm.set_residue_selector(rs.ResidueIndexSelector(start))
    modifyvariant_nterm.apply(pose)

    modifyvariant_cterm = sm.ModifyVariantTypeMover()
    modifyvariant_cterm.set_additional_type_to_add("CUTPOINT_LOWER")
    modifyvariant_cterm.set_residue_selector(rs.ResidueIndexSelector(end))
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
        focus_chain: int,
        distance: float = 4.0,
        include_focus_chain: bool = True,
        ) -> rs.NeighborhoodResidueSelector:
    """Grab the residue neighbors that are within some distance cutoff
    of the focus_chain

    PARAMS
    ------
    :focus_chain: The chain we want to use as our origin
    :distance: How far away from the origin, do we check to see
    :include_focus_chain: Whether to include the focus (focus_chain) or not
        in the resulting selection

    RETURNS
    -------
    :neighborhoodSel: The correct setup you would want to use for 
        downstream searches
    """
    # setup our neighborhood selector
    neighborSel = rs.NeighborhoodResidueSelector()
    neighborSel.set_focus_selector(
            rs.ChainSelector(focus_chain),
            )
    neighborSel.set_distance(distance)
    neighborSel.set_include_focus_in_subset(include_focus_chain)
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
        cartesian: bool = False,
        ) -> core.select.movemap.MoveMapFactory:
    mmf = core.select.movemap.MoveMapFactory()
    # Set all movements specified
    mmf.all_bb(bb)
    mmf.all_chi(chi)
    mmf.all_bondangles(bondangles and cartesian)
    mmf.all_bondlengths(bondlengths and cartesian)
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
    if specific_bondangles_selector and cartesian:
        mmf.add_bondangles_action(
                core.select.movemap.move_map_action.mm_enable,
                specific_bondangles_selector,
                )
    if specific_bondlengths_selector and cartesian:
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

def setup_taskfactory(
        task_operations: List[opt] = [],
        allow_dchiral: bool = False,
        ) -> TaskFactory:
    """Setup a taskfactory for use in design or relax

    PARAMS
    ------
    :task_operations: A list of task operations that will be passed to the tf.
        The order of the list matters, as the first item will be passed first.
        This is important because taskoperations can be overwritten.
        These should be specified in your specific script
    :allow_dchiral: Allow the design of dchiral amino acids

    RETURNS
    --------
    :tf: A setup taskfactory ready for use in FastDesign or FastRelax
    """
    # init our taskfactory
    tf = TaskFactory()

    # Good defaults to include
    # Carry over command line arguments like -ex1 -ex2
    tf.push_back(opt.InitializeFromCommandline())
    # Include the current rotamer from the pose into the first pack search
    tf.push_back(opt.IncludeCurrent())

    if allow_dchiral:
        dchiral_aminos = CustomBaseTypePackerPalette()
        # Setup dchiral_aminos
        dchiral_aminos.parse_additional_residue_types(
                "DALA,DLEU,DILE,DVAL,DTHR,DSER,DASP,DASN,DGLN,DGLU,DMET,DARG,DHIS,DLYS,DPHE,DTRP,DTYR,DPRO"
                )

        tf.set_packer_palette(dchiral_aminos)

    # Apply our list of operations
    for taskop in task_operations:
        tf.push_back(taskop)
    return tf

def hydrogen_bond_check(
        pose: core.pose.Pose,
        resi: int,
        ) -> bool:
    """Check if the current residue in our pose is forming any
    mainchain hydrogen bonds. Return True yes, False if not

    PARAMS
    ------
    :pose: Non empty pose
    :resi: The index of the residue to check (1-index)

    RETURNS
    -------
    True/False
    """
    # init a hbond set
    hbset = core.scoring.hbonds.HBondSet(
            pose = pose, 
            bb_only = True,
            )

    # Now check if it is either involved in donor/acceptor interactions
    if not (hbset.don_bbg_in_bb_bb_hbond(resi) or hbset.acc_bbg_in_bb_bb_hbond(resi)):
        return True
    else:
        return False

def mutate_residue(
        pose: core.pose.Pose,
        resi: int,
        residue_type: str,
        heterochiral: bool = True,
        ) -> int:
    """Check if the current residue in our pose is forming any
    mainchain hydrogen bonds. Return True yes, False if not

    PARAMS
    ------
    :pose: Non empty pose
    :resi: The index of the residue to check (1-index)
    :residue_type: the canonical residue type
    :heterochiral: if true then check for to see if the phi is pos/neg
        mutate to the correct enantiomer

    RETURNS
    -------
    Mutates in place the pose to our desired residue_type.
    """
    # first if check for heterochiral bool
    if heterochiral:
        # This is within the D chiral structural space, so lets
        # mutate to the correct enantiomer
        if pose.phi(resi) > 0:
            residue_type = "D"+residue_type
    
    # Now setup the mutate mover
    mutateRes = sm.MutateResidue()
    mutateRes.set_preserve_atom_coords(True)
    mutateRes.set_target(resi)
    mutateRes.set_res_name(residue_type)

    mutateRes.apply(pose)
    return 0

def search_proline_mutate(
        pose: core.pose.Pose,
        residue_selection: rs,
        number_prolines: int = 1,
        ) -> Tuple[bool, List[int]]:
    """Search the rama of residues that could be a proline
    if they are not invovled in hydrogen bonding then mutate until
    number_prolines is satisfied. For macrocycles having at least one
    proline can help with stability, given its restricted DoFs.

    PARAMS
    ------
    :pose: Our filled pose object
    :residue_selection: The residues we would like to scan
    :number_prolines: The number of desired prolines at the end

    RETURNS
    -------
    True/False: whether we were able to satifisy our number_prolines criteria
        or not
    """
    # init our count
    mutated_prolines = 0
    proline_positions = list()

    # grab the proline positions
    dpro_selector = rs.BinSelector()
    dpro_selector.set_bin_name("DPRO")
    dpro_selector.set_bin_params_file_name("PRO_DPRO")
    dpro_selector.initialize_and_check()
    pro_selector = rs.BinSelector()
    pro_selector.set_bin_name("LPRO")
    pro_selector.set_bin_params_file_name("PRO_DPRO")
    pro_selector.initialize_and_check()

    # only select the positions that are in our residue_selection
    specific_dprolines = rs.AndResidueSelector()
    specific_dprolines.add_residue_selector(residue_selection)
    specific_dprolines.add_residue_selector(dpro_selector)

    specific_lprolines = rs.AndResidueSelector()
    specific_lprolines.add_residue_selector(residue_selection)
    specific_lprolines.add_residue_selector(pro_selector)

    specific_prolines = rs.OrResidueSelector()
    specific_prolines.add_residue_selector(specific_lprolines)
    specific_prolines.add_residue_selector(specific_dprolines)


    # Get the resi numbers for looping
    specific_proline_resi = core.select.get_residues_from_subset(specific_prolines.apply(pose))

    # There are no positions that match our criteria within our selection. An empty list is returned, so exit
    if not specific_proline_resi:
        assert specific_proline_resi == {}, f"proline resi are not empty. There is a bug: {len(specific_proline_resi)}"
        return (False, proline_positions)

    for ir in specific_proline_resi:
        if not hydrogen_bond_check(pose, ir):
            mutate_residue(pose, ir, "PRO")
            mutated_prolines+=1
            proline_positions.append(ir)
            if mutated_prolines == number_prolines: return (True, proline_positions)
    # If we get to here, then we did not mutate enough positions to meet our criteria therefore the filter fails.
    # However, our pose was still possibly mutated if you have more then 1 desired proline. (keep in mind)
    return (False, proline_positions)

def generate_resfile_operation(
        filename: str,
        res_selector: rs,
        ) -> opt:
    """Generate a resfile operations task. Essentially this is useful for specifying
    L- & D-chiral amino acids (combine with PhiSelector)

    PARAMS
    ------
    :filename: The resfile
    :res_selector: The selection this resfile will apply to.

    RETURNS
    -------
    opt object that operates on that selection
    """
    resfile = opt.ReadResfile()
    resfile.filename(filename)
    resfile.set_residue_selector(res_selector)
    return resfile

def design_peptide(
        pose: core.pose.Pose,
        tf: TaskFactory,
        mmf: MoveMapFactory,
        scorefxn: ScoreFunction,
        cartesian: bool = False,
        script: str = "monomer",
        ) -> int:
    """Design our macrocycle pose, in the context of the target.
    TaskFactory and MoveMapFactory are alrady specified before this.

    PARAMS
    ------
    :pose: Our complex pose
    :tf: Our already specified TaskFactory object
    :mmf: Our mmf which determines what is allowed to move.
    :scorefxn: Our design score function (either cartesian or not)
    :cartesian: Allow cartesian minimization or not.
    :allow_dchiral: If true then allow dchiral design

    RETURNS
    -------
    Our inplace designed + relaxed pose
    """
    # init our design mover
    designMover = FastDesign(scorefxn_in = scorefxn, script_file="MonomerDesign2019" if script == "monomer" else "InterfaceDesign2019")
    # Set Cartesian non-ideal minimization
    designMover.cartesian(cartesian)
    # Set the minimize method
    designMover.min_type("lbfgs_armijo_nonmonotone")
    # Set the number of iterations
    designMover.max_iter(200) # Might need to be edited
    # Enable Design
    designMover.set_enable_design(True)
    
    # Set these based on cartesian passed setting
    designMover.minimize_bond_angles(cartesian)
    designMover.minimize_bond_lengths(cartesian)

    # set our passed movemap and taskfactories
    designMover.set_movemap_factory(mmf)
    designMover.set_task_factory(tf)

    # now apply to our pose
    designMover.apply(pose)
    return 0

def relax_selection(
        currpose: core.pose.Pose,
        res_selection: rs,
        scorefxn: ScoreFunction,
        script: str = "monomer",
        cartesian: bool = False,
        ) -> int:
    """Complex relaxation for the specific areas

    PARAMS
    ------
    :currpose: The current non-empty pose
    :res_selection: residue selection which is specified
    :scorefxn: The applied scorefunction
    :script: Monomer or Complex
    :cartesian: Allow cartesian 

    RETURNS
    -------
    Inplace relaxation of selection
    """
    # setup the mmf
    mmf = setup_mmf(
            specific_bb_selector = res_selection,
            specific_chi_selector = res_selection,
            specific_bondangles_selector = res_selection,
            specific_bondlengths_selector = res_selection,
            cartesian = cartesian,
            )

    # Now set up our FastRelax Mover
    fr = protocols.relax.FastRelax(scorefxn_in = scorefxn, script_file = "MonomerDesign2019" if script == "monomer" else "InterfaceDesign2019")
    fr.set_movemap_factory(mmf)
    fr.set_enable_design(False)

    # Now apply our relax to the pose
    fr.apply(currpose)
    return 0


def relax_sidechains(
        pose: core.pose.Pose,
        res_selection: rs,
        scorefxn: ScoreFunction,
        script: str = "monomer",
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
    fr = protocols.relax.FastRelax(scorefxn_in = scorefxn, script_file = "MonomerDesign2019" if script == "monomer" else "InterfaceDesign2019")
    fr.set_movemap_factory(mmf)
    fr.set_enable_design(False)
#     fr.set_scorefxn(scorefxn)

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
    emopts = core.scoring.methods.EnergyMethodOptions(scorefxn.energy_method_options())
    emopts.hbond_options().decompose_bb_hb_into_pair_energies(True)
    scorefxn.set_energy_method_options(emopts)
    return scorefxn

