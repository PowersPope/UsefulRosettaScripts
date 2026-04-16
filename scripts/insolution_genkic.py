#!/usr/bin/env python
#
#

# Packages
import argparse
from typing import List, Tuple
import datetime

# decorators
from utils.helper import flip_trans_omegabonds, peptide_stub_mover
from utils.decorators import timeit

# Rosetta imports
from pyrosetta import init
import pyrosetta.rosetta as rosetta
import pyrosetta.io as io
import pyrosetta.rosetta.numeric as numeric
import pyrosetta.rosetta.core as core
import pyrosetta.rosetta.protocols as protocols
import pyrosetta.rosetta.protocols.generalized_kinematic_closure as genkic
import pyrosetta.rosetta.protocols.simple_moves as sm
import pyrosetta.rosetta.core.select.residue_selector as select
from pyrosetta.rosetta.protocols.cyclic_peptide import PeptideInternalHbondsFilter
import pyrosetta.rosetta.core.io.silent as silent


# Class Generation
class BackboneGeneration:
    def __init__(
        self,
        debug: bool = False,
        randomize_root: bool = True,
        time_test: bool = False,
        empty_pose_test: bool = False,
        thioether: bool = False,
        homochiral: bool = False,
        thioether_chloroacetyl_homochiral: bool = False,
        thioether_chloroacetyl_dchiral: bool = False,
        lariat_sidechain_index: int = 0,
        ):
        """
        Setup and init our Pyrosetta instance to generate macrocycle backbones + Design

        PARAMS
        ------
        :debug: Output Debug Print Messages
        :randomize_root: This should be set if your root residue (which is an anchor) needs
            to be randomized. If you had a target, you would not randomize this
        :time_test: Comparison for XML vs PyRosetta. Since XML writes to disk every struct
        :empty_pose_test: Instead of using pose_from_sequence, we will create an empty
            pose and then append to it.
        :thioether: True -> thioether peptides are generated, False -> N-C peptide bond peptides
        :homochiral: True -> poly-alanine peptides, False -> poly-glycine peptides
        :thioether_chloroacetyl_homochiral: Whether the N-term residue is alanine (True) or glycine (false)
        :lariat_sidechain_index: The position of our lariat sidechain atom
        """
        # Set up our instance of PyRosetta
        if debug:
#             init("-in:file:fullatom true -out:level 2000 -score:weights ref2015 -out:file:silent_struct_type binary -overwrite")
            init("-in:file:fullatom true -score:weights ref2015 -out:file:silent_struct_type binary -overwrite -precompute_ig")
        else:
            init("-mute all -in:file:fullatom true -score:weights ref2015 -out:file:silent_struct_type binary -overwrite -precompute_ig")

        # init class global variables
        self.scorefxn = core.scoring.get_score_function(is_fullatom = True)
        self.scorefxn.set_weight(core.scoring.hbond_lr_bb, 1.0)
        self.scorefxn.set_weight(core.scoring.hbond_sr_bb, 1.0)
#         self.scorefxn.set_weight(core.scoring.coordinate_constraint, 1.0)
        self.scorefxn.set_weight(core.scoring.atom_pair_constraint, 1.0)
        self.scorefxn.set_weight(core.scoring.dihedral_constraint, 1.0)
        self.scorefxn.set_weight(core.scoring.angle_constraint, 1.0)
        self.scorefxn.set_weight(core.scoring.chainbreak, 1.0)
        emopts = core.scoring.methods.EnergyMethodOptions( self.scorefxn.energy_method_options() ) # make a copy
        emopts.hbond_options().decompose_bb_hb_into_pair_energies( True )
        self.scorefxn.set_energy_method_options( emopts )

        # ParsedProtocol Scorefxns
        self.sfxn_highhbond = self.scorefxn.clone()
        self.sfxn_highhbond.set_weight(core.scoring.hbond_lr_bb, 10.0)
        self.sfxn_highhbond.set_weight(core.scoring.hbond_sr_bb, 10.0)
        self.sfxn_highhbond.set_energy_method_options( emopts )
        # Cart Sfxn
        self.sfxn_highhbond_cart = self.sfxn_highhbond.clone()
        self.sfxn_highhbond_cart.set_energy_method_options( emopts )
        if self.sfxn_highhbond_cart.get_weight(core.scoring.cart_bonded) == 0.0:
            self.sfxn_highhbond_cart.set_weight(core.scoring.cart_bonded, 0.5)
            self.sfxn_highhbond_cart.set_weight(core.scoring.pro_close, 0.0)


        # randomize our root residue option
        self.randomize_root = randomize_root
        self.time_test = time_test
        self.DEBUG = debug
        self.empty_pose_test = empty_pose_test
        self.lariat_sidechain_index = lariat_sidechain_index
        self.thioether = thioether
        self.homochiral = homochiral
        self.thioether_chloroacetyl_homochiral = thioether_chloroacetyl_homochiral
        self.thioether_chloroacetyl_dchiral = thioether_chloroacetyl_dchiral



    def grab_pivot_anchor_and_cyclization_points(self, 
                            pose: core.pose.Pose,
                            ) -> tuple[list[int], int, int, int]:
        """
        Determine the pivot and root residues, this is done randomly, to allow different loops

        PARAMS
        ------
        :pose: Our peptide pose

        RETURNS
        -------
        :pivote_res: A list of pivot residues [first, middle, last]
        :anchor_res: The root of our peptide
        :cyclization_point_end: The cyclization points of our peptide.
            This is based on if it is n-c or thioether (this is generally
            the first residue)
        :cyclization_point_start: The cycliztion starting position
            This is based on if it is n-c or thioether (this is generally
            the last residue)
        """
        # Extract out the length of our peptide
        peptide_length = pose.total_residue()

        # Setup our method
        cyclization_point_end, cyclization_point_start = 1, peptide_length

        if self.thioether:
            cyclization_point_end, cyclization_point_start = self.find_first_and_last_thioether_lariat_residues(pose)

        # Anchor residue range and selection
        anchor_res_min = cyclization_point_end + 2 if self.thioether else 3
        anchor_res_max = cyclization_point_start - 2 if self.thioether else peptide_length - 2
        anchor_res = numeric.random.rg().random_range(anchor_res_min, anchor_res_max)
        first_loop_res, last_loop_res = anchor_res + 1, anchor_res - 1

        assert cyclization_point_end < last_loop_res
        assert first_loop_res < cyclization_point_start

        # Randomly pick a residue to be the middle pivot residue. Cant be first, anchor, or last res
        middle_loop_res = numeric.random.rg().random_range(cyclization_point_end, cyclization_point_start-3)
        if middle_loop_res == last_loop_res: middle_loop_res += 3
        elif middle_loop_res == anchor_res: middle_loop_res += 2
        elif middle_loop_res == first_loop_res: middle_loop_res += 1
        if middle_loop_res > peptide_length: middle_loop_res -= peptide_length

        # set our pivot residues
        pivot_res = [first_loop_res, middle_loop_res, last_loop_res]

        return pivot_res, anchor_res, cyclization_point_end, cyclization_point_start

    def get_nonroot_residues(self, peptide_length: int, root: int) -> List[int]:
        """
        Grab the residue indexes in the order that they need to be in for GenKIC as in
        from root -> peptide_length -> 1 -> root (where root is excluded)

        PARAMS
        ------
        :peptide_length: Amount of amino acids in our peptide
        :root: The root residue idx

        RETURNS
        -------
        :non_root: A list of residues that are not the root residue. Also,
            in the order that needs to be specified for genkic
        """
        # init list
        non_root = [i for i in range(root+1, peptide_length+1)] # post root
        pre_root = [i for i in range(1, root)]
        non_root.extend(pre_root)
        return non_root

    def remove_term_variants(
            self,
            pose: core.pose.Pose,
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
        modifyvariant_nterm.set_residue_selector(select.ResidueIndexSelector(1))
        modifyvariant_nterm.apply(pose)

        modifyvariant_cterm = sm.ModifyVariantTypeMover()
        modifyvariant_cterm.set_additional_type_to_remove("CUTPOINT_LOWER")
        modifyvariant_cterm.set_additional_type_to_remove("UPPER_TERMINUS_VARIANT")
        modifyvariant_cterm.set_additional_type_to_remove("DIMETHYLATED_CTERMINUS_VARIANT")
        modifyvariant_cterm.set_additional_type_to_remove("CTERM_AMIDATION")
        modifyvariant_cterm.set_additional_type_to_remove("CTERM_CONNECT")
        modifyvariant_cterm.set_residue_selector(select.ResidueIndexSelector(pose.total_residue()))
        modifyvariant_cterm.apply(pose)
        return 0

    def apply_hbond_filter(self, pose: io.Pose) -> Tuple[bool, int]:
        """
        Count internal bb-bb filter where we get if it meets our cutoff

        PARAMS
        ------
        :pose: Glycine minimized pose

        RETURNS
        -------
        :filter_check: The bool output of if our peptide passes
        :hbond_count: The number of bb-bb hbonds interaction
        """
        # Get the length of our pose
        pep_len = pose.total_residue()
        root = self.foldtree_define(pep_len)

        ### init variables
        cutoff_hbonds = pep_len // 3  # we want this many internal hbonds (6: 1, 7:2, 8:2, 9:3)

        ### Setup Filters for our backbones
        # Number of internal Hbonds
        hbond_filter = PeptideInternalHbondsFilter()
        hbond_filter.set_hbond_cutoff(cutoff_hbonds)
        hbond_filter.set_hbond_types(
            backbone_backbone_setting = True,
            backbone_sidechain_setting = False,
            sidechain_sidechain_setting = False,
        )
        hbond_filter.set_scorefxn(self.scorefxn)

        # bool out
        filter_check = hbond_filter.apply(pose)
        hbond_count = hbond_filter.report_sm(pose)
        pose.scores["bb-hbonds"] = hbond_count
        return filter_check

    def randomize_mainchain(self, pose: core.pose.Pose, n_to_c_cyclize: bool) -> int:
        """Randomize the main chain torsions of our pose"""
        self.debug: print("---Randomizing mainchain torsions.---")
        nres = pose.total_residue()
        rama = core.scoring.ScoringManager.get_instance().get_Ramachandran()
        ramaprepro = core.scoring.ScoringManager.get_instance().get_RamaPrePro()

        for i in range(1, nres+1):
            if pose.residue_type(i).is_alpha_aa() or pose.residue_type(i).is_peptoid():
                # Extract out the correct number of torsions vector1_double
                rand_torsions = rosetta.utility.vector1_double(len(pose.residue(i).mainchain_torsions()) - 1)
                # Use a fake alanine for following_rsd for the cterm residue unless you're doing
                # cterm cyclization and thus you can use its upper-connected residue.
                if i == nres and not n_to_c_cyclize:
                    following_rsd = core.chemical.ResidueTypeFinder(
                            core.chemical.ChemicalManager.get_instance().residue_type_set("fa_standard")
                            ).residue_base_name("ALA").get_representative_type()
                else:
                    following_rsd = pose.residue_type(pose.residue(i).residue_connection_partner(
                        pose.residue(i).upper_connect().index()
                        ))
                # torsion values based on the ramaprepro (preferenced by pose conf, resi, and resi+1)
                ramaprepro.random_mainchain_torsions( 
                                                     pose.conformation(), 
                                                     pose.residue_type(i), 
                                                     following_rsd, 
                                                     rand_torsions,
                                                     )
                for j in range(1, len(rand_torsions)+1): 
                    pose.set_torsion( core.id.TorsionID( i, core.id.TorsionType.BB, j ), rand_torsions[j] )
                if i!=nres: pose.set_omega(i, 180.0)
                if pose.residue_type(i).is_oligourea(): pose.set_mu(i, 180.0)
        return 0
    
#     def randomize_mainchain(self, pose: core.pose.Pose) -> int:
#         """Randomize the main chain torsions of our pose"""
#         self.debug: print("---Randomizing mainchain torsions.---")
#         nres = pose.total_residue()
#         rama = core.scoring.ScoringManager.get_instance().get_Ramachandran()
#         ramaprepro = core.scoring.ScoringManager.get_instance().get_RamaPrePro()
# 
#         for i in range(1, nres+1):
#             if pose.residue_type(i).is_alpha_aa() or pose.residue_type(i).is_peptoid():
#                 rand_torsions = pose.residue(i).mainchain_torsions()
#                 # Use a fake alanine for following_rsd for the cterm residue unless you're doing
#                 # cterm cyclization and thus you can use its upper-connected residue.
#                 if i == nres:
#                     following_rsd = core.chemical.ResidueTypeFinder(
#                             core.chemical.ChemicalManager.get_instance().residue_type_set("fa_standard")
#                             ).residue_base_name("ALA").get_representative_type()
#                 else:
#                     following_rsd = pose.residue_type_ptr(
#                         pose.residue(i).residue_connection_partner( pose.residue(i).upper_connect().index() )
#                     )
#                 ramaprepro.random_mainchain_torsions( pose.conformation(), pose.residue_type_ptr(i), following_rsd, rand_torsions )
#                 covered_torsions = ramaprepro.get_mainchain_torsions_covered( pose.conformation(), pose.residue_type_ptr(i), following_rsd )
#                 assert len(pose.residue(i).mainchain_torsions()) - 1 == len(rand_torsions) - 1, "This is wrong"
#                 for j in range(1, len(rand_torsions)):
#                     if covered_torsions[j] == 0.0:
#                         # Copy torsion values for torsions that are not covered by the mainchain potential -- i.e. don't set these to 0.
#                         rand_torsions[j] = pose.residue(i).mainchain_torsions()[j] 
#                     else: # Using classic rama tables for sampling:
#                         if pose.residue(i).backbone_aa() != core.chemical.aa_unk:
#                             rama.random_phipsi_from_rama(pose.residue(i).backbone_aa(), rand_torsions[1], rand_torsions[2])
#                         else:
#                             rama.random_phipsi_from_rama( pose.residue_type(i).aa(), rand_torsions[1], rand_torsions[2])
#                 assert len(rand_torsions) - 1 == len(pose.residue(i).mainchain_torsions()) - 1, "Wrong here"
#                 for j in range(1, len(rand_torsions)): 
#                     pose.set_torsion( core.id.TorsionID( i, core.id.BB, j ), rand_torsions[j] )
#                 if i!=nres: pose.set_omega(i, 180.0)
#                 if pose.residue_type(i).is_oligourea(): pose.set_mu(i, 180.0);
#             else: #If this is not a recognized type:
#                 for j in range(1, len(pose.residue(i).mainchain_torsions())): # Loop through all mainchain torsions.
#                     if i==nres and j==len(pose.residue(i).mainchain_torsions())-1: continue; # Skip the last mainchain torsion (not a DOF).
#                     setting = 180.0
#                     if j!=len(pose.residue(i).mainchain_torsions())-1:
#                         setting = numeric.random.rg().uniform()*360.0 - 180.0
#                     pose.set_torsion( core.id.TorsionID(i, core.id.BB, j), setting );
#         return 0

    def anchor_randomizebyrama(self,
                               pose: io.Pose,
                               anchor_resi: int,
                               ) -> int:
        """
        Generate randomize Anchor position by Rama Prepro Mover, since GenKIC does not alter
        it by default

        PARAMS
        ------
        :pose: Our Pose with the residue to randomize
        :anchor_resi: The residue index position of our anchor (1-indexing)
        """
        randomizeBB = protocols.backbone_moves.RandomizeBBByRamaPrePro()
        randomizeBB.set_residue_selector(
            select.ResidueIndexSelector(anchor_resi)
        )
        randomizeBB.apply(pose)
        return 0

    def setup_pp(self, pose: core.pose.Pose) -> protocols.rosetta_scripts.ParsedProtocol:
        # PP
        pp = protocols.rosetta_scripts.ParsedProtocol()

        # Setup and Update OH Mover
        if self.thioether:
            termini_close = protocols.simple_moves.DeclareBond()
            termini_close_update = self.declare_thioether_bond_mover(pose, termini_close)
            pp.add_step(termini_close_update, "Update_cyclization_point_polymer_dependent1", None)

        else:
            termini_close = protocols.simple_moves.DeclareBond()
            termini_close_update = self.declare_terminal_bond_mover(pose, termini_close)
            pp.add_step(termini_close_update, "Update_cyclization_point_polymer_dependent1", None)

        # Check oversat
        oversat1 = protocols.cyclic_peptide.OversaturatedHbondAcceptorFilter()
        oversat1.set_scorefxn( self.scorefxn )
        oversat1.set_hbond_energy_cutoff( -0.1 ) 
        pp.add_step(None, "Oversaturated_Hbond_Acceptors", oversat1)

#         # Set up a movemap
#         mm = core.kinematics.MoveMap()
#         mm.set_bb(True)
#         mm.set_chi(True)
# 
#         # Set up minimover
#         minmover = protocols.minimization_packing.MinMover()
#         minmover.movemap(mm)
#         minmover.score_function(self.sfxn_highhbond)
#         minmover.min_type("lbfgs_armijo_nonmonotone")
#         minmover.tolerance(1.0e-7)
# #         frlx = protocols.relax.FastRelax( self.sfxn_highhbond, 3 )
#         pp.add_step(minmover, "High_Hbond_MinMover", None)
# 
#         # Cartesian relax
# #         frlx_cart = protocols.relax.FastRelax( self.sfxn_highhbond_cart, 3 )
# #         frlx_cart.minimize_bond_angles(True)
#         # Set up minimover
#         minmover_cart = protocols.minimization_packing.MinMover()
#         minmover.movemap(mm)
#         minmover.score_function(self.sfxn_highhbond_cart)
#         minmover.min_type("lbfgs_armijo_nonmonotone")
#         minmover.tolerance(1.0e-7)
#         pp.add_step(minmover_cart, "High_Hbond_FastRelax_Cartesian", None)

        # Update OH again
        if self.thioether:
            termini_close2 = protocols.simple_moves.DeclareBond()
            termini_close2_updated = self.declare_thioether_bond_mover(pose, termini_close2)
            pp.add_step(termini_close2_updated, "Update_cyclization_point_polymer_dependent2", None)
        else:
            termini_close2 = protocols.simple_moves.DeclareBond()
            termini_close2_update = self.declare_terminal_bond_mover(pose, termini_close2)
            pp.add_step(termini_close2_update, "Update_cyclization_point_polymer_dependent2", None)

        # Oversat check final
        oversat2 = protocols.cyclic_peptide.OversaturatedHbondAcceptorFilter()
        oversat2.set_scorefxn( self.scorefxn )
        oversat2.set_hbond_energy_cutoff( -0.1 ) 
        pp.add_step(None, "Oversaturated_Hbond_Acceptors", oversat2)

        return pp

    def apply_genkic(self, pose: io.Pose, lariat_sample_cis_: bool = False) -> io.Pose:
        """
        Apply the generalized kinematic loop closure to a pose. This will designate the appropriate
        size genkic to apply

        PARAMS
        ------
        :pose: Our input glycine pose that has been set up
        :lariat_sample_cis_: Determines if our CA-N-CO-CP2 bond is cis (true) or trans (false)

        RETURNS
        -------
        :genkic_pose: A stochastically sampled backbone given sequence RAMA preferences
        :genkic_succ: Whether genkic found a successful solution
        """
#         # Get the length of our pose
#         if thioether:
#             pep_len = self.find_last_disulf_res(pose)
#         else:
#             pep_len = pose.total_residue()
#         root = self.foldtree_define(pep_len)

        # Calculate which residues to perturb and set as pivots
        pivot_res, root, cyclization_point_end, cyclization_point_start = self.grab_pivot_anchor_and_cyclization_points(
                pose,
                )

        ft = core.kinematics.FoldTree()
        ft.clear()
        ft.add_edge(root, 1, -1)
        ft.add_edge(root, pose.total_residue(), -1)
        ft.reorder(root)
        pose.fold_tree(ft)

        # free_residues, pivot_res = self.residues_to_perturb(pep_len, root, thioether)

        # Calculate residues to include in GenKIC
#         non_root_residues = self.get_nonroot_residues(pep_len, root)

        # setup our ParsedProtocol
        pp = self.setup_pp(pose)

        # init the genkic class object
        GenKIC = genkic.GeneralizedKIC()
        GenKIC.set_closure_attempts(500)
        GenKIC.set_min_solution_count(1)
        GenKIC.set_selector_type("lowest_energy_selector")
        GenKIC.set_selector_scorefunction(self.sfxn_highhbond)
        GenKIC.set_preselection_mover(pp)
        GenKIC.set_correct_polymer_dependent_atoms(True)
        # Add bb randomization for Anchor (rama prepro) if doing selection
        if self.randomize_root:
            if self.DEBUG: print("RANDOMIZE ROOT RESIDUE (THIS IS ONLY DONE FOR IN SOLUTION GENERATION)\nApplying To Pose Brefore GenKIC Closure...")
            self.anchor_randomizebyrama(pose, root)

        #Define the loop residues
        for i in range(pivot_res[0], cyclization_point_start+1): GenKIC.add_loop_residue(i)
        for i in range(cyclization_point_end, pivot_res[-1]+1): GenKIC.add_loop_residue(i)

        # Define tail residues
        if self.thioether:
            if cyclization_point_end > 1:
                for i in range(cyclization_point_end-1, 0, -1): GenKIC.add_tail_residue(i)
            if cyclization_point_start < pose.total_residue():
                for i in range(cyclization_point_start+1, pose.total_residue()+1): GenKIC.add_tail_residue(i)

        self.debug: print("--Before Set Pivots--")
        GenKIC.set_pivot_atoms(
            rsd1 = pivot_res[0],
            at1 = "CA",
            rsd2 = pivot_res[1],
            at2 = "CA",
            rsd3 = pivot_res[2],
            at3 = "CA",
        )

        # Add close bond logic
        if self.thioether:
            # Set the omega bond to trans
            # This only produces CA-N-CO-CP2 bonds without randomization
            # trans should be energetically favored, so setting bond to trans here if option specified
            if not lariat_sample_cis_:
                self.debug: print("---Set Omega bonds to trans---")
                GenKIC.add_perturber( protocols.generalized_kinematic_closure.perturber.set_dihedral )
                omega_vect = rosetta.utility.vector1_core_id_NamedAtomID()
                omega_vect.append( core.id.NamedAtomID( "N", cyclization_point_end) )
                omega_vect.append( core.id.NamedAtomID( "CO",cyclization_point_end) )
                GenKIC.add_atomset_to_perturber_atomset_list( omega_vect )
                GenKIC.add_value_to_perturber_value_list( 180.0 )

            # Setup atomsets of thieother to perturb
            self.debug: print("---Setup Atomset of thioether to perturb---")
            GenKIC.add_perturber( protocols.generalized_kinematic_closure.perturber.randomize_dihedral )
            curtype = pose.residue_type( cyclization_point_start )
            first_sc_at = curtype.first_sidechain_atom()
            curat = curtype.residue_connect_atom_index( curtype.n_possible_residue_connections() ) 
            while curat != first_sc_at:
                otherat = curtype.icoor(curat).stub_atom1().atomno()
                sc_vect = rosetta.utility.vector1_core_id_NamedAtomID()
                sc_vect.append( core.id.NamedAtomID( curtype.atom_name(curat), cyclization_point_start ) )
                sc_vect.append( core.id.NamedAtomID( curtype.atom_name(otherat), cyclization_point_start ) )
                GenKIC.add_atomset_to_perturber_atomset_list( sc_vect )
                curat = curtype.icoor(curat).stub_atom1().atomno()

            # Randomize backbone of upper res
            self.debug: print("---Ranomize Upper---")
            GenKIC.add_perturber( protocols.generalized_kinematic_closure.perturber.randomize_dihedral )
            bb_vect = rosetta.utility.vector1_core_id_NamedAtomID()
            if curtype.is_alpha_aa() or curtype.is_oligourea() or curtype.is_peptoid() or curtype.is_beta_aa():
                bb_vect.append ( core.id.NamedAtomID( "N", cyclization_point_start ) )
                bb_vect.append ( core.id.NamedAtomID( "CA", cyclization_point_start ) )
            elif curtype.is_gamma_aa():
                bb_vect.append ( core.id.NamedAtomID( "N", cyclization_point_start ) )
                bb_vect.append ( core.id.NamedAtomID( "C4", cyclization_point_start ) )
            GenKIC.add_atomset_to_perturber_atomset_list( bb_vect )

            # Randomize Lower res backbone (Only going to make for alpha AA and peptoids at the moment)
            self.debug: print("---Ranomize Lower---")
            GenKIC.add_perturber( protocols.generalized_kinematic_closure.perturber.randomize_dihedral )
            if pose.residue_type( cyclization_point_end ).is_alpha_aa() or pose.residue_type( cyclization_point_end ).is_peptoid():
                bb_vect = rosetta.utility.vector1_core_id_NamedAtomID()
                bb_vect.append( core.id.NamedAtomID( "CA", cyclization_point_end ) )
                bb_vect.append( core.id.NamedAtomID( "C", cyclization_point_end ) )
                GenKIC.add_atomset_to_perturber_atomset_list( bb_vect )

            # Setup our close bond logic
            GenKIC.close_bond(
                rsd1 = cyclization_point_start,
                at1 = pose.residue_type(cyclization_point_start).get_disulfide_atom_name(),
                rsd2 = cyclization_point_end,
                at2 = "CP2",
                bondlength = 1.827,  # taken from protocols/cyclic_peptide/crosslinker/thioether_util.hh
                bondangle1 = 1.960/(3.14159265358979323846264338327950288*180.0), # taken from numeric/NumericTraits.hh
                bondangle2 = 1.781/(3.14159265358979323846264338327950288*180.0),
                torsion = 180.,
                rsd1_before = 0,
                at1_before = "",
                rsd2_after = 0,
                at2_after = "",
                randomize_this_torsion = False,
                randomize_flanking_torsions = False,
            )
        else:
            GenKIC.close_bond(
                rsd1 = pose.total_residue(),
                at1 = "C",
                rsd2 = 1,
                at2 = "N",
                bondlength = 1.32,
                bondangle1 = 114.,
                bondangle2 = 123.,
                torsion = 180.,
                rsd1_before = 0,
                at1_before = "",
                rsd2_after = 0,
                at2_after = "",
                randomize_this_torsion = False,
                randomize_flanking_torsions = False,
            )
        
        # Apply backbone pertrubing effect
        self.debug: print("---Apply backbone perturbing while loop---")
        i = root+1
        while i != root:
            if i > cyclization_point_start: i = cyclization_point_end
            if i == root: 
                i+=1
                continue # This is just in cae the while doesnt close
            if self.thioether:
                if i == cyclization_point_start or i == cyclization_point_end: 
                    i+=1
                    continue

            GenKIC.add_perturber(
                    genkic.perturber.perturber_effect.randomize_backbone_by_rama_prepro
                    )
            GenKIC.add_residue_to_perturber_residue_list(i)
            i += 1

        # Add bump check filter
        GenKIC.add_filter(genkic.filter.filter_type.loop_bump_check)

        # Set rama check filters
        for i in range(1, pose.total_residue()+1):
            if i != pivot_res[0] and i != pivot_res[1] and i != pivot_res[2]: continue
#             if self.thioether:
#                 if i == cyclization_point_end or i == cyclization_point_start: continue
            GenKIC.add_filter(genkic.filter.filter_type.rama_prepro_check)
            GenKIC.set_filter_resnum(i)
            GenKIC.set_filter_rama_cutoff_energy(2)
            if i == pivot_res[0]: GenKIC.set_filter_attach_boinc_ghost_observer(True)

        self.debug: print("Pose before apply:", pose.sequence(), pose.total_residue())
        GenKIC.apply(pose)
        genkic_succ = GenKIC.last_run_successful()

        # Clear our genkic just incase
        GenKIC.clear_info()

        return genkic_succ

    def declare_terminal_bond_mover(self, pose: io.Pose, termini: sm.DeclareBond) -> int:
       """
       Fix terminal bond

       PARAMS
       ------
       :pose: Pose object
       :termini: DeclareBond mover
       """
       # Fix termini, though this isn't that important
       termini.set(
           res1 = pose.total_residue(),
           atom1 = 'C',
           res2 = 1,
           atom2 ='N',
           add_termini = True,
       )
       return termini.clone()

    def declare_terminal_bond(self, pose: io.Pose) -> int:
       """
       Fix terminal bond

       PARAMS
       ------
       :pose: Pose object
       """
       # Fix termini, though this isn't that important
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

    def foldtree_define(self, N: int) -> int:
        """
        Define what the root of our fold tree is

        PARAMS
        ------
        :N (int): Desired size of input structure

        RETURNS
        -------
        :root (int): The root residue index (zero indexed)
        """
        if N % 2 == 0:
                return int(N/2)
        else:
            return int((N-1)/2)

    def find_first_and_last_polymer_residues(self, pose:core.pose.Pose) -> tuple[int, int]:
        """Fine the first and last polymer residues in a pose

        RETURNS
        -------
        :first_polymer_res: First Polymer in pose
        :last_polymer_res: Last Polymer in pose
        """
        first_polymer_res, last_polymer_res = 0, 0
        for ir in range(1, pose.total_residue()+1):
            if pose.residue_type(ir).is_polymer():
                last_polymer_res = ir
                if first_polymer_res == 0: 
                    first_polymer_res = ir

        return first_polymer_res, last_polymer_res

    def find_first_and_last_thioether_lariat_residues(self, pose: core.pose.Pose) -> tuple[int,int]:
        """Given a pose, find the first and last thioether lariat bond-forming residues.

        @details First residue is 1 by definition (chloroacetyl goes where??); last is
        TYPICALLY C-term in the classic peptidream approach but does not have to be.
        n.b. Suga has methods for incorporating additional cysteines into bicyclic
        peptides, but for the moment "closest to the opposite terminus" is sufficient.
        """
        firstres = 1
        if self.lariat_sidechain_index != 0:
            lastres = self.lariat_sidechain_index
        else:
            lastres = self.find_last_disulf_res(pose)

        return firstres, lastres

    def set_up_terminal_thioether_lariat_variants(self, pose: core.pose.Pose) -> int:
        """Given a pose, add sidechain conjugation variant tytpes to the C-terminal
        cysteine and add a special acetyl terminus to the N-terminal residue

        RETURNS
        -------
        :cys_index: The thioether cyteien index
        """
        self.debug: print("-- Begin Pose Setup with Variant Types ---")
        firstres, lastres = self.find_first_and_last_thioether_lariat_residues(pose)
        self.debug: print("FirstRes, LastRes:", firstres, lastres)

        first_polymer, last_polymer = self.find_first_and_last_polymer_residues(pose)
        self.debug: print("FirstPolyer, LastPolymer:", firstres, lastres)

        core.pose.add_upper_terminus_type_to_pose_residue(pose, last_polymer)
        protocols.cyclic_peptide.crosslinker.set_up_thioether_variants(pose, firstres, lastres)
        self.debug: print("-- Finished Pose Setup with Variant Types ---")

        return lastres


    def add_cutpoints_upper_lower(self, pose: io.Pose) -> int:
        """
        Cyclize the pose using a Rosetta based method

        PARAMS
        ------
        :pose: Un cyclized pose
        """
        # Modify the variant types
        modifyvariant_nterm = sm.ModifyVariantTypeMover()
        modifyvariant_nterm.set_additional_type_to_add("CUTPOINT_UPPER")
        modifyvariant_nterm.set_residue_selector(select.ResidueIndexSelector(1))
        modifyvariant_nterm.apply(pose)

        modifyvariant_cterm = sm.ModifyVariantTypeMover()
        modifyvariant_cterm.set_additional_type_to_add("CUTPOINT_LOWER")
        modifyvariant_cterm.set_residue_selector(select.ResidueIndexSelector(pose.total_residue()))
        modifyvariant_cterm.apply(pose)

        pose.update_residue_neighbors()
        return 0

    def full_cyclize_native(self, pose: io.Pose) -> int:
        """
        Cyclize the native PDB pose using a Rosetta based method

        PARAMS
        ------
        :pose: Un cyclized pose
        """
        # Modify the variant types
        modifyvariant_nterm = sm.ModifyVariantTypeMover()
        modifyvariant_nterm.set_additional_type_to_remove("CUTPOINT_UPPER")
        modifyvariant_nterm.set_additional_type_to_remove("LOWER_TERMINUS_VARIANT")
        modifyvariant_nterm.set_additional_type_to_remove("ACETYLATED_NTERMINUS_VARIANT")
        modifyvariant_nterm.set_additional_type_to_remove("ACETYLATED_NTERMINUS_CONNECTION_VARIANT")
        modifyvariant_nterm.set_residue_selector(select.ResidueIndexSelector(1))
        modifyvariant_nterm.apply(pose)

        modifyvariant_cterm = sm.ModifyVariantTypeMover()
        modifyvariant_cterm.set_additional_type_to_remove("CUTPOINT_LOWER")
        modifyvariant_cterm.set_additional_type_to_remove("UPPER_TERMINUS_VARIANT")
        modifyvariant_cterm.set_additional_type_to_remove("DIMETHYLATED_CTERMINUS_VARIANT")
        modifyvariant_cterm.set_residue_selector(select.ResidueIndexSelector(pose.total_residue()))
        modifyvariant_cterm.apply(pose)

        # Fix termini, though this isn't that important
        declarebond = sm.DeclareBond()
        declarebond.set(
            res1 = pose.total_residue(),
            atom1 = 'C',
            res2 = 1,
            atom2 ='N',
            add_termini = True,
        )
        declarebond.apply(pose)
        return 0

    def generate_initial_empty_glycine_thioether_pose(self, N: int) -> core.pose.Pose:
        """
        Generate a single cyclic peptide pose of glycines and a cysteine at (cysteine_position) that is of size N

        PARAMS
        ------
        :N (int): Desired size of input structure

        RETURNS
        -------
        :thio_pose (core.pose.Pose): A Rosetta Pose instance of desired N-mer size
            FoldTree has been updated and is uses the middle residue as a root
        """
        
        thio_pose = core.pose.Pose()
        rts = core.chemical.ChemicalManager.get_instance().residue_type_set( 'fa_standard' )

        # Build up Glycines depending on where cysteine is placed
        # Remember Poses are 1-index not zero
        stubmover = protocols.cyclic_peptide.PeptideStubMover()
        stubmover.set_reset_mode(True)
        stubmover.reset_mover_data()
        for i in range(1, N):
            if i == 1:
                if self.thioether_chloroacetyl_homochiral:
                    stubmover.add_residue( "Append", "ALA", 0, False, "", 1, 0, None, "" )
                elif self.thioether_chloroacetyl_dchiral:
                    stubmover.add_residue( "Append", "DALA", 0, False, "", 1, 0, None, "" )
                else:
                    stubmover.add_residue( "Append", "GLY", 0, False, "", 1, 0, None, "" )
            else:
                if self.homochiral:
                    stubmover.add_residue( "Append", "ALA", 0, False, "", 1, 0, None, "" )
                else:
                    stubmover.add_residue( "Append", "GLY", 0, False, "", 1, 0, None, "" )
        stubmover.add_residue( "Append", "CYS", 0, False, "", 1, 0, None, "" )
        stubmover.apply(thio_pose)

        for i in range(1, thio_pose.total_residue()+1):
            thio_pose.set_omega(i, 180.0)

        # Set the fold tree
        root = self.foldtree_define(N)
        ft = core.kinematics.FoldTree()
        ft.clear()
        ft.add_edge(root, 1, -1)
        ft.add_edge(root, N, -1)
        ft.reorder(root)
        thio_pose.fold_tree(ft)

        core.pose.remove_variant_type_from_pose_residue(
                thio_pose,
                core.chemical.VariantType.UPPER_TERMINUS_VARIANT,
                thio_pose.total_residue(),
                )

        core.pose.remove_lower_terminus_type_from_pose_residue(
                thio_pose,
                thio_pose.total_residue(),
                )


        self.debug: print("--- Pose Made ---")
        self.debug: print("Seq:", thio_pose.sequence())
        self.debug: print("Seq:", thio_pose.fold_tree())

        return thio_pose

    def generate_initial_empty_glycine_pose(self, N: int) -> core.pose.Pose:
        """
        Generate a single cyclic peptide pose of glycines that is of size N

        PARAMS
        ------
        :N (int): Desired size of input structure

        RETURNS
        -------
        :gly_pose (core.pose.Pose): A Rosetta Pose instance that is all glycines of desired N-mer size
            FoldTree has been updated and is uses the middle residue as a root
            Cutpoints are added to the ends and a peptide bond is declared.
        """

        if self.empty_pose_test:
            empty_pose = core.pose.Pose()
            gly_pose = peptide_stub_mover(
                    pose = empty_pose,
                    residue_type = "GLY",
                    sequence_len = N,
                    )
            flip_trans_omegabonds(gly_pose, N)
        else:
            gly_pose = io.pose_from_sequence(
                seq="A"*N if self.homochiral else "G"*N,
                res_type="fa_standard",
                auto_termini=False,
            )

        # Set the fold tree
        root = self.foldtree_define(N)
        ft = core.kinematics.FoldTree()
        ft.clear()
        ft.add_edge(1, root, -1)
        ft.add_edge(root, N, -1)
        ft.reorder(root)
        gly_pose.fold_tree(ft)

        # Cyclize the pose (but don't declare a bond yet)
        self.add_cutpoints_upper_lower(gly_pose)
        return gly_pose

    def minimize_pose(
        self,
        pose: io.Pose,
        ) -> int:
        """
        Apply a minmover to the pose for scoring

        PARAMS
        ------
        :pose: Pose with desired motif residues and the scaffold is poly-glycine

        RETURNS
        -------
        :pose: minimized version of this pose
        """
        pose_before = pose.clone()

        # Set up a movemap
        mm = core.kinematics.MoveMap()
        mm.set_bb(True)
        mm.set_chi(True)

        # Set up minimover
        minmover = protocols.minimization_packing.MinMover()
        minmover.movemap(mm)
        minmover.score_function(self.scorefxn)
        minmover.min_type("lbfgs_armijo_nonmonotone")
        minmover.tolerance(1.0e-7)
        minmover.apply(pose)

        rmsd = core.scoring.bb_rmsd_including_O(pose_before, pose)
        pose.scores["MinMover_bbheavy_RMSD"] = rmsd
        self.scorefxn(pose)
        return 0

    def find_last_disulf_res(self, pose: core.pose.Pose) -> int:
        """Return the index of the last resiude that can form a disulfide

        PARAMS
        ------
        :pose: Our pose with the generated appropriate sequence

        RETURNS
        -------
        The index of our last residue that can form a disulfide
        """
        for i in range(pose.total_residue(), 0, -1):
            if pose.residue_type(i).is_sidechain_thiol() or pose.residue_type(i).is_disulfide_bonded():
                return i

    def declare_thioether_bond_mover(
            self, pose: core.pose, termini_close: protocols.simple_moves.DeclareBond,
            ) -> protocols.simple_moves.DeclareBond:

        firstres, lastres = self.find_first_and_last_thioether_lariat_residues(pose)

        # Get the name of the first sidechain connection atom:
        restype = pose.residue_type( firstres )
        ntermres_sc_connection_id = restype.n_possible_residue_connections()
        ntermres_sc_connection_atom_index = restype.residue_connect_atom_index( ntermres_sc_connection_id )
        ntermres_sc_connection_atom = restype.atom_name( ntermres_sc_connection_atom_index )

        # Get the name of the second sidechain connection atom:
        restype2 = pose.residue_type( lastres )
        sidechainres_sc_connection_id = restype2.n_possible_residue_connections()
        sidechainres_sc_connection_atom_index = restype2.residue_connect_atom_index( sidechainres_sc_connection_id )
        sidechainres_sc_connection_atom = restype2.atom_name( sidechainres_sc_connection_atom_index )

        self.debug: print("Setting up thioether lariat covalent bond from " + restype.base_name() + str(firstres) + ", atom " + ntermres_sc_connection_atom + " to residue " + restype2.base_name() + str(lastres) + ", atom " + sidechainres_sc_connection_atom + "." )

#         termini_close.set( sidechainres, sidechainres_sc_connection_atom, ntermres, ntermres_sc_connection_atom, false );
        termini_close.set( lastres, sidechainres_sc_connection_atom, firstres, ntermres_sc_connection_atom, False );
        return termini_close.clone()


    def declare_thioether_constraints(self, pose: core.pose.Pose) -> None:
        
        """Given a pose, find the first and last thioether lariat bond-forming residues.
        First residue is 1 by definition (chloroacetyl goes where??); last is
        TYPICALLY C-term in the classic peptidream approach but does not have to be.
        n.b. Suga has methods for incorporating additional cysteines into bicyclic
        peptides, but for the moment "closest to the opposite terminus" is sufficient.
        """
        n_index, c_index = self.find_first_and_last_thioether_lariat_residues(pose)
        protocols.cyclic_peptide.crosslinker.set_up_thioether_constraints(
            pose,
            n_index,
            c_index,
        )

    @timeit
    def generate_ensemble(self, size: int, nstruct: int, 
                          nofilter: bool,
                          ) -> int:
        """
        Helper function for running full process of pose gen + genkic/filter + minimize output
        for a given size 

        PARAMS
        ------
        :size: The desired output size of our pose
        :nstruct: Number of desired outputs
        :nofilter: Argparse argument for if you should filter or not based on hbonds

        RETURNS
        -------
        :ensemble: No actual return, but instead a silent file of starting poses for design
        """
        print("-"*4, "GenKIC, Size:", size, "Nstruct:", nstruct, "DEBUG:", self.DEBUG, "-"*4)

        if self.thioether:
            # Generate our initial thioether pose
            pose = self.generate_initial_empty_glycine_thioether_pose(size)
            self.thioether_sidechain_index = self.set_up_terminal_thioether_lariat_variants(pose)

            # Declare our thioether bond
            termini = protocols.simple_moves.DeclareBond()
            termini_updated = self.declare_thioether_bond_mover(pose, termini)
            termini_updated.apply(pose)

            # Correct virtual atoms if necessary
            self.debug: print("Thioether Sidechain Index:", self.thioether_sidechain_index)
            if self.thioether_sidechain_index != 0:
                protocols.cyclic_peptide.crosslinker.correct_thioether_virtuals(pose, 1, self.thioether_sidechain_index)
            # Apply thioether constraints
            self.declare_thioether_constraints(pose)
            pose.update_residue_neighbors()
        else:
            # Generate our initial pose
            pose = self.generate_initial_empty_glycine_pose(size)

            # Declare our terminal bond
            self.declare_terminal_bond(pose)


        # Setup SilentFile Output
        opts = silent.SilentFileOptions()
        opts.in_fullatom(True)
        opts.set_binary_output(True)
        silentFile = silent.SilentFileData(opts)
        if self.time_test:
            if self.thioether:
                out_name = f"genkicbb_thioether_size{s}_{datetime.date.today().strftime('%m%d%Y')}_timetest.silent"
            else:
                out_name = f"genkicbb_size{s}_{datetime.date.today().strftime('%m%d%Y')}_timetest.silent"
        elif self.thioether:
            out_name = f"genkicbb_thioether_size{s}_{datetime.date.today().strftime('%m%d%Y')}_timetest.silent"
        else:
            out_name = f"genkicbb_size{s}_{datetime.date.today().strftime('%m%d%Y')}.silent"

        # Apply genkic to our pose
        success = 0
        while success < nstruct:
        # for n in range(nstruct):
            # We clone our pose here, so that the bb torsion randomization is stochastic
            genkic_pose = pose.clone()
            self.randomize_mainchain(genkic_pose, not self.thioether)
            # Apply Genkic to pose and clone it
            genkic_succ = self.apply_genkic(genkic_pose, True)
#             print("Genkic success:", genkic_succ, "Number of Success:", success)
            if not genkic_succ:
                continue
            # Minimize with MinMover for small energetic improvement
            self.minimize_pose(genkic_pose)
            self.scorefxn(genkic_pose)
            # Calculate bb-bb hbond amount + Get True/False output
            hbond_filter = self.apply_hbond_filter(genkic_pose)
            # We want to filter by bb-bb hbonds (HbondFilter: True + use-hbond-filter: True)
            if hbond_filter and not nofilter:
                # Write to silent file
                silentStruct = silent.BinarySilentStruct(
                    opts, genkic_pose, f"size{size}_bb_{str(success+1).zfill(6)}",
                )
                if self.DEBUG:
                    genkic_pose.dump_pdb(f"debug_{size}_{str(success+1).zfill(6)}.pdb")
                success+=1
                if success % 1000 == 0:
                    print("Successfully Generated Filtered Poly-glycine Backbones:", success)
                if self.time_test:
                    silentFile.write_silent_struct(silentStruct, out_name)
                else:
                    silentFile.add_structure(silentStruct)
                # We don't want to filter by number of bb-bb hbonds. So Flag not set use-hbond-filter: False
            elif nofilter:
                # Write to silent file
                silentStruct = silent.BinarySilentStruct(
                    opts, genkic_pose, f"size{size}_bb_{str(success+1).zfill(6)}",
                )
                if self.DEBUG:
                    genkic_pose.dump_pdb(f"debug_{size}_{str(success+1).zfill(6)}.pdb")
                success+=1
                if success % 1000 == 0:
                    print("Successfully Generated Non-Filtered Poly-glycine Backbones:", success)
                if self.time_test:
                    silentFile.write_silent_struct(silentStruct, out_name)
                else:
                    silentFile.add_structure(silentStruct)

        if not self.time_test:
            silentFile.write_all(out_name)
        print("SilentFile Finished:", out_name)

        return 0

if __name__ == "__main__":
    p = argparse.ArgumentParser("Generate Macrocycle Backbones",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                )
    p.add_argument('-s', '--size', nargs="+", type=int, default=[6, 7, 8, 9], help="List of number of \
            residue macrocycles you'd like to generate")
    p.add_argument('-n', '--nstruct', type=int, default=10_000, help="Number of structures to generate per size.")
    p.add_argument('--debug', action="store_true", help="Dump PDBs for testing, unmute Rosetta, print helpful trace messages")
    p.add_argument('--nofilter', action="store_true", help="Dont apply strcit filter on hbonds")
    p.add_argument('--sample-root', action="store_true", help="If generating samples in-solution then this should be set \
            as your root residue is generally not an anchor.")
    p.add_argument("--time-test", action="store_true", help="This is only for a time comparison between XML and PyRosetta \
            as XML writes to disk for every structure, but here we do not. This can be done for time comparisons. Dont set \
            as it will make the process slower.")
    p.add_argument("--empty-pose-test", action="store_true", help="Run an empty pose test, since this is more representative of the XML script")
    p.add_argument("--thioether", action="store_true", help="Generate Thioether Cyclic Peptides")
    p.add_argument("--homochiral", action="store_true", help="This makes an all alanine peptide, as we want the peptide to be of the same chirality.")
    p.add_argument("--thioether-chloroacetyl-homochiral", action="store_true", help="This makes the choloracetyl residue be an alanine. \
            if not passed then it is a glycine.")
    p.add_argument("--thioether-chloroacetyl-dchiral", action="store_true", help="This makes the choloracetyl residue be an alanine. \
            if not passed then it is a glycine.")
    p.add_argument("--lariat-sidechain-index", type=int, default=0, help="Setup a special placement for the lariat n-c positions")
    args = p.parse_args()

    bbgen = BackboneGeneration(args.debug, args.sample_root, args.time_test, args.empty_pose_test,
                               args.thioether, args.homochiral, args.thioether_chloroacetyl_homochiral,
                               args.lariat_sidechain_index,
                               )
    for s in args.size:
        bbgen.generate_ensemble(s, args.nstruct, args.nofilter)

