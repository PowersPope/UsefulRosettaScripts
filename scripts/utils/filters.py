#!/usr/bin/env python
#
# @author: Andrew Powers (apowers4@uoregon.edu / apowers@flatironinstitute.org)
# @brief: Modular filters that could be used
#

# Import Package
# from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
import pyrosetta.rosetta.core.select.residue_selector as residue_selector
import pyrosetta.rosetta.core.simple_metrics.metrics as metrics
import pyrosetta.rosetta.protocols as protocols
import pyrosetta.rosetta.protocols.cyclic_peptide as cp
import pyrosetta.rosetta.core as core

from pyrosetta.rosetta.core.scoring import ScoreFunction
from pyrosetta.rosetta.core.select.residue_selector import ResidueSelector

from typing import Tuple, Union

# Functions
def score_selection(
    refpose: core.pose.Pose,
    cursel: ResidueSelector,
    scorefxn: ScoreFunction,
    ) -> float:
    """Score the specific selection passed within our pose, and return that score value

    PARAMS
    ------
    :refpose: Our filled pose
    :cursel: The selection that is within our pose
    :scorefxn: the scorefunction that is used to generate our score

    RETURNS
    -------
    :score: The current score
    """
    # get the vector1 bool
    sel_vector1 = cursel.apply(refpose)
    # Now compute the score
    score = scorefxn.get_sub_score(refpose, sel_vector1)
    return score

def compare_rmsd_pose(
        refpose: core.pose.Pose,
        currpose: core.pose.Pose,
        refsel: ResidueSelector,
        cursel: ResidueSelector,
        superimpose: bool = False,
        calculate: bool = False,
        filtername: str = "relax_after_design_",
        ) -> float:
    """Generate a RMSD metric with all bb heavy atoms including O

    PARAMS
    ------
    :refpose: The refernce pose that you want to compare against
    :currpose: The current pose you want to compare
    :refsel: the selection in the reference
    :curesel: the selection in the current pose
    :superimpose: True run superimpose before calc, False dont
    :calculate: RMSD metric calc only and store in value, dont apply
    :filtername: The name of our filter that will show up in the score file

    RETURNS
    -------
    Calculated RMSD Metric
    """
    # Clone our input pose, so that it isnt altered
    currpose_clone = currpose.clone()

    # filter xml
    rmsdMetric = metrics.RMSDMetric()
    rmsdMetric.set_comparison_pose(refpose)
    rmsdMetric.set_residue_selector_reference(refsel)
    rmsdMetric.set_residue_selector(cursel)
    rmsdMetric.set_rmsd_type(
            core.scoring.rmsd_atoms.rmsd_protein_bb_heavy_including_O
                             )
    rmsdMetric.set_run_superimpose(superimpose)

    if calculate:
        return rmsdMetric.calculate(currpose_clone)
    else:
        # apply the filter
        rmsdMetric.apply(currpose_clone, filtername)
        return rmsdMetric.calculate(currpose_clone)


def determine_internal_bb_hbonds(
        currpose: core.pose.Pose,
        cursel: ResidueSelector,
        scorefxn: ScoreFunction,
        filtername: str = "design_bb_hbonds",
        ) -> int:
    """Count the number of bb internal 

    PARAMS
    ------
    :currpose: The current pose
    :cursel: The selection of the pose
    :scorefxn: The scorefunction used to determine hbonds
    :filtername: The name of our filter that will show up in the score file

    RETURNS
    -------
    Number of internal bb hbonds
    """
    # setup the bb internal hbonds
    bb_hbond_filter = cp.PeptideInternalHbondsFilter()
    # Only set bb hbonds
    bb_hbond_filter.set_hbond_types(
            backbone_backbone_setting = True,
            backbone_sidechain_setting = False,
            sidechain_sidechain_setting = False,
            )
    # Set the scorefxn
    bb_hbond_filter.set_scorefxn(scorefxn)
    # Set selection
    bb_hbond_filter.set_residue_selector(cursel)
    bb_hbond_filter.apply(currpose)

    currpose.scores[filtername] = bb_hbond_filter.report_sm(currpose)
    return bb_hbond_filter.report_sm(currpose)

def oversat_filter(
        currpose: core.pose.Pose,
        scorefxn: ScoreFunction,
        accp_selection: residue_selector,
        donor_selection: residue_selector,
        mainchain_only: bool,
        max_oversat_num: int = 0,
        ) -> bool:
    """Determine if there are any oversaturated hbond acceptors in our pose

    PARAMS
    ------
    :currpose: Our current pose filled object
    :max_oversat_num: The number of allowed oversat hbond acceptors

    RETURNS
    -------
    True if there is oversat, false if not
    """
    # init our filter
    oversat = cp.OversaturatedHbondAcceptorFilter()
    # Now apply appropriate settings
    oversat.set_acceptor_selector(accp_selection)
    oversat.set_consider_mainchain_only(mainchain_only)
    oversat.set_donor_selector(donor_selection)
    oversat.set_max_allowed_oversaturated(max_oversat_num)
    oversat.set_scorefxn(scorefxn)
    oversat.set_user_defined_name("oversat-check")
    currpose.scores["oversat-check"] = oversat.report_sm(currpose)
    return oversat.apply(currpose)



def shapeComp(
        currpose: core.pose.Pose,
        jump_selection: int,
        filtername: str = "shapeComp",
        ) -> float:
    """Calculate the shape complementarity between a jump

    PARAMS
    ------
    :currpose: The filled current pose
    :jump_selection: The jump you want to check between

    RETURNS
    -------
    ShapeComp metric
    """
    # Setup a jump selector
    jump_sel = core.select.jump_selector.JumpIndexSelector()
    jump_sel.jump(jump_selection)
    # init our shape comp
    sc = protocols.simple_filters.ShapeComplementarityFilter()
    sc.set_jump_selector(jump_sel)
    sc.filtered_sc(0.5)
    sc.write_int_area(True)
    sc_score = sc.score(currpose)
    currpose.scores[filtername] = sc_score
    return sc_score

def count_nonpolar_interactions(
        currpose: core.pose.Pose,
        interfaceA: ResidueSelector,
        interfaceB: ResidueSelector,
        scorefxn: ScoreFunction,
        apolar_res: str = "PHE,ILE,LEU,MET,PRO,THR,VAL,TRP,TYR,DPH,DIL,DLE,DME,DPR,DTH,DVA,DTR,DTY",
        threshold: int = 2,
        filtername: str = "interface_hydrophobic_filter",
        ) -> Tuple[Union[bool,float]]:
    """Count the number of interactions between two selections

    PARAMS
    ------
    :currpose: filled pose with target and receptor
    :interfaceA: The target/peptide interactions
    :interfaceB: The non target interface residues
    :scorefxn: An already defined scorefunction to be used
    :apolar_res: The apolar residues that you are going to check if they are making contacts. A list of comma separated name3
    :threshold: The number of residues that need to match our score_cut, to return true when applied.
    :filtername: Set the name of the filter on the score header

    RETURNS
    -------
    0: Bool of if it passes our threshold
    1: Number of non-polar interactions 
    """
    # init our interface hydrophobic filter
    interfaceHydrophobic = protocols.simple_filters.InterfaceHydrophobicResidueContactsFilter(
            hydrophobic_residue_contacts_threshold = threshold,
            target_selector = interfaceB,
            binder_selector = interfaceA,
            scorefxn = scorefxn,
            score_cut = -0.5,
            apolar_res = apolar_res,
            )
    # name our filter
#     interfaceHydrophobic.set_user_defined_name(filtername)
    # Extract count
    non_polar_count = interfaceHydrophobic.score(currpose)
    currpose.scores[filtername] = non_polar_count
    # return the bool, but also cash to our pose score
    non_polar_bool = interfaceHydrophobic.apply(currpose)
    return non_polar_bool, non_polar_count


def count_polar_interactions(
        currpose: core.pose.Pose,
        interfaceA: ResidueSelector,
        interfaceB: ResidueSelector,
        calculate_only: bool = False,
        filtername_prefix: str = "polar_interactions_count_",
        ) -> int:
    """Count the number of non-covalent polar interactions across the two selections

    PARAMS
    ------
    :currpose: filled pose with target and receptor
    :interfaceA: The target/peptide interactions
    :interfaceB: The non target interface residues
    :calculate_only: As we dont want this added to the pose, but we want the output
    :filtername_prefix: Set the prefix name of the filter on the score header

    RETURNS
    -------
    Number of polar interactions if calculate_only else 0 for successful call
    """
    # Setup our filter
    hbondMetric = core.simple_metrics.per_residue_metrics.HbondMetric()
    # Setup our selectors
    hbondMetric.set_include_self(False)
    hbondMetric.set_residue_selector(interfaceA)
    hbondMetric.set_residue_selector2(interfaceB)
    # Apply it, it will get stashed in our pose's score
    if calculate_only:
        polar_interactions = hbondMetric.calculate(currpose)
        return polar_interactions
    hbondMetric.apply(currpose, prefix=filtername_prefix)
    return 0

def count_hbondset_across_interface( pose: core.pose.Pose, 
                                    debug: bool = False,
                                    ) -> Tuple[int, dict[str,int]]:
    """Using HBondSet count the number of hbonds across the interface
    of a protein-peptide complex (can work for protein-protein complex too). 
    This function assumes the pose being passed only has two chains and one interface.

    PARAMS
    ------
    :pose: Pose that contains two chains and one interface. Written for a
        receptor and peptide pairing, but can be used for protein-protein
        complexes as well.
    :debug: Output interaction trace information

    RETURNS
    -------
    :num_interacts: The total number of hbonds across the interface
    :hbond_counts: A dictionary containing bb-bb, sc-sc, and bb-sc hbond
        counts across the interface.
    """
    # Setup our HBondSet and fill with relevent information from
    # our pose
    hbset = core.scoring.hbonds.HBondSet()
    core.scoring.hbonds.fill_hbond_set(pose, False, hbset)

    # Set up our type of interaction map
    hbond_counts = {"bb-bb": 0, "sc-sc": 0, "bb-sc": 0}

    # Loop over all hbonds within our filled HBondSet Object
    for i in range(1, hbset.nhbonds()+1):
        hb = hbset.hbond(i)
        donor_res = hb.don_res()
        acceptor_res = hb.acc_res()

        # Only select hbonds that are across the interface (i.e. different chains)
        if pose.residue(donor_res).chain() != pose.residue(acceptor_res).chain():
            assert pose.residue(donor_res).atom_is_polar_hydrogen(hb.don_hatm()), "The donor atom is not polar. The way HBondSet was filled might be wrong."
            # Determine if the donor or acceptor atom is a part of the backbone
            donor_bb = pose.residue(donor_res).atom_is_backbone(hb.don_hatm())
            acceptor_bb = pose.residue(acceptor_res).atom_is_backbone(hb.acc_atm())

            if debug: print("Donor Res:", donor_res, "Acceptor Res:", acceptor_res)

            # Map our type of hbond
            if donor_bb and acceptor_bb:
                hbond_counts["bb-bb"] += 1
            elif not donor_bb and not acceptor_bb:
                hbond_counts["sc-sc"] += 1
            else:
                hbond_counts["bb-sc"] += 1

    # Sum the total number of hbonds across the interface
    num_interacts = sum(hbond_counts.values())
    return num_interacts, hbond_counts


def grab_name3_sequence(
        currpose: core.pose.Pose,
        selection: ResidueSelector,
        filtername: str = "selection_sequence",
        ) -> str:
    """Grab the three letter code in a - separated string of the selection

    PARAMS
    ------
    :currpose: The filled posed
    :selection: The residue/chain subset we want the residue names of
    :filtername: The three letter sequence filtername in the pose score

    RETURNS
    -------
    :selection_residues_name: A string of joined three letter residue names
    """
    # Apply our selection to get the residue index
    resi_sel = selection.selection_positions(currpose)

    # init our return variable list 
    selection_residues = list()

    # now loop over and grab names
    for resi in resi_sel:
        selection_residues.append(currpose.residue(resi).name3())

    selection_residues_name = "-".join(selection_residues)
    currpose.scores[filtername] = selection_residues_name
    return selection_residues_name

def score_selection_outofcontext(
        currpose: core.pose.Pose,
        selection: ResidueSelector,
        scorefxn: ScoreFunction,
        filtername: str = "peptide_score_outofcontext",
        ) -> float:
    """Take our selection out of the context it is in (such as being bound) and
    score it alone

    PARAMS
    ------
    :currpose: Our filled multichain pose
    :selection: The chain/selection we would like to score by itself
    :scorefxn: Scorefunction, we want to use

    RETURNS
    -------
    :score: The out of context score of our selection
    """
    # init our filter
    scorePose = protocols.fold_from_loops.filters.ScorePoseSegmentFromResidueSelectorFilter()
    scorePose.in_context(False)
    scorePose.scorefxn(scorefxn)
    scorePose.residue_selector(selection)
    scorePose.set_user_defined_name(filtername)
    # now apply to our pose
    score = scorePose.compute(currpose)
    # Add to pose score
    currpose.scores[filtername] = score
    return score

def specificsite_hbond_filter(
        currpose: core.pose.Pose,
        res_selection: ResidueSelector,
        resnum: int,
        scorefxn: ScoreFunction,
        hbond_amount: int = 1,
        ) -> bool:
    """Check that a pose contains a hydrogen bond within our desired selection between
    chains.

    PARAMS
    ------
    :currpose: Our designed pose
    :res_selection: The residue selection that will be used to see if they contain
        a hbond with the residue specified in resnum
    :resnum: The rosetta residue index being checked to contain a hbond
    :scorefxn: Scorefxn to use for checking present Hbonds
    :hbond_amount: Number of hbonds present that we want

    RETURNS
    -------
    :contains_bond: Bool of if there is a desired hbond amount
    """
    # Setup our filter
    resHbondCounter = protocols.protein_interface_design.filters.HbondsToResidueFilter()
    resHbondCounter.set_backbone(True)
    resHbondCounter.set_bb_bb(True)
    resHbondCounter.set_energy_cutoff(-0.25)
    resHbondCounter.set_from_other_chains(True)
    resHbondCounter.set_from_same_chain(False)
    resHbondCounter.set_partners(hbond_amount)
    resHbondCounter.set_scorefxn(scorefxn)
    resHbondCounter.set_selector(res_selection)
    resHbondCounter.set_resnum(resnum)
    resHbondCounter.set_sidechain(True)
    contains_bond = resHbondCounter.apply(currpose)
    return contains_bond

def interface_analyzer(
        currpose: core.pose.Pose,
        scorefxn: ScoreFunction,
        pack_separated: bool = True,
        pack_before_separate: bool = True,
        pack_rounds: int = 1,
        packstat: bool = True,
        sasa_separated: bool = True,
        interface_energy: bool = True,
        interface_shapecomp: bool = True,
        dSASA: bool = True,
        hbond_sasaE: bool = True,
        dHbond_unsat: bool = True,
        jump_sel: int = 1,
        prefix_scorefile: str = "intAnalyzer_"
        ) -> int:
    """Perform a multitude of interface analyzer tasks that will be computed, if specified

    PARAMS
    ------
    :currpose: The current filled pose
    :selectionA: A selection within the pose
    :scorefxn: A predefined scorefunction
    :pack_separated: After separating repack the interface
    :pack_before_separate: Pack interface before separating
    :pack_rounds: Number of packing rounds
    :packstat: Run packstat
    :sasa_separated: Calculate the sasa of separated interfaces
    :interface_energy: Calculate the interface energy score
    :interface_shapecomp: Calculate the shapeComp
    :dSASA: Calculate the change in sasa
    :hbond_sasaE: Calculate the change in sasa energy for hbonds
    :dHbond_unsat: Calculate the number of unsat hbonds
    :jump_sel: The specific jump of our interface
    :prefix_scorefile: The prefix of the score headers in our scorefile

    RETURNS
    -------
    All filters specified will be applied to our poses scorefile.
    """
    # init our interfaceAnalyzer filter
    interfaceAnalyzer = protocols.analysis.InterfaceAnalyzerMover()
    # set our scorefunction
    interfaceAnalyzer.set_scorefunction(scorefxn)
    interfaceAnalyzer.set_pack_separated(pack_separated)
    interfaceAnalyzer.set_interface_jump(jump_sel)
    interfaceAnalyzer.set_calc_dSASA(dSASA)
    interfaceAnalyzer.set_calc_hbond_sasaE(hbond_sasaE)
    interfaceAnalyzer.set_compute_interface_delta_hbond_unsat(dHbond_unsat)
    interfaceAnalyzer.set_compute_interface_energy(interface_energy)
    interfaceAnalyzer.set_compute_interface_sc(interface_shapecomp)
    interfaceAnalyzer.set_compute_packstat(packstat)
    interfaceAnalyzer.set_compute_separated_sasa(sasa_separated)
    interfaceAnalyzer.set_pack_input(pack_before_separate)
    interfaceAnalyzer.set_pack_rounds(pack_rounds)
    interfaceAnalyzer.set_scorefile_reporting_prefix(prefix_scorefile)


    # Apply and add all score outputs to our pose
    interfaceAnalyzer.apply(currpose)
    interfaceAnalyzer.add_score_info_to_pose(currpose)
    return 0
