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
def compare_rmsd_pose(
        refpose: core.pose.Pose,
        currpose: core.pose.Pose,
        refsel: ResidueSelector,
        cursel: ResidueSelector,
        filtername: str = "relax_after_design_",
        ) -> float:
    """Generate a RMSD metric with all bb heavy atoms including O

    PARAMS
    ------
    :refpose: The refernce pose that you want to compare against
    :currpose: The current pose you want to compare
    :refsel: the selection in the reference
    :curesel: the selection in the current pose
    :filtername: The name of our filter that will show up in the score file

    RETURNS
    -------
    Calculated RMSD Metric
    """
    # filter xml
    rmsdMetric = metrics.RMSDMetric()
    rmsdMetric.set_comparison_pose(refpose)
    rmsdMetric.set_residue_selector_reference(refsel)
    rmsdMetric.set_residue_selector(cursel)
    rmsdMetric.set_rmsd_type(
            core.scoring.rmsd_atoms.rmsd_protein_bb_heavy_including_O
                             )

    # apply the filter
    rmsdMetric.apply(currpose, filtername)
    return rmsdMetric.calculate(currpose)


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
    return oversat.apply(currpose)



def shapeComp(
        currpose: core.pose.Pose,
        jump_selection: int,
        ) -> float:
    """Calculate the shape complementarity between a jump

    PARAMS
    ------
    :currpose: The filled current pose
    :jump_selection: Which jump you want to check between

    RETURNS
    -------
    ShapeComp metric
    """
    # init our shape comp
    sc = protocols.simple_filters.ShapeComplementarityFilter()
    sc.set_jump_selector(jump_selection)
    sc.filtered_sc(0.5)
    sc.write_int_area(True)
    sc_score = sc.compute(currpose)
    return sc_score

def count_nonpolar_iteractions(
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
    interfaceHydrophobic.set_user_defined_name(filtername)
    # Extract count
    non_polar_count = interfaceHydrophobic.score(currpose)
    # return the bool, but also cash to our pose score
    return interfaceHydrophobic.apply(currpose), non_polar_count


def count_polar_interactions(
        currpose: core.pose.Pose,
        interfaceA: ResidueSelector,
        interfaceB: ResidueSelector,
        filtername_prefix: str = "polar_interactions_count_",
        ) -> int:
    """Count the number of non-covalent polar interactions across the two selections

    PARAMS
    ------
    :currpose: filled pose with target and receptor
    :interfaceA: The target/peptide interactions
    :interfaceB: The non target interface residues
    :filtername_prefix: Set the prefix name of the filter on the score header

    RETURNS
    -------
    Number of polar interactions 
    """
    # Setup our filter
    hbondMetric = core.simple_metrics.per_residue_metrics.HbondMetric()
    # Setup our selectors
    hbondMetric.set_include_self(False)
    hbondMetric.set_residue_selector(interfaceA)
    hbondMetric.set_residue_selector2(interfaceB)
    # Apply it, it will get stashed in our pose's score
    hbondMetric.apply(currpose, prefix=filtername_prefix)
    return 0

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
    resi_sel = selection.apply(currpose)

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
    scorePose = protocols.fold_from_loops.ScorePoseSementFromResidueSelectorFilter()
    scorePose.in_context(False)
    scorePose.scorefxn(scorefxn)
    scorePose.residue_selector(selection)
    scorePose.set_user_defined_name(filtername)
    # now apply to our pose
    score = scorePose.compute(currpose)
    # Add to pose score
    pose.scores[filtername] = score
    return score



