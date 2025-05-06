#!/usr/bin/env python
#
# @author: Andrew Powers (apowers4@uoregon.edu / apowers@flatironinstitute.org)
# @brief: Modular filters that could be used
#

# Import Package
# from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
import pyrosetta.rosetta.core.select.residue_selector as residue_selector
import pyrosetta.rosetta.core.simple_metrics.metrics as metrics
import pyrosetta.rosetta.protocols.cyclic_peptide as cp
from pyrosetta.rosetta.core.scoring import ScoreFunction
import pyrosetta.rosetta.core as core

# Functions
def compare_rmsd_pose(
        refpose: core.pose.Pose,
        currpose: core.pose.Pose,
        refsel: residue_selector,
        cursel: residue_selector,
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
        cursel: residue_selector,
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

