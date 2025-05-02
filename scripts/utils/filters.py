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
        ) -> float:
    """Generate a RMSD metric

    PARAMS
    ------
    :refpose: The refernce pose that you want to compare against
    :currpose: The current pose you want to compare
    :refsel: the selection in the reference
    :curesel: the selection in the current pose

    RETURNS
    -------
    Calculated RMSD Metric
    """
    # filter xml
    rmsdMetric = metrics.RMSDMetric()
    rmsdMetric.set_comparison_pose(refpose)
    rmsdMetric.set_residue_selector_reference(refsel)
    rmsdMetric.set_residue_selector(cursel)

    # apply the filter
    rmsdMetric.apply(currpose, "relax_after_design")
    return rmsdMetric.calculate(currpose)


def determine_internal_bb_hbonds(
        currpose: core.pose.Pose,
        cursel: residue_selector,
        scorefxn: ScoreFunction,
        ) -> int:
    """Count the number of bb internal 

    PARAMS
    ------
    :currpose: The current pose
    :cursel: The selection of the pose
    :scorefxn: The scorefunction used to determine hbonds

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

    return bb_hbond_filter.report_sm(currpose)


