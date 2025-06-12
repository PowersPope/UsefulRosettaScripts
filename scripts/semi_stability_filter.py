#!/usr/bin/env python
#
# @brief: External filter for assessing how well a peptide remains the lowest conformation
#           A mini-simple_cycpep_predict script is run to determine good ones
#

#### Import
# Base packages
import argparse
import time

# Import specific rosetta functions/classes
from pyrosetta import init
from pyrosetta.rosetta.core.import_pose.pose_stream import SilentFilePoseInputStream
# Import pyrosetta modules
import pyrosetta.io as io
import pyrosetta.rosetta.core as core
import pyrosetta.rosetta.core.io.silent as silent
import pyrosetta.rosetta.core.select.residue_selector as sel
import pyrosetta.rosetta.core.pack.task.operation as opts

# Import helper functions that I have created 
import utils.helper as helpfuncs
import utils.filters as filterfuncs

def main():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--infile", type=str, help="insilent file that was outputted by cyclicCAE for design")
    p.add_argument("--cartesian", action="store_true", help="Perform cartesian minimize/design (changes weights to ref2015_cart)")
    p.add_argument("--debug", action="store_true", help="DEBUG. We will not mute any rosetta output")
    p.add_argument("--peptide-chain", default=1, type=int, help="The particular chain of the cyclic peptide we want to test.")
    p.add_argument("--outsilentfile", type=str, default="out.silent", help="Specify the out silentfile name")
    p.add_argument("--cutoff", type=float, default=-4.0, help="Relax Score cutoff. It has to be at least this low or lower. If not then GenKIC isnt attempted.")
    p.add_argument("--tag", type=str, default=None, help="Tag to grab (not necessary)")
    args = p.parse_args()

    # Setup our scorefunction appropriately
    if args.cartesian and args.debug:
        init(extra_options="-in:file:fullatom true -out:file:silent_struct_type binary -ex1 -ex2 -score:weights ref2015_cart -detect_disulf false")
    elif args.cartesian:
        init(extra_options="-in:file:fullatom true -mute all -out:file:silent_struct_type binary -ex1 -ex2 -score:weights ref2015_cart -detect_disulf false")
    elif args.debug:
        init(extra_options="-in:file:fullatom true -out:file:silent_struct_type binary -ex1 -ex2 -score:weights ref2015 -detect_disulf false")
    else:
        init(extra_options="-in:file:fullatom true -mute all -out:file:silent_struct_type binary -ex1 -ex2 -score:weights ref2015 -detect_disulf false")

    # test
    t1 = time.time()

    # setup our out silentfile
    silentClass = helpfuncs.SilentFileWrite(args.outsilentfile)

    # Load in our silent file
    silopts = silent.SilentFileOptions()
    silopts.in_fullatom(True)
    silopts.set_binary_output(True)
    silentfile = silent.SilentFileData(silopts)
    if args.tag != None:
        silentfile._read_file(args.infile)
    else:
        silentfile = SilentFilePoseInputStream(args.infile)

    # init our scorefunction for design
    scorefxn = helpfuncs.init_scorefunction()
    print('Time %.3f' % float(time.time() - t1))

    # init an empty pose
    pose = core.pose.Pose()

    if args.tag != None:
        # Then tag here
        struct = silentfile.get_structure(args.tag)
        struct.fill_pose(pose)

        #setup variables
        pose_clone = pose.clone()
        peptide_sel = sel.ChainSelector(args.peptide_chain)

        # Grab out the peptide chain and only have that selection scored by rosetta
        selection_score, selection_pose = helpfuncs.generate_clean_conf(
            pose_clone,
            peptide_sel,
            scorefxn,
        )
        # Relax rthe structure
        helpfuncs.relax_selection(
                currpose = selection_pose, 
                res_selection = peptide_sel,
                scorefxn = scorefxn,
                cartesian = args.cartesian,
                )

        # Recompute the score for out of context relaxed peptide
        relax_score = scorefxn(selection_pose)
        print("Peptide Alone Score:", selection_score, "Relaxed Peptide Score:", relax_score)

        # lets test our conformation perturber
        selection_pose_clone = selection_pose.clone()
        is_stable = True
        for n in range(100):
            genkic_out = helpfuncs.apply_genkic(
                    pose = selection_pose_clone,
                    scorefxn = scorefxn,
                    randomize_root = True,
                    )
            helpfuncs.relax_selection(
                    currpose = genkic_out,
                    res_selection = sel.TrueResidueSelector(),
                    scorefxn = scorefxn,
                    cartesian = args.cartesian,
                    )
            bb_score = scorefxn(genkic_out)        

            # RMSD check 
            rmsd_small = filterfuncs.compare_rmsd_pose(
                    selection_pose, 
                    genkic_out,
                    peptide_sel,
                    peptide_sel,
                    True,
                    "test_%i_rmsd" % n
                    )
            print("RMSD Change %.4f" % rmsd_small)
            print("Relax Score: %.4f and GenKIC Score: %.4f" % (relax_score, bb_score))
            low_score_and_large_rmsd = rmsd_small >= 1.0 and bb_score <= (relax_score+2)
            if bb_score < relax_score or low_score_and_large_rmsd:
                print("Moving on found lower score during GenKIC...")
                print("Relax Score: %.4f and GenKIC Score: %.4f" % (relax_score, bb_score))
                is_stable = False
                break
        print("This peptide is stable? %s\n" % is_stable)
        if is_stable:
            # Write a pose to the file, but check to make sure it is available
            silentClass.write_when_not_busy(pose, core.pose.tag_from_pose(pose))

    else:
        struct_num = 0
        # Loop through our silentfile, until we have no more structs
        while silentfile.has_another_pose():
            struct_num += 1

            # Add in a pose
            silentfile.fill_pose(pose)

            #setup variables
            pose_clone = pose.clone()
            peptide_sel = sel.ChainSelector(args.peptide_chain)

            # Grab out the peptide chain and only have that selection scored by rosetta
            selection_score, selection_pose = helpfuncs.generate_clean_conf(
                pose_clone,
                peptide_sel,
                scorefxn,
            )
            # We want things that do not move dramatically when relaxed.
            pre_relax_pose = selection_pose.clone()

            # Relax the structure
            helpfuncs.relax_selection(
                    currpose = selection_pose, 
                    res_selection = peptide_sel,
                    scorefxn = scorefxn,
                    cartesian = args.cartesian,
                    )

            # Recompute the score for out of context relaxed peptide
            relax_score = scorefxn(selection_pose)

            # Compute RMSD
            rmsd_relax = filterfuncs.compare_rmsd_pose(
                    selection_pose, 
                    pre_relax_pose,
                    peptide_sel,
                    peptide_sel,
                    True,
                    "pre_relax_rmsd",
                    )

            print("Silent Struct:", struct_num, 
                  "Peptide Alone Score:", selection_score, 
                  "Relaxed Peptide Score:", relax_score,
                  "Relax RMSD BB Heavy Change:", rmsd_relax)
            if relax_score > args.cutoff or rmsd_relax > 0.7:
                print("Score is too high or RMSD change was too large! Moving on to next struct...")
                continue

            # lets test our conformation perturber
            selection_pose_clone = selection_pose.clone()
            is_stable = True
            time_start = time.time()
            for n in range(100):
                genkic_out = helpfuncs.apply_genkic(
                        pose = selection_pose_clone,
                        scorefxn = scorefxn,
                        randomize_root = True,
                        )
                helpfuncs.relax_selection(
                        currpose = genkic_out,
                        res_selection = sel.TrueResidueSelector(),
                        scorefxn = scorefxn,
                        cartesian = args.cartesian,
                        )
                bb_score = scorefxn(genkic_out)        

                # RMSD check 
                rmsd_small = filterfuncs.compare_rmsd_pose(
                        selection_pose, 
                        genkic_out,
                        peptide_sel,
                        peptide_sel,
                        True,
                        "test_%i_rmsd" % n
                        )
                print("RMSD Change %.4f for struct %i" % (rmsd_small, n))
                print("Relax Score: %.4f and GenKIC Score: %.4f" % (relax_score, bb_score))
                low_score_and_large_rmsd = rmsd_small >= 1.0 and bb_score <= (relax_score+2)
                if bb_score < relax_score:
                    print("Moving on found lower score during GenKIC...")
                    print("Relax Score: %.4f and GenKIC Score: %.4f" % (relax_score, bb_score))
                    is_stable = False
                    break
            print("This peptide is stable? %s\n" % is_stable)
            if is_stable:
                # Write a pose to the file, but check to make sure it is available
                silentClass.write_when_not_busy(pose, core.pose.tag_from_pose(pose))
                print("Total Time For GenKIC with nstruct of 100 was %.3f seconds" % float(time.time() - time_start))


if __name__ == "__main__":
    main()
