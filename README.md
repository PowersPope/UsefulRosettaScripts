# UsefulRosettaScripts

A bunch of PyRosetta scripts that I have made for cyclic peptides over the course of my PhD.
Anyone is allowed to use them. 
Please let me know if you have any questions, though I will try my best to make the scripts self explanatory. 


## Repo Layout

- `scripts/` will house all of the `.py` or `.xml` scripts.
- `scripts/utils/` will house utility Rosetta functions and filters that can be used on their own.
    Additionally, they will be used in some scripts found in `scripts/` as well.
- `examples/` will house all of the simple example use cases for the scripts. (Though this will be added later.)


## Current Scripts
(This is old and needs to be updated, but I don't have the time at the moment! Yay for dissertation writing!)

- `calculate_ddg.py`: This can be used for a some calculating the ddG of binding on some dimeric complex. I 
it to handle my Cyclic Peptide designs, but I made that a flag and should work for general peptides and proteins as well.
General Workflow of the script: Score the complex, translate chainB away, repack the interface, rescore the system, calculate the ddg.
This can be performer `nstruct` number of times allowing you to get a mean ddG and a variance output.
Currently only calls in `ref2015_cart` scorefunction, but could be changed in the code. I will add in a flag option for score functions later.
```
python calculate_ddg.py -h
usage: calculate_ddg.py [-h] --path_to_pdb_dir PATH_TO_PDB_DIR [--relax_pdb] [--max_cpus MAX_CPUS] [--score_per_residue] [--cycpep] [--nstruct NSTRUCT]

optional arguments:
  -h, --help            show this help message and exit
  --path_to_pdb_dir PATH_TO_PDB_DIR
                        Path to directory that houses, bound target + peptide complexes
  --relax_pdb           Pass flag if you want to relax your pdb before you compute the ddg of binding
  --max_cpus MAX_CPUS   Set number of CPUs for Multiprocessing (default: os.cpu_count())
  --score_per_residue   Determines if you get per residue scores or just binding information
  --cycpep              Specify if working with cyclic peptides
  --nstruct NSTRUCT     Number of attempts to compute the binding ddG (default: 10)
```

- `insolution_genkic.py`: Example script on how to use genkic within PyRosetta. 
As I needed to make myself a script for comparison to models, and I noticed there was no good GenKIC PyRosetta example script on GitHub or in the PyRosetta workshop notebooks.
This uses `ref2015` weights by default:
```
usage: Generate Macrocycle Backbones [-h] [-s SIZE [SIZE ...]] [-n NSTRUCT] [--debug]
                                     [--nofilter] [--sample-root] [--time-test]

optional arguments:
  -h, --help            show this help message and exit
  -s SIZE [SIZE ...], --size SIZE [SIZE ...]
                        List of number of residue macrocycles you'd like to generate
                        (default: [6, 7, 8, 9])
  -n NSTRUCT, --nstruct NSTRUCT
                        Number of structures to generate per size. (default: 10000)
  --debug               Dump PDBs for testing, unmute Rosetta, print helpful trace messages
                        (default: False)
  --nofilter            Dont apply strcit filter on hbonds (default: False)
  --sample-root         If generating samples in-solution then this should be set as your
                        root residue is generally not an anchor. (default: False)
  --time-test           This is only for a time comparison between XML and PyRosetta as XML
                        writes to disk for every structure, but here we do not. This can be
                        done for time comparisons. Dont set as it will make the process
                        slower. (default: False)
```

## Author
Andrew Powers - Hosseinzadeh Lab, University of Oregon

Email: apowers4@uoregon.edu
