# Description

PepBD_PottsModel is a Fortran program that expresses peptide affinity for a receptor as a function of one- and two-body energies. Used to find the optimal amino acid sequence for a fixed system conformation. 

# Required input files
Details on files are provided below

- compiled executable
- input.txt
- sys.pdb
- the lib folder, placed at the location ~/PepBD/

# Installation and Compilation
Requirements: Intel ifort compiler (tested on versions 2017 - 2024)

1. Navigate to the `src` directory
2. Run the command `bash ../compile.sh`

# Running PepBD_PottsModel
Run from the terminal using `{PATH_TO_src}/main`. `input.txt` and `sys.pdb` must be in the same directory where the command is run.
The run time usually takes a few minutes. 

# Outputs of PepBD_PottsModel
1. `SingleEnergy.txt`: contains the one-body energy of each amino acid type at each location in the peptide.
2. `PairwiseEnergy.txt`: contains the two-body energy of all possible amino acid pairs at all possible locations in the peptide

# Details on required input files

## input.txt
An example is provided in `ExampleInputs` of the repository

- PDBFILE (required): the name of the pdb file containing the peptide-receptor complex. Usually `sys.pdb`.
- LIB_PATH (required): path to the lib folder. Usually `~/PepBD/`.
- PEP_RES (required): Number of residues in the peptide.
- RECEPTOR_NAME (required): Type of receptor. Options are nucleic, peptide, or molecule. Use molecule for designing plastic-binding peptides
- WEIGHTING_FACTOR (default=0.01): Factor that multiplies peptide-peptide interactions when calculating the PepBD score
  
## sys.pdb
Contains the peptide-receptor complex. An example is provided in `ExampleInputs` of the repository.

Since PepBD_Potts takes longer to run as the system size increases, receptor atoms that are 10 Angstroms from the peptide are typically removed, e.g. using VMD. sys.pdb should be formatted such that

- Only lines for atomic coordinates remain
- The peptide appears first, followed by the receptor
