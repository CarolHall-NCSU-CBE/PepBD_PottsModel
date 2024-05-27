# PepBD_PottsModel

REQUIRED FILES FOR RUNNING PROGRAM

The following files are required to run PepBD. Each is described in more detail below. All files should be placed in the directory for running PepBD, unless otherwise noted.

    all Fortran 90 files
    input.txt (details provided below)
    a pdb file of the peptide-receptor (details provided below)
    the uncompressed lib folder

REQUIRED PEPBD FILE - INPUT FILE

The name of this file must be "input.txt". An example input file is provided in the repository. The variables that can be specified in the input file are the following:

    PDBFILE (required): the name of the pdb file containing the peptide-receptor complex
    LIB_PATH (required): path to the lib folder 
    PEP_RES (required): Number of residues in the peptide.
    RECEPTOR_NAME (required): Type of receptor. Options are nucleic, peptide, or molecule
    SOLVATION_FLAG (default=0): Either 0 or 1; if 1, then will include generalizd Born (GB) solvation energy in one- and two-body energy calculations
    WEIGHTING_FACTOR (default=0.01): Factor that multiplies peptide-peptide interactions when calculating the PepBD score
    CG_RBORN (default=3.0): radius of bead placed at beta-carbon when calculating GB solvation energy (see paper for details on this)

REQUIRED PEPBD FILE - SYSTEM PDB FILE

Contains the peptide-receptor complex. Since PepBD runs slower as the system size increases, receptor atoms that are far from the peptide are removed from the pdb file. Common cutoff distances used are 8-10 Angstroms.

After reducing the system size, prepare the pdb file so it can be understood by PepBD by doing the following

    Remove any lines that do not have atomic coordinates, including TER and END
    If needed, reorder the atoms so the peptide appears first in the file, then is followed by the receptor
        If this is done, then fix the residue and atom numbering to match the new order. The tleap Amber of module is useful for this.

An example pdb file is provided in the repository.

PROGRAM COMPILATION

The intel compiler needs to be used (issues may occur if the program is compiled with gfortran). To compile, run the below lines:

ifort -c -O2 datatypes.f90
ifort -c -O2 sys_vars.f90
ifort -c -O2 input.f90
ifort -c -O2 math.f90
ifort -c -O2 utilities.f90
ifort -c -O2 pdbfile.f90
ifort -c -O2 database.f90 
ifort -c -O2 surface_area.f90
ifort -c -O2 energy_calculation.f90
ifort -c -O2 advanced_function.f90
ifort -c -O2 optimization_techniques.f90
ifort -c -O2 EvalEnergies.f90
ifort -O2 -o main *.o -mkl=sequential

GENERATED DATA FILES
The outputs from PepBD are the following: 

1. SingleEnergy.txt: ontains the one-body energy for each amino acid at each peptide residue. 
    -An example is provided in the repository
	
2. PairwiseEnergy.txt: contains the two-body energy between all possible amino acid pairs at all possible positions in the peptide.
	-An example is provided in the repository
    -Note that the energies in this file are scaled by WEIGHTING_FACTOR specified in input.txt
